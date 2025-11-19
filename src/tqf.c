#include "tqf.h"

void fmpz_tqf_init(fmpz_tqf_t C)
{
    fmpz_factor_init(C->coeffs[0]);
    fmpz_factor_init(C->coeffs[1]);
    fmpz_factor_init(C->coeffs[2]);
}

void fmpz_tqf_set(fmpz_tqf_t C, fmpz_factor_t a, fmpz_factor_t b, fmpz_factor_t c)
{
    fmpz_factor_set(C->coeffs[0], a);
    fmpz_factor_set(C->coeffs[1], b);
    fmpz_factor_set(C->coeffs[2], c);
}

void fmpz_tqf_init_set(fmpz_tqf_t C, fmpz_factor_t a, fmpz_factor_t b, fmpz_factor_t c)
{
    fmpz_tqf_init(C);
    fmpz_tqf_set(C, a, b, c);
}

void fmpz_tqf_clear(fmpz_tqf_t C)
{
    fmpz_factor_clear(C->coeffs[0]);
    fmpz_factor_clear(C->coeffs[1]);
    fmpz_factor_clear(C->coeffs[2]);
}

void fmpz_tqf_print(const fmpz_tqf_t C)
{
    fmpz_t a, b, c;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    fmpz_factor_expand(a, C->coeffs[0]);
    fmpz_factor_expand(b, C->coeffs[1]);
    fmpz_factor_expand(c, C->coeffs[2]);

    flint_printf("(");
    fmpz_print(a);
    flint_printf(", ");
    fmpz_print(b);
    flint_printf(", ");
    fmpz_print(c);
    flint_printf(")");

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
}

static inline void reduce_square_append(fmpz_tqf_t A, fmpz_t* t, fmpz** p, ulong* e, slong l)
{
    fmpz_t tmp_mpz;
    
    fmpz_init(tmp_mpz);
    if (e[l] > 1) // embed square part in linear transformation
    {
        fmpz_pow_ui(tmp_mpz, p[l], e[l] / 2);
        fmpz_mul(t[(l + 1) % 3], t[(l + 1) % 3], tmp_mpz);
        fmpz_mul(t[(l + 2) % 3], t[(l + 2) % 3], tmp_mpz);
        e[l] %= 2; 
    }

    if (e[l] > 0)
    {
        _fmpz_factor_append(A->coeffs[l], p[l], e[l]);
    }

    fmpz_clear(tmp_mpz);
}

static inline int reduce_squares(fmpz_t* t, fmpz** p, ulong* e, slong null_idx)
{
    ulong nb_odd = 0;
    fmpz_t tmp_mpz;
    
    fmpz_init(tmp_mpz);
    for (slong l = 0; l < 3; l++)
    {
        if (l == null_idx) continue;

        if (e[l] > 1) // embed square part in linear transformation
        {
            fmpz_pow_ui(tmp_mpz, p[l], e[l] / 2);
            fmpz_mul(t[(l + 1) % 3], t[(l + 1) % 3], tmp_mpz);
            fmpz_mul(t[(l + 2) % 3], t[(l + 2) % 3], tmp_mpz);
            e[l] %= 2; 
        }
        nb_odd += e[l];
    }
    fmpz_clear(tmp_mpz);

    return nb_odd;
}

static inline void reduce_append_squares_medium_case(fmpz_tqf_t A, fmpz_t* t, fmpz** p, ulong* e, slong null_idx)
{
    slong l, non_null_idx;
    int nb_odd;
    
    non_null_idx = (null_idx + 1) % 3;  // we know tha both other ps are the same
    nb_odd = reduce_squares(t, p, e, null_idx);

    if (nb_odd == 2)
    {   
        #pragma GCC unroll 3
        for (l = 0; l < 3; l++)
        {
            if (l != null_idx) 
            {
                e[l]++;
            }
        }
        _fmpz_factor_append(A->coeffs[null_idx], p[non_null_idx], 1);
        reduce_squares(t, p, e, null_idx);
    }

    #pragma GCC unroll 3
    for (l = 0; l < 3; l++)
    {
        if (l != null_idx && e[l] != 0) 
        {
            _fmpz_factor_append(A->coeffs[l], p[non_null_idx], e[l]);
        }
    }
}

void fmpz_tqf_reduce(fmpz_tqf_t A, const fmpz_tqf_t C, fmpz_t t[3])
{
    slong l;
    slong done_idx, i1, i2; // first idx of coeff to be completely treated
    ulong tmp, nb_odd;
    slong idx[3]; // the current idx in the factorization struct of each coeff
    fmpz* p[3];   // the current primes   considered in the form
    ulong e[3], em;   // the current exponent considered in the form
    int   c[3];   // different comparison between curent lowest prime of the 3 coeffs of the form
    char  m[3] = {0};   // marks which index was treated during the main loop
    fmpz_t g;

    fmpz_init(g);

    A->coeffs[0]->sign = C->coeffs[0]->sign;
    A->coeffs[1]->sign = C->coeffs[1]->sign;
    A->coeffs[2]->sign = C->coeffs[2]->sign;

    A->coeffs[0]->num = 0;
    A->coeffs[1]->num = 0;
    A->coeffs[2]->num = 0;

    fmpz_set_ui(t[0], 1);
    fmpz_set_ui(t[1], 1);
    fmpz_set_ui(t[2], 1);

    idx[0] = 0;
    idx[1] = 0;
    idx[2] = 0;

    p[0] = C->coeffs[0]->p; e[0] = C->coeffs[0]->exp[0];
    p[1] = C->coeffs[1]->p; e[1] = C->coeffs[1]->exp[0];
    p[2] = C->coeffs[2]->p; e[2] = C->coeffs[2]->exp[0];

    c[0] = fmpz_cmp(p[0], p[1]);
    c[1] = fmpz_cmp(p[1], p[2]);
    c[2] = fmpz_cmp(p[2], p[0]);

    // main loop processes primes in ascending order
    while ((idx[0] < C->coeffs[0]->num) && 
           (idx[1] < C->coeffs[1]->num) && 
           (idx[2] < C->coeffs[2]->num))
    {
        if (m[0] == 1 || m[1] == 1) c[0] = fmpz_cmp(p[0], p[1]);
        if (m[1] == 1 || m[2] == 1) c[1] = fmpz_cmp(p[1], p[2]);
        if (m[2] == 1 || m[0] == 1) c[2] = fmpz_cmp(p[2], p[0]);

        // easy cases : p[i] alone and smallest
        // p[0] alone 
        if (c[0] < 0 && c[2] > 0) // c[0] is alone and smallest
        {
            reduce_square_append(A, t, p, e, 0);
            m[0] = 1; m[1] = 0; m[2] = 0;
        }

        // p[1] alone 
        else if (c[0] > 0 && c[1] < 0) // c[1] is alone and smallest
        {
            reduce_square_append(A, t, p, e, 1);
            m[0] = 0; m[1] = 1; m[2] = 0;
        }

        // p[2] alone 
        else if (c[1] > 0 && c[2] < 0) // c[0] is alone and smallest
        {
            reduce_square_append(A, t, p, e, 2);
            m[0] = 0; m[1] = 0; m[2] = 1;
        }

        // medium case : two p are the same 
        // p[0] == p[1]
        else if (c[0] == 0 && c[1] != 0)
        {
            reduce_append_squares_medium_case(A, t, p, e, 2);
            m[0] = 1; m[1] = 1; m[2] = 0;
        }

        // p[1] == p[2]
        else if (c[1] == 0 && c[2] != 0)
        {
            reduce_append_squares_medium_case(A, t, p, e, 0);
            m[0] = 0; m[1] = 1; m[2] = 1;
        }

        // p[2] == p[0]
        else if (c[2] == 0 && c[0] != 0)
        {
            reduce_append_squares_medium_case(A, t, p, e, 1);
            m[0] = 1; m[1] = 0; m[2] = 1;
        }
        // hard case
        // we have (c[0] == 0 && c[1] == 0)
        else // p[0] = p[1] = p[2] -> we have a common gcd to reduce first
        {
            tmp = 0;
            #pragma GCC unroll 3
            for (l = 1; l < 3; l++)
            {
                if (e[l] < e[tmp])
                {
                    tmp = l;
                }
            }

            em = e[tmp];
            e[0] -= em;
            e[1] -= em;
            e[2] -= em;

            // now e[tmp] = 0, lets reduce squares now
            nb_odd = reduce_squares(t, p, e, tmp);
            // if nb_odd = 2 (which is the max), then we still have two non-comprime
            // in this case we multiply everything by p and embed both squares
            if (nb_odd == 2)
            {
                e[0]++; e[1]++; e[2]++;
                reduce_squares(t, p, e, tmp);
            }

            // finally add everything to the reduced form
            #pragma GCC unroll 3
            for (l = 0; l < 3; l++)
            {
                if (e[l] != 0) 
                {
                    _fmpz_factor_append(A->coeffs[l], p[0], e[l]);
                }
            }

            m[0] = 1; m[1] = 1; m[2] = 1;
        }

        // final actualization
        #pragma GCC unroll 3
        for (l = 0; l < 3; l++)
        {
            if (m[l] == 1)
            {
                idx[l]++;
                p[l] = C->coeffs[l]->p + idx[l]; e[l] = C->coeffs[l]->exp[idx[l]];
            }
        }
    }   
    
    // now only 2 coefficients are left
    done_idx = 0;
    #pragma GCC unroll 3
    for (l = 0; l < 3; l++)
    {
        if (idx[l] == C->coeffs[l]->num)
        {
            done_idx = l;
            break;
        }
    }

    i1 = (done_idx + 1) % 3;
    i2 = (done_idx + 2) % 3;
    m[done_idx] = 0;

    // second loop processes primes in ascending order
    while ((idx[i1] < C->coeffs[i1]->num) && 
           (idx[i2] < C->coeffs[i2]->num))
    {
        c[i1] = fmpz_cmp(p[i1], p[i2]); // this is technically an extra comparison more than necessary...

        // easy cases : p[i] alone and smallest
        // p[i1] alone 
        if (c[i1] < 0) // c[i1] is alone and smallest
        {
            reduce_square_append(A, t, p, e, i1);
            m[i1] = 1; m[i2] = 0;
        }

        else if (c[i1] > 0) // c[i2] is alone and smallest
        {
            reduce_square_append(A, t, p, e, i2);
            m[i2] = 1; m[i1] = 0;
        }

        // medium case : two p are the same 
        // p[0] == p[1]
        else // c[i1] == 0
        {
            reduce_append_squares_medium_case(A, t, p, e, done_idx);
            m[i1] = 1; m[i2] = 1;
        }

        // final actualization
        #pragma GCC unroll 3
        for (l = 0; l < 3; l++)
        {
            if (m[l] == 1)
            {
                idx[l]++;
                p[l] = C->coeffs[l]->p + idx[l]; e[l] = C->coeffs[l]->exp[idx[l]];
            }
        }
    }

    // now only 1 coefficients is left
    #pragma GCC unroll 3
    for (l = 0; l < 3; l++)
    {
        if (idx[l] != C->coeffs[l]->num)
        {
            while (idx[l] < C->coeffs[l]->num)
            {
                reduce_square_append(A, t, p, e, l);
                idx[l]++;
                p[l] = C->coeffs[l]->p + idx[l]; 
                e[l] = C->coeffs[l]->exp[idx[l]];
            }
            break;
        }
    }

    // finally reduce translation
    fmpz_gcd3(g, t[0], t[1], t[2]);
    fmpz_divexact(t[0], t[0], g);
    fmpz_divexact(t[1], t[1], g);
    fmpz_divexact(t[2], t[2], g);

    fmpz_clear(g);
}

int fmpz_tqf_certif(fmpz_t k, const fmpz_t a, const fmpz_t b, const fmpz_factor_t c)
{
    slong i;
    fmpz* roots;
    fmpz_t m_ab, tmp;
    int err;

    fmpz_init(tmp);
    fmpz_init(m_ab);
    
    err = 1;
    roots = _fmpz_vec_init(c->num);
    fmpz_mul(m_ab, a, b);
    fmpz_neg(m_ab, m_ab);

    #pragma GCC unroll 3
    for (i = 0; i < c->num; i++)
    {
        fmpz_mod(tmp, m_ab, c->p + i);
        err = fmpz_sqrtmod(roots + i, tmp, c->p + i);

        if (err != 1) 
        {
            goto CLEAN;
        }

        // find smallest root
        fmpz_sub(tmp, c->p + i, roots + i);
        if (fmpz_cmp(tmp, roots + i) < 0)
        {
            fmpz_set(roots + i, tmp);
        }
    }
    err = fmpz_multi_CRT(k, c->p, roots, c->num, 0);

    CLEAN:;
    fmpz_clear(tmp);
    fmpz_clear(m_ab);
    _fmpz_vec_clear(roots, c->num);

    return err;
}

// this is some ugly bit tricks to compute x mod 4 with x being an fmpz
// depending on the storage type one has to fix some bits
#define FMPZ_LIMB_MOD4(f, l) (fmpz_sgn(f) == -1 ? ((4 - l) & 3) : (l))
#define FMPZ_MOD4(f) (COEFF_IS_MPZ(*(f)) ? FMPZ_LIMB_MOD4(f, (FMPZ_TO_ZZ(*(f))->ptr[0] & 3)) : ((*(f)) & WORD(3)))
#define FMOD4(a4, b4, c4, x4, y4, z4)  (((a4)*(x4)*(x4) + (b4)*(y4)*(y4) + (c4)*(z4)*(z4)) & 3)
#define INDEX2(a4, b4, c4, abc4, x4, y4, z4) (((FMOD4(a4, b4, c4, x4, y4, z4)) >> (((abc4) & 1) ^ 1)) & 1)

static inline ulong _fmpz_tqf_vec3_index2(const fmpz* v, ulong a4, ulong b4, ulong c4, ulong abc4)
{
    ulong x4, y4, z4;

    x4 = FMPZ_MOD4(v + 0);
    y4 = FMPZ_MOD4(v + 1);
    z4 = FMPZ_MOD4(v + 2);

    return INDEX2(a4, b4, c4, abc4, x4, y4, z4);
}

int _fmpz_tqf_test_sol(const fmpz* a, const fmpz* b, const fmpz* c, const fmpz* S)
{
    fmpz_t x2, y2, z2;
    int ret;

    fmpz_init(x2);
    fmpz_init(y2);
    fmpz_init(z2);

    fmpz_mul(x2, S + 0, S + 0);
    fmpz_mul(y2, S + 1, S + 1);
    fmpz_mul(z2, S + 2, S + 2);

    fmpz_mul(x2, x2, a);
    fmpz_mul(y2, y2, b);
    fmpz_mul(z2, z2, c);

    fmpz_add(x2, x2, y2);
    fmpz_add(x2, x2, z2);

    ret = fmpz_is_zero(x2);

    fmpz_clear(x2);
    fmpz_clear(y2);
    fmpz_clear(z2);

    return ret;
}

int fmpz_tqf_solve_reduced(fmpz_tqf_t R, fmpz_t z1, fmpz_t z2, fmpz_t z3)
{
    fmpz_t a, b, c;
    fmpz_t k1, k2, k3;
    int err;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(k1);
    fmpz_init(k2);
    fmpz_init(k3);
    err = 1;

    // compute the expansions
    fmpz_factor_expand(a, R->coeffs[0]);
    fmpz_factor_expand(b, R->coeffs[1]);
    fmpz_factor_expand(c, R->coeffs[2]);
    
    // take care of special cases
    if ((fmpz_cmp_si(a,  1) == 0 && fmpz_cmp_si(b, -1) == 0) || 
        (fmpz_cmp_si(a, -1) == 0 && fmpz_cmp_si(b,  1) == 0))
    {
        fmpz_set_si(z1, 1);
        fmpz_set_si(z2, 1);
        fmpz_set_si(z3, 0);
        return err;
    }

    if ((fmpz_cmp_si(a,  1) == 0 && fmpz_cmp_si(c, -1) == 0) || 
             (fmpz_cmp_si(a, -1) == 0 && fmpz_cmp_si(c,  1) == 0))
    {
        fmpz_set_si(z1, 1);
        fmpz_set_si(z2, 0);
        fmpz_set_si(z3, 1);
        return err;
    }
    
    if ((fmpz_cmp_si(b,  1) == 0 && fmpz_cmp_si(c, -1) == 0) || 
             (fmpz_cmp_si(b, -1) == 0 && fmpz_cmp_si(c,  1) == 0))
    {
        fmpz_set_si(z1, 0);
        fmpz_set_si(z2, 1);
        fmpz_set_si(z3, 1);
        return err;
    }

    // compute solubility certificate
    err = fmpz_tqf_certif(k3, a, b, R->coeffs[2]);
    if (err != 1)
    {
        goto CLEAN;
    }  

    err = fmpz_tqf_certif(k2, c, a, R->coeffs[1]);
    if (err != 1)
    {
        goto CLEAN;
    }  

    err = fmpz_tqf_certif(k1, b, c, R->coeffs[0]);
    if (err != 1)
    {
        goto CLEAN;
    }  

    // past this point we expect a solution and no error should happen
    fmpz_t tmp, u, v, ap, bp;
    fmpz S[3];
    fmpz_t alpha, beta, gamma, bc;
    fmpz_mat_t M;
    fmpz_lll_t LLL;
    fmpz* curr_row;
    ulong a4, b4, c4, abc4;
    slong i, j, k;

    fmpz_mat_init(M, 3, 3);
    fmpz_lll_context_init(LLL, 0.99, 0.501, Z_BASIS, APPROX);

    fmpz_init(S + 0);
    fmpz_init(S + 1);
    fmpz_init(S + 2);

    fmpz_init(tmp);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(ap);
    fmpz_init(bp);

    fmpz_init(alpha);
    fmpz_init(beta);
    fmpz_init(gamma);
    fmpz_init(bc);

    fmpz_mul(bc, b, c);

    fmpz_xgcd(tmp, u, v, b, c);
    fmpz_xgcd(tmp, ap, bp, a, bc);


    // alpha = (bp*c*k1) % a
    fmpz_mul(alpha, bp, c);
    fmpz_mod(alpha, alpha, a);
    fmpz_mul(alpha, alpha, k1);
    fmpz_mod(alpha, alpha, a);

    // beta = (u*ap*b*k3) % (b*c)
    fmpz_mul(beta, u, ap);
    fmpz_mod(beta, beta, bc);
    fmpz_mul(beta, beta, b);
    fmpz_mod(beta, beta, bc);
    fmpz_mul(beta, beta, k3);
    fmpz_mod(beta, beta, bc);

    // gamma = (v*ap*c*k2) % (b*c)
    fmpz_mul(gamma, v, ap); // gamma is wrong !!!!!!!!!!!!!!!!!!
    fmpz_mod(gamma, gamma, bc);
    fmpz_mul(gamma, gamma, c);
    fmpz_mod(gamma, gamma, bc);
    fmpz_mul(gamma, gamma, k2);
    fmpz_mod(gamma, gamma, bc);

    // cleanup ks and other temps
    fmpz_clear(k1);
    fmpz_clear(k2);
    fmpz_clear(k3);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(ap);
    fmpz_clear(bp);

    // v1 = (b*c, 0, 0)
    fmpz_set   (fmpz_mat_entry(M, 0, 0), bc);
    fmpz_set_ui(fmpz_mat_entry(M, 0, 1), 0);
    fmpz_set_ui(fmpz_mat_entry(M, 0, 2), 0);

    // v2 = (a*beta, a, 0)
    fmpz_mul   (fmpz_mat_entry(M, 1, 0), a, beta);
    fmpz_set   (fmpz_mat_entry(M, 1, 1), a);
    fmpz_set_ui(fmpz_mat_entry(M, 1, 2), 0);

    // v3 = (alpha*beta + gamma, alpha, 1)
    fmpz_mul(tmp, alpha, beta);
    fmpz_add   (fmpz_mat_entry(M, 2, 0), tmp, gamma);
    fmpz_set   (fmpz_mat_entry(M, 2, 1), alpha);
    fmpz_set_ui(fmpz_mat_entry(M, 2, 2), 1);

    // cleanup rest of temps
    fmpz_clear(alpha);
    fmpz_clear(beta);
    fmpz_clear(gamma);
    fmpz_clear(tmp);
    fmpz_clear(bc);

    // compute reduced coeffs
    a4 = FMPZ_MOD4(a);
    b4 = FMPZ_MOD4(b);
    c4 = FMPZ_MOD4(c);
    abc4 = (a4*b4*c4) & 3;

    // compute odd index2 vector
    #pragma GCC unroll 3
    for (j = 0; j < 3; j++)
    {
        if (_fmpz_tqf_vec3_index2(fmpz_mat_row(M, j), a4, b4, c4, abc4) == 1)
        {
            break;
        }
    }

    // compute index2 lattice on place
    #pragma GCC unroll 3
    for (i = 0; i < 3; i++)
    {
        if (i != j)
        {
            curr_row = fmpz_mat_row(M, i);
            if (_fmpz_tqf_vec3_index2(curr_row, a4, b4, c4, abc4) == 1)
            {
                _fmpz_vec_sub(curr_row, curr_row, fmpz_mat_row(M, j), 3);
            }
        }
    }
    curr_row = fmpz_mat_row(M, j);
    _fmpz_vec_add(curr_row, curr_row, curr_row, 3);

    // compute scaling factors
    fmpz_abs(S + 0, a);
    fmpz_abs(S + 1, b);
    fmpz_abs(S + 2, c);
    fmpz_sqrt(S + 0, S + 0);
    fmpz_sqrt(S + 1, S + 1);
    fmpz_sqrt(S + 2, S + 2);

    // scale columns
    #pragma GCC unroll 3
    for (i = 0; i < 3; i++)
    {
        #pragma GCC unroll 3
        for (j = 0; j < 3; j++)
        {
            curr_row = fmpz_mat_entry(M, i, j);
            fmpz_mul(curr_row, curr_row, S + j);
        }
    }

    // do the actual reduction
    fmpz_lll_mpf(M, NULL, LLL);

    // unscale columns
    #pragma GCC unroll 3
    for (i = 0; i < 3; i++)
    {
        #pragma GCC unroll 3
        for (j = 0; j < 3; j++)
        {
            curr_row = fmpz_mat_entry(M, i, j);
            fmpz_divexact(curr_row, curr_row, S + j);
        }
    }

    // iterate over small vectors and test for solution
    #pragma GCC unroll 3
    for (i = -1; i < 2; i++)
    {
        #pragma GCC unroll 3
        for (j = -1; j < 2; j++)
        {
            #pragma GCC unroll 3
            for (k = -1; k < 2; k++)
            {
                if (i == 0 && j == 0 && k == 0)
                {
                    continue;
                }

                // we reuse the scaling vector for tmp vector
                if (i > 0)
                {
                    _fmpz_vec_set(S, fmpz_mat_row(M, 0), 3);
                }
                else if (i < 0)
                {
                    _fmpz_vec_neg(S, fmpz_mat_row(M, 0), 3);
                }
                else
                {
                    fmpz_zero(S + 0);
                    fmpz_zero(S + 1);
                    fmpz_zero(S + 2);
                }

                if (j > 0)
                {
                    _fmpz_vec_add(S, S, fmpz_mat_row(M, 1), 3);
                }
                else if (j < 0)
                {
                    _fmpz_vec_sub(S, S, fmpz_mat_row(M, 1), 3);
                }

                if (k > 0)
                {
                    _fmpz_vec_add(S, S, fmpz_mat_row(M, 2), 3);
                }
                else if (k < 0)
                {
                    _fmpz_vec_sub(S, S, fmpz_mat_row(M, 2), 3);
                }

                if (_fmpz_tqf_test_sol(a, b, c, S))
                {
                    // just reuse a to store gcd
                    fmpz_gcd3(a, S + 0, S + 1, S + 2);

                    fmpz_divexact(z1, S + 0, a);
                    fmpz_divexact(z2, S + 1, a);
                    fmpz_divexact(z3, S + 2, a);

                    goto SUCCESS;
                }
            }
        }
    }

    flint_printf("tqf failed to find a solution despite solubility certif !\n");
    fmpz_tqf_print(R); flint_printf("\n");
    
    flint_abort();

    SUCCESS:;
    
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);

    fmpz_mat_clear(M);

    fmpz_clear(S + 0);
    fmpz_clear(S + 1);
    fmpz_clear(S + 2);

    return 1;

    CLEAN:;
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(k1);
    fmpz_clear(k2);
    fmpz_clear(k3);

    return err;
}

void fmpz_tqf_parametrize(qfb_t q1, qfb_t q2, qfb_t q3, const fmpz_tqf_t R, 
    const fmpz_t x0, const fmpz_t y0, const fmpz_t z0)
{
    fmpz_t a, b, c, ax0, e, z02, by0, tmp, ay0, eax0, eay0;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(ax0);
    fmpz_init(e);
    fmpz_init(z02);
    fmpz_init(by0);
    fmpz_init(tmp);
    fmpz_init(ay0);
    fmpz_init(eax0);
    fmpz_init(eay0);

    fmpz_factor_expand(a, R->coeffs[0]);
    fmpz_factor_expand(b, R->coeffs[1]);
    fmpz_factor_expand(c, R->coeffs[2]);
    fmpz_mul(z02, z0, z0);
    
    fmpz_invmod(tmp, y0, z02);
    fmpz_mul(e, x0, tmp);
    fmpz_mod(e, e, z02);

    fmpz_mul(ax0, a, x0);
    fmpz_mul(ay0, a, y0);
    fmpz_mul(by0, b, y0);
    fmpz_mul(eax0, e, ax0);
    fmpz_mul(eay0, e, ay0);

    // q1 = a*x0*U**2 + 
    //      (2 * (e*a*x0 + b*y0) // z0) * U*V + 
    //      ((a*x0*e**2 + 2*e*b*y0 - b*x0) // z0**2) * V**2

    fmpz_set(q1->a, ax0);

    fmpz_add(q1->b, eax0, by0);
    fmpz_add(q1->b, q1->b, q1->b);
    fmpz_divexact(q1->b, q1->b, z0);
    
    fmpz_mul(q1->c, eax0, e);
    fmpz_mul(tmp, e, by0);
    fmpz_add(tmp, tmp, tmp);
    fmpz_add(q1->c, q1->c, tmp);
    fmpz_mul(tmp, b, x0);
    fmpz_sub(q1->c, q1->c, tmp);
    fmpz_divexact(q1->c, q1->c, z02);

    // q2 = -a*y0*U**2 + 
    //      (2 * (a*x0 - a*y0*e) // z0) * U*V + 
    //      ((2*a*x0*e - a*y0*e**2 + b*y0) // z0**2) * V**2

    fmpz_neg(q2->a, ay0);

    fmpz_sub(q2->b, ax0, eay0);
    fmpz_add(q2->b, q2->b, q2->b);
    fmpz_divexact(q2->b, q2->b, z0);

    fmpz_add(q2->c, eax0, eax0);
    fmpz_mul(tmp, eay0, e);
    fmpz_sub(q2->c, q2->c, tmp);
    fmpz_add(q2->c, q2->c, by0);
    fmpz_divexact(q2->c, q2->c, z02);

    // q3 = a*z0*U**2 + 
    //      2*a*e*U*V + 
    //      ((a*e**2 + b) // z0) * V**2

    fmpz_mul(q3->a, a, z0);

    fmpz_mul(tmp, a, e);
    fmpz_add(q3->b, tmp, tmp);

    fmpz_mul(q3->c, tmp, e); // tmp still has a*e
    fmpz_add(q3->c, q3->c, b);
    fmpz_divexact(q3->c, q3->c, z0);

    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(ax0);
    fmpz_clear(e);
    fmpz_clear(z02);
    fmpz_clear(by0);
    fmpz_clear(tmp);
    fmpz_clear(ay0);
    fmpz_clear(eax0);
    fmpz_clear(eay0);
}

void test_tqf()
{
    fmpz c[3];
    fmpz sol[3];
    fmpz_factor_t af, bf, cf;

    fmpz_tqf_t TBF, TBF_r;
    fmpz_t tr[3];

    int err;
    ulong l, nb;

    flint_rand_t rstate;

    flint_rand_init(rstate);

    fmpz_init(c + 0);
    fmpz_init(c + 1);
    fmpz_init(c + 2);

    fmpz_init(sol + 0);
    fmpz_init(sol + 1);
    fmpz_init(sol + 2);

    fmpz_init(tr[0]); 
    fmpz_init(tr[1]); 
    fmpz_init(tr[2]);

    // init qbf
    fmpz_tqf_init(TBF);
    fmpz_tqf_init(TBF_r);

    // init factors
    fmpz_factor_init(af);
    fmpz_factor_init(bf);
    fmpz_factor_init(cf);

    flint_printf("Testing from bitlength %ld to %ld with %ld run per sample\n", 
        TQL_TEST_MIN_BIT_LEN, TQL_TEST_MAX_BIT_LEN, TQL_TEST_NB_PER_SAMPLE);

    for (l = TQL_TEST_MIN_BIT_LEN; l <= TQL_TEST_MAX_BIT_LEN; l++)
    {
        flint_printf("Test with bitlength %ld\n", l);

        nb = 0;
        while (nb < TQL_TEST_NB_PER_SAMPLE)
        {
            // parameters for c_0 x^2 + c_1 y^2 + c_2 z^2 = 0
            do 
            {
                fmpz_randbits(c + 0, rstate, l);
                fmpz_randbits(c + 1, rstate, l);
                fmpz_randbits(c + 2, rstate, l);

                fmpz_mul_si(c + 0, c + 0, (n_randlimb(rstate) % 2) * 2 - 1);
                fmpz_mul_si(c + 1, c + 1, (n_randlimb(rstate) % 2) * 2 - 1);
                fmpz_mul_si(c + 2, c + 2, (n_randlimb(rstate) % 2) * 2 - 1);

            } while (fmpz_is_zero(c + 0) || fmpz_is_zero(c + 1) || fmpz_is_zero(c + 2) ||
                     ((fmpz_sgn(c + 0) == fmpz_sgn(c + 1)) && (fmpz_sgn(c + 1) == fmpz_sgn(c + 2))));
            

            // factor coeffs
            fmpz_factor(af, c + 0);
            fmpz_factor(bf, c + 1);
            fmpz_factor(cf, c + 2);

            // set tqf
            fmpz_tqf_set(TBF, af, bf, cf);

            // reduce qbf
            fmpz_tqf_reduce(TBF_r, TBF, tr);
            
            // solving qbf
            err = fmpz_tqf_solve_reduced(TBF_r, sol + 0, sol + 1, sol + 2);
            
            // translating
            fmpz_mul(sol + 0, sol + 0, tr[0]);
            fmpz_mul(sol + 1, sol + 1, tr[1]);
            fmpz_mul(sol + 2, sol + 2, tr[2]);

            // print solution
            if (err == 1)
            {   
                flint_printf("\nTBF : "); fmpz_tqf_print(TBF); flint_printf("\n");
                flint_printf("TBF reduced : "); fmpz_tqf_print(TBF_r); flint_printf("\n");
                flint_printf("tr : "); 
                fmpz_print(tr[0]); flint_printf(" "); 
                fmpz_print(tr[1]); flint_printf(" ");
                fmpz_print(tr[2]); flint_printf("\n");

                flint_printf("Found solution : ("); 
                fmpz_print(sol + 0); flint_printf(", ");
                fmpz_print(sol + 1); flint_printf(", ");
                fmpz_print(sol + 2); flint_printf(")\n");

                if (!_fmpz_tqf_test_sol(c + 0, c + 1, c + 2, sol))
                {
                    flint_printf("Solution FAIL !!\n");
                    flint_abort();
                }
                nb++;
            }

        }
    }

    flint_rand_clear(rstate);

    fmpz_clear(c + 0);
    fmpz_clear(c + 1);
    fmpz_clear(c + 2);

    fmpz_clear(sol + 0);
    fmpz_clear(sol + 1);
    fmpz_clear(sol + 2);

    fmpz_clear(tr[0]); 
    fmpz_clear(tr[1]); 
    fmpz_clear(tr[2]);

    // init qbf
    fmpz_tqf_clear(TBF);
    fmpz_tqf_clear(TBF_r);

    // init factors
    fmpz_factor_clear(af);
    fmpz_factor_clear(bf);
    fmpz_factor_clear(cf);
}
