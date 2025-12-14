#include "bkf.h"
#include <stdlib.h>

void bkf_power_precompute(fmpz_t U4, fmpz_t U3V1, fmpz_t U2V2, fmpz_t U1V3, fmpz_t V4, long U, long V)
{
    fmpz_t U2, V2;

    fmpz_init(U2);
    fmpz_init(V2);

    ulong pU = labs(U);
    ulong pV = labs(V);

    fmpz_ui_mul_ui(U2, pU, pU);
    fmpz_ui_mul_ui(V2, pV, pV);
    fmpz_mul(U2V2, U2, V2);
    fmpz_mul(U4, U2, U2);
    fmpz_mul(V4, V2, V2);

    fmpz_mul_si(U3V1, U2, U);
    fmpz_mul_si(U3V1, U3V1, V);

    fmpz_mul_si(U1V3, V2, V);
    fmpz_mul_si(U1V3, U1V3, U);

    fmpz_clear(U2);
    fmpz_clear(V2);
}

int bkf_power_bounded(fmpz_t test, const bkf_t C, const fmpz_t U4, const fmpz_t U3V1, const fmpz_t U2V2, const fmpz_t U1V3, const fmpz_t V4)
{   
    fmpz_t T_P, T_I, tmp;
    int ret;

    fmpz_init(T_P);
    fmpz_init(T_I);
    fmpz_init(tmp);

    // evaluate the bkf at U, V
    fmpz_mul(T_P, U4, C->a);

    fmpz_mul(T_I, U3V1, C->b);

    fmpz_mul(tmp, U2V2, C->c);
    fmpz_add(T_P, T_P, tmp);

    fmpz_mul(tmp,  U1V3, C->d);
    fmpz_add(T_I, T_I, tmp);

    fmpz_mul(tmp,  V4, C->e);
    fmpz_add(T_P, T_P, tmp);

    fmpz_add(test, T_P, T_I);

    ret = 0;

    if (fmpz_is_square(test) && !fmpz_is_zero(test))
    {
        ret = 1;
        goto CLEAR;
    }

    // evaluate the bkf at U, -V
    fmpz_sub(test, T_P, T_I);

    if (fmpz_is_square(test) && !fmpz_is_zero(test))
    {
        ret = -1;
    }

    CLEAR:;
    fmpz_clear(T_P);
    fmpz_clear(T_I);
    fmpz_clear(tmp);

    return ret;
}

int bkf_local_solubility(bkf_t C, const fmpz_t delta_curve, const fmpz_factor_t delta_curve_f)
{
    if (!bkf_real_solubility(C))
    {
        return 0;
    }

    fmpz_t p, tmp, disc;
    slong i;

    fmpz_init_set_ui(p, 2);

    if (!bkf_Qp_solubility(C, p))
    {
        fmpz_clear(p);
        return 0;
    }

    fmpz_clear(p);
    bkf_disc(disc, C);
    fmpz_init(tmp);
    if (!fmpz_divides(tmp, disc, delta_curve))
    {
        fmpz_clear(tmp);
        return 0;
    }

    fmpz_clear(tmp);
    fmpz_factor_init(C->disc);
    if (fmpz_factor_over_base(C->disc, disc, delta_curve_f) != 1)
    {
        fmpz_factor_clear(C->disc);
        return 0;
    }

    for (i = 0; i < C->disc->num; i++)
    {
        if (fmpz_cmp_ui(C->disc->p + i, 2) == 0)
        {
            continue;
        }

        if (!bkf_Qp_solubility(C, C->disc->p + i))
        {
            fmpz_factor_clear(C->disc);
            return 0;
        }
    }

    return 1;
}

int bkf_real_solubility(const bkf_t C)
{
    if (fmpz_sgn(C->a) > 0)
    {
        return 1;
    }

    fmpz_poly_t P;
    fmpz_poly_factor_t f;
    acb_ptr v;
    slong i, j;

    v = _acb_vec_init(4);
    fmpz_poly_init2(P, 5);
    fmpz_poly_factor_init(f);

    fmpz_poly_set_coeff_fmpz(P, 0, C->e);
    fmpz_poly_set_coeff_fmpz(P, 1, C->d);
    fmpz_poly_set_coeff_fmpz(P, 2, C->c);
    fmpz_poly_set_coeff_fmpz(P, 3, C->b);
    fmpz_poly_set_coeff_fmpz(P, 4, C->a);
    
    fmpz_poly_factor_squarefree(f, P);

    for (i = 0; i < f->num; i++)
    {
        arb_fmpz_poly_complex_roots(v, f->p + i, 0, 256);
 
        for (j = 0; j < fmpz_poly_degree(f->p + i); j++)
        {
            if (acb_is_real(v + j))
            {
                return 1;
            }
        }
    }   

    _acb_vec_clear(v, 4);
    fmpz_poly_clear(P);
    fmpz_poly_factor_clear(f);

    return 0;
}

void bkf_disc(fmpz_t disc, const bkf_t C)
{
    fmpz_t ae, bd, c2, I, tmp, J, D;

    // I = 12ae - 3bd + c**2
    // J = 72ace + 9bcd - 27ad**2 - 27eb**2 - 2c**3
    //   = c*(72ae + 9bd - 2c**2) - 27a*d**2 - 27*e*b**2
    fmpz_init(ae);
    fmpz_init(bd);
    fmpz_init(c2);
    fmpz_init(I);
    fmpz_init(tmp);
    fmpz_init(J);
    fmpz_init(D);

    fmpz_mul(ae, C->a, C->e);
    fmpz_mul(bd, C->b, C->d);
    fmpz_mul(c2, C->c, C->c);

    fmpz_mul_ui(I, ae, 12);
    fmpz_submul_ui(I, bd, 3);
    fmpz_add(I, I, c2);

    fmpz_mul_ui(tmp, ae, 72);
    fmpz_addmul_ui(tmp, bd, 9);
    fmpz_submul_ui(tmp, c2, 2);
    fmpz_mul(J, C->c, tmp);

    fmpz_mul(tmp, C->d, C->d);
    fmpz_mul(tmp, tmp, C->a);
    fmpz_submul_ui(J, tmp, 27);

    fmpz_mul(tmp, C->b, C->b);
    fmpz_mul(tmp, tmp, C->e);
    fmpz_submul_ui(J, tmp, 27);

    fmpz_pow_ui(D, I, 3);
    fmpz_mul_ui(D, D, 4);
    fmpz_mul(tmp, J, J);
    fmpz_sub(D, D, tmp);

    fmpz_divexact_ui(disc, D, 27);

    fmpz_clear(ae);
    fmpz_clear(bd);
    fmpz_clear(c2);
    fmpz_clear(I);
    fmpz_clear(tmp);
    fmpz_clear(J);
    fmpz_clear(D);
}

int bkf_Qp_solubility(const bkf_t C, const fmpz_t p)
{
    fmpz_t zero;

    fmpz_init_set_ui(zero, 0);

    if (bkf_Zp_solubility(C, p, zero, 0))
    {
        return 1;
    }

    bkf_t rev;

    bkf_init_set(rev, C->e, C->d, C->c, C->b, C->a);

    if (bkf_Zp_solubility(rev, p, zero, 1))
    {
        fmpz_clear(zero);
        bkf_clear(rev);
        return 1;
    }

    fmpz_clear(zero);
    bkf_clear(rev);
    return 0;
}

static inline ulong ord(const fmpz_t a, const fmpz_t p)
{
    fmpz_t tmp;
    ulong cnt;

    cnt = 0;

    fmpz_init_set(tmp, a);

    while (fmpz_divides(tmp, tmp, p))
    {
        cnt++;
    }

    fmpz_clear(tmp);

    return cnt;
}

static inline int lemma7(const bkf_t C, const fmpz_t x, ulong n)
{
    fmpz_t tmp, tmp2;
    ulong val2, gx_odd;
    ulong l, m;

    fmpz_init_set(tmp, C->a);
    fmpz_init(tmp2);

    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->b);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->c);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->d);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->e);

    // this corresponds to being a 2-adic square
    val2 = fmpz_val2(tmp);
    if (val2 % 2 == 0)
    {
        fmpz_fdiv_r_2exp(tmp2, tmp, 3);

        if (fmpz_cmp_ui(tmp2, 1) == 0)
        {
            return 1;
        }
    }

    fmpz_fdiv_q_2exp(tmp, tmp, val2);
    gx_odd = fmpz_fdiv_ui(tmp, 4);
    l = fmpz_val2(tmp);

    fmpz_mul_ui(tmp, C->a, 4);
    fmpz_mul(tmp, tmp, x);
    fmpz_addmul_ui(tmp, C->b, 3);
    fmpz_mul(tmp, tmp, x);
    fmpz_addmul_ui(tmp, C->c, 2);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->d);

    m = fmpz_val2(tmp);
    fmpz_clear(tmp);
    fmpz_clear(tmp2);

    if (l >= m + n && n > m) return 1;
    if (n > m && l == m + n - 1 && l % 2 == 0) return 1;
    if (n > m && l == m + n - 2 && gx_odd == 1 && l % 2 == 0) return 1;
    if (m >= n && l >= 2*n) return 0;
    if (m >= n && l == 2*n - 2 && gx_odd == 1) return 0;

    return -1;
}

static inline int lemma6(const bkf_t C, const fmpz_t x, const fmpz_t p, ulong n)
{
    fmpz_t tmp;
    ulong l, m;

    fmpz_init_set(tmp, C->a);
    
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->b);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->c);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->d);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->e);

    if (fmpz_jacobi(tmp, p) >= 0)
    {
        return 1;
    }

    l = ord(tmp, p);

    fmpz_mul_ui(tmp, C->a, 4);
    fmpz_mul(tmp, tmp, x);
    fmpz_addmul_ui(tmp, C->b, 3);
    fmpz_mul(tmp, tmp, x);
    fmpz_addmul_ui(tmp, C->c, 2);
    fmpz_mul(tmp, tmp, x);
    fmpz_add(tmp, tmp, C->d);

    m = ord(tmp, p);
    fmpz_clear(tmp);

    if (l >= m + n && n > m) return 1;
    if (l >= 2*n && m >= n) return 0;

    return -1;
}

int bkf_Zp_solubility(const bkf_t C, const fmpz_t p, const fmpz_t x, ulong k)
{
    int code;
    fmpz_t t, tmp;

    if (fmpz_cmp_ui(p, 2) == 0)
    {
        code = lemma7(C, x, k);
    }
    else
    {
        code = lemma6(C, x, p, k);
    }

    if (code == 1)
    {
        return 1;
    }
    
    if (code == -1)
    {
        return 0;   
    }

    fmpz_init(t);
    fmpz_init(tmp);
    
    for (fmpz_set_ui(t, 0); fmpz_cmp(t, p) < 0; fmpz_add_ui(t, t, 1))
    {
        fmpz_pow_ui(tmp, p, k);
        fmpz_mul(tmp, tmp, t);
        fmpz_add(tmp, tmp, x);

        if (bkf_Zp_solubility(C, p, tmp, k + 1))
        {
            fmpz_clear(t);
            fmpz_clear(tmp);

            return 1;
        }
    }

    fmpz_clear(t);
    fmpz_clear(tmp);

    return 0;
}

void bkf_reduce_naive(bkf_t C)
{
    fmpz_t tmp, u;
    fmpz_factor_t f;
    slong i;

    fmpz_init(tmp);
    fmpz_init(u);
    fmpz_factor_init(f);

    fmpz_gcd3(tmp, C->a, C->b, C->c);
    fmpz_gcd3(tmp, tmp, C->d, C->e);
    fmpz_factor(f, tmp);

    fmpz_set_ui(u, 1);

    for (i = 0; i < f->num; i++)
    {
        if (f->exp[i] > 1)
        {
            fmpz_pow_ui(tmp, f->p + i, 2*(f->exp[i] / 2));
            fmpz_mul(u, u, tmp);
        }
    }

    fmpz_divexact(C->a, C->a, u);
    fmpz_divexact(C->b, C->b, u);
    fmpz_divexact(C->c, C->c, u);
    fmpz_divexact(C->d, C->d, u);
    fmpz_divexact(C->e, C->e, u);
    
    fmpz_clear(u);
    fmpz_clear(tmp);
}