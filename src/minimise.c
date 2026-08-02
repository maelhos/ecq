#include "minimise.h"
#include "reduce.h"

void bkf_invariants(fmpz_t I, fmpz_t J, const bkf_t C)
{
    fmpz_t tmp;

    fmpz_init(tmp);

    // I = 12ae - 3bd + c^2
    fmpz_mul(I, C->a, C->e);
    fmpz_mul_ui(I, I, 12);
    fmpz_mul(tmp, C->b, C->d);
    fmpz_submul_ui(I, tmp, 3);
    fmpz_addmul(I, C->c, C->c);

    // J = (72ae + 9bd - 2c^2) c - 27 (a d^2 + b^2 e)
    fmpz_mul(J, C->a, C->e);
    fmpz_mul_ui(J, J, 72);
    fmpz_mul(tmp, C->b, C->d);
    fmpz_addmul_ui(J, tmp, 9);
    fmpz_mul(tmp, C->c, C->c);
    fmpz_submul_ui(J, tmp, 2);
    fmpz_mul(J, J, C->c);
    fmpz_mul(tmp, C->d, C->d);
    fmpz_mul(tmp, tmp, C->a);
    fmpz_submul_ui(J, tmp, 27);
    fmpz_mul(tmp, C->b, C->b);
    fmpz_mul(tmp, tmp, C->e);
    fmpz_submul_ui(J, tmp, 27);

    fmpz_clear(tmp);
}

// M <- M * [[e00, e01], [e10, e11]], matching the convention of bkf_apply_mat
static void mat_mul_right(fmpz_mat_t M, const fmpz_t e00, const fmpz_t e01,
                                        const fmpz_t e10, const fmpz_t e11)
{
    fmpz_mat_t E;

    fmpz_mat_init(E, 2, 2);
    fmpz_set(fmpz_mat_entry(E, 0, 0), e00);
    fmpz_set(fmpz_mat_entry(E, 0, 1), e01);
    fmpz_set(fmpz_mat_entry(E, 1, 0), e10);
    fmpz_set(fmpz_mat_entry(E, 1, 1), e11);
    fmpz_mat_mul(M, M, E);
    fmpz_mat_clear(E);
}

static void mat_mul_right_si(fmpz_mat_t M, slong e00, slong e01, slong e10, slong e11)
{
    fmpz_t a, b, c, d;

    fmpz_init_set_si(a, e00); fmpz_init_set_si(b, e01);
    fmpz_init_set_si(c, e10); fmpz_init_set_si(d, e11);
    mat_mul_right(M, a, b, c, d);
    fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(d);
}

// records the substitution x -> k x; the coefficient arithmetic is done by the caller
static void mat_x_scale(fmpz_mat_t M, const fmpz_t k)
{
    fmpz_t zero, one;

    fmpz_init(zero);
    fmpz_init_set_ui(one, 1);
    mat_mul_right(M, k, zero, zero, one);
    fmpz_clear(zero);
    fmpz_clear(one);
}

// C(x, z) <- C(x + alpha z, z)
static void q_xshift(bkf_t C, const fmpz_t alpha, fmpz_mat_t M)
{
    fmpz_t t, zero, one;

    if (fmpz_is_zero(alpha))
    {
        return;
    }

    fmpz_init(t);
    fmpz_init(zero);
    fmpz_init_set_ui(one, 1);

    // e += alpha (d + alpha (c + alpha (b + alpha a)))
    fmpz_mul(t, alpha, C->a);
    fmpz_add(t, t, C->b);
    fmpz_mul(t, t, alpha);
    fmpz_add(t, t, C->c);
    fmpz_mul(t, t, alpha);
    fmpz_add(t, t, C->d);
    fmpz_mul(t, t, alpha);
    fmpz_add(C->e, C->e, t);

    // d += alpha (2c + alpha (3b + 4 alpha a))
    fmpz_mul(t, alpha, C->a);
    fmpz_mul_ui(t, t, 4);
    fmpz_addmul_ui(t, C->b, 3);
    fmpz_mul(t, t, alpha);
    fmpz_addmul_ui(t, C->c, 2);
    fmpz_mul(t, t, alpha);
    fmpz_add(C->d, C->d, t);

    // c += alpha (3b + 6 alpha a)
    fmpz_mul(t, alpha, C->a);
    fmpz_mul_ui(t, t, 6);
    fmpz_addmul_ui(t, C->b, 3);
    fmpz_mul(t, t, alpha);
    fmpz_add(C->c, C->c, t);

    // b += 4 alpha a
    fmpz_mul(t, alpha, C->a);
    fmpz_mul_ui(t, t, 4);
    fmpz_add(C->b, C->b, t);

    mat_mul_right(M, one, alpha, zero, one);

    fmpz_clear(t);
    fmpz_clear(zero);
    fmpz_clear(one);
}

// C(x, z) <- C(x, z + gamma x)
static void q_zshift(bkf_t C, const fmpz_t gamma, fmpz_mat_t M)
{
    fmpz_t t, zero, one;

    if (fmpz_is_zero(gamma))
    {
        return;
    }

    fmpz_init(t);
    fmpz_init(zero);
    fmpz_init_set_ui(one, 1);

    // a += gamma (b + gamma (c + gamma (d + gamma e)))
    fmpz_mul(t, gamma, C->e);
    fmpz_add(t, t, C->d);
    fmpz_mul(t, t, gamma);
    fmpz_add(t, t, C->c);
    fmpz_mul(t, t, gamma);
    fmpz_add(t, t, C->b);
    fmpz_mul(t, t, gamma);
    fmpz_add(C->a, C->a, t);

    // b += gamma (2c + gamma (3d + 4 gamma e))
    fmpz_mul(t, gamma, C->e);
    fmpz_mul_ui(t, t, 4);
    fmpz_addmul_ui(t, C->d, 3);
    fmpz_mul(t, t, gamma);
    fmpz_addmul_ui(t, C->c, 2);
    fmpz_mul(t, t, gamma);
    fmpz_add(C->b, C->b, t);

    // c += gamma (3d + 6 gamma e)
    fmpz_mul(t, gamma, C->e);
    fmpz_mul_ui(t, t, 6);
    fmpz_addmul_ui(t, C->d, 3);
    fmpz_mul(t, t, gamma);
    fmpz_add(C->c, C->c, t);

    // d += 4 gamma e
    fmpz_mul(t, gamma, C->e);
    fmpz_mul_ui(t, t, 4);
    fmpz_add(C->d, C->d, t);

    mat_mul_right(M, one, zero, gamma, one);

    fmpz_clear(t);
    fmpz_clear(zero);
    fmpz_clear(one);
}

// C(x, z) <- C(-z, x) : brings a multiple root at infinity down to 0
static void q_invert(bkf_t C, fmpz_mat_t M)
{
    fmpz_swap(C->a, C->e);
    fmpz_swap(C->b, C->d);
    fmpz_neg(C->b, C->b);
    fmpz_neg(C->d, C->d);

    mat_mul_right_si(M, 0, -1, 1, 0);
}

/* Assuming p | I and p | J, the unique alpha mod p modulo which the quartic has a
   root of multiplicity at least 3. alpha = -1 means that root is at infinity. */
static void root_p(fmpz_t alpha, const bkf_t C, const fmpz_t p)
{
    fmpz_t b2, ac, psem, rsem, t, tmp;

    if (fmpz_divisible(C->a, p) && fmpz_divisible(C->b, p))
    {
        fmpz_set_si(alpha, -1);
        return;
    }
    if (fmpz_divisible(C->e, p) && fmpz_divisible(C->d, p))
    {
        fmpz_zero(alpha);
        return;
    }
    if (fmpz_cmp_ui(p, 2) == 0)
    {
        fmpz_one(alpha);
        return;
    }

    fmpz_init(b2); fmpz_init(ac); fmpz_init(psem);
    fmpz_init(rsem); fmpz_init(t); fmpz_init(tmp);

    if (fmpz_cmp_ui(p, 3) == 0)
    {
        if (fmpz_divisible(C->a, p))
        {
            fmpz_mul(alpha, C->b, C->e);
        }
        else
        {
            fmpz_mul(alpha, C->a, C->d);
        }
        fmpz_neg(alpha, alpha);
        fmpz_mod(alpha, alpha, p);
        goto cleanup;
    }

    fmpz_mul(b2, C->b, C->b);
    fmpz_mul(ac, C->a, C->c);

    // p_seminv = 3b^2 - 8ac
    fmpz_mul_ui(psem, b2, 3);
    fmpz_submul_ui(psem, ac, 8);
    fmpz_mod(psem, psem, p);

    if (fmpz_is_zero(psem))
    {
        // quadruple root : alpha = -b / (4a)
        fmpz_mul_ui(t, C->a, 4);
        fmpz_invmod(t, t, p);
        fmpz_neg(alpha, C->b);
        fmpz_mul(alpha, alpha, t);
        fmpz_mod(alpha, alpha, p);
        goto cleanup;
    }

    if (fmpz_divisible(C->a, p))
    {
        // fourth root at infinity : alpha = -c / (3b)
        fmpz_mul_ui(t, C->b, 3);
        fmpz_invmod(t, t, p);
        fmpz_neg(alpha, C->c);
        fmpz_mul(alpha, alpha, t);
        fmpz_mod(alpha, alpha, p);
        goto cleanup;
    }

    // triple root : alpha = (3 r_seminv - b p_seminv) / (4 a p_seminv)
    fmpz_mul_ui(t, C->a, 4);
    fmpz_mul(t, t, psem);
    fmpz_invmod(t, t, p);

    fmpz_mul(rsem, C->b, b2);
    fmpz_mul(tmp, C->a, C->a);
    fmpz_mul(tmp, tmp, C->d);
    fmpz_addmul_ui(rsem, tmp, 8);
    fmpz_mul(tmp, ac, C->b);
    fmpz_submul_ui(rsem, tmp, 4);

    fmpz_mul_ui(alpha, rsem, 3);
    fmpz_submul(alpha, C->b, psem);
    fmpz_mul(alpha, alpha, t);
    fmpz_mod(alpha, alpha, p);

cleanup:
    fmpz_clear(b2); fmpz_clear(ac); fmpz_clear(psem);
    fmpz_clear(rsem); fmpz_clear(t); fmpz_clear(tmp);
}

static int all_divisible(const bkf_t C, const fmpz_t k)
{
    return fmpz_divisible(C->a, k) && fmpz_divisible(C->b, k)
        && fmpz_divisible(C->c, k) && fmpz_divisible(C->d, k)
        && fmpz_divisible(C->e, k);
}

static void divide_all(bkf_t C, const fmpz_t k)
{
    fmpz_divexact(C->a, C->a, k);
    fmpz_divexact(C->b, C->b, k);
    fmpz_divexact(C->c, C->c, k);
    fmpz_divexact(C->d, C->d, k);
    fmpz_divexact(C->e, C->e, k);
}

/* One minimisation step at p, assuming p^4 | I and p^6 | J, dividing those by p^4
   and p^6. Returns 0 only for p = 2, where the step can genuinely fail. */
static int minim_p(bkf_t C, const fmpz_t p, fmpz_mat_t M)
{
    fmpz_t p2, p3, p4, p6, alpha, gamma, t, b0, sixteen;
    bkf_t C0;
    int ret = 1;
    int p_is_2 = (fmpz_cmp_ui(p, 2) == 0);
    int p_is_3 = (fmpz_cmp_ui(p, 3) == 0);

    fmpz_init(p2); fmpz_init(p3); fmpz_init(p4); fmpz_init(p6);
    fmpz_init(alpha); fmpz_init(gamma); fmpz_init(t); fmpz_init(b0);
    fmpz_init_set_ui(sixteen, 16);
    bkf_init(C0);

    fmpz_mul(p2, p, p);
    fmpz_mul(p3, p2, p);
    fmpz_mul(p4, p2, p2);
    fmpz_mul(p6, p4, p2);

    // trivial case: p^2 divides every coefficient
    if (all_divisible(C, p2))
    {
        divide_all(C, p2);
        goto cleanup;
    }

    if (all_divisible(C, p))
    {
        // locate the multiple root using the quotient by p
        bkf_set(C0, C->a, C->b, C->c, C->d, C->e);
        divide_all(C0, p);

        if (fmpz_divisible(C0->a, p) && fmpz_divisible(C0->b, p))
        {
            q_invert(C, M);
        }
        else if (!fmpz_divisible(C0->e, p) || !fmpz_divisible(C0->d, p))
        {
            root_p(alpha, C0, p);
            q_xshift(C, alpha, M);
        }

        fmpz_divexact(b0, C->b, p);
        if (!fmpz_divisible(C->a, p2) && !fmpz_divisible(b0, p))
        {
            // shift the fourth root to infinity
            fmpz_divexact(t, C->a, p);
            fmpz_invmod(gamma, b0, p);
            fmpz_mul(gamma, gamma, t);
            fmpz_neg(gamma, gamma);
            fmpz_mod(gamma, gamma, p);
            q_zshift(C, gamma, M);
            fmpz_divexact(b0, C->b, p);
        }

        if (p_is_2 && !fmpz_divisible(C->e, sixteen))
        {
            ret = 0;
            goto cleanup;
        }

        fmpz_set(C->b, b0);
        fmpz_divexact(C->c, C->c, p2);
        fmpz_divexact(C->d, C->d, p3);
        fmpz_divexact(C->e, C->e, p4);
        mat_x_scale(M, p);
        goto cleanup;
    }

    // not all coefficients divisible by p
    if (fmpz_divisible(C->a, p) && fmpz_divisible(C->b, p))
    {
        q_invert(C, M);
    }
    else if (!fmpz_divisible(C->e, p) || !fmpz_divisible(C->d, p))
    {
        root_p(alpha, C, p);
        q_xshift(C, alpha, M);
    }

    if (!fmpz_divisible(C->a, p) && !fmpz_divisible(C->b, p))
    {
        fmpz_invmod(gamma, C->b, p);
        fmpz_mul(gamma, gamma, C->a);
        fmpz_neg(gamma, gamma);
        fmpz_mod(gamma, gamma, p);
        q_zshift(C, gamma, M);
    }

    if (fmpz_divisible(C->a, p))
    {
        // triple root case
        fmpz_zero(alpha);

        if (p_is_3)
        {
            fmpz_t I, J;
            slong vpi;

            fmpz_init(I); fmpz_init(J);
            bkf_invariants(I, J, C);
            vpi = fmpz_is_zero(I) ? 1000 : (slong) fmpz_remove(t, I, p);

            if (vpi == 4)
            {
                fmpz_neg(t, C->e);
                fmpz_divexact_ui(t, t, 27);
            }
            else
            {
                fmpz_neg(t, C->c);
                fmpz_divexact_ui(t, t, 9);
            }

            if (!fmpz_divisible(t, p))
            {
                fmpz_invmod(alpha, C->b, p);
                fmpz_mul(alpha, alpha, t);
                fmpz_mod(alpha, alpha, p);
                fmpz_mul_ui(alpha, alpha, 3);
            }

            fmpz_clear(I); fmpz_clear(J);
        }
        else
        {
            fmpz_neg(t, C->c);
            fmpz_divexact(t, t, p);
            if (!fmpz_divisible(t, p))
            {
                fmpz_mul_ui(alpha, C->b, 3);
                fmpz_invmod(alpha, alpha, p);
                fmpz_mul(alpha, alpha, t);
                fmpz_mod(alpha, alpha, p);
                fmpz_mul(alpha, alpha, p);
            }
        }

        q_xshift(C, alpha, M);

        fmpz_mul(C->a, C->a, p2);
        fmpz_divexact(C->c, C->c, p2);
        fmpz_divexact(C->d, C->d, p4);
        fmpz_divexact(C->e, C->e, p6);
        mat_x_scale(M, p2);
        goto cleanup;
    }

    // quadruple root case
    if (p_is_3)
    {
        fmpz_divexact_ui(b0, C->b, 3);
        if (!fmpz_divisible(b0, p))
        {
            fmpz_invmod(alpha, C->a, p);
            fmpz_mul(alpha, alpha, b0);
            fmpz_neg(alpha, alpha);
            fmpz_mod(alpha, alpha, p);
            fmpz_mul_ui(alpha, alpha, 3);
            q_xshift(C, alpha, M);
        }
    }

    if (p_is_2 && !fmpz_divisible(C->e, sixteen))
    {
        ret = 0;
        goto cleanup;
    }

    fmpz_divexact(C->b, C->b, p);
    fmpz_divexact(C->c, C->c, p2);
    fmpz_divexact(C->d, C->d, p3);
    fmpz_divexact(C->e, C->e, p4);
    mat_x_scale(M, p);

cleanup:
    fmpz_clear(p2); fmpz_clear(p3); fmpz_clear(p4); fmpz_clear(p6);
    fmpz_clear(alpha); fmpz_clear(gamma); fmpz_clear(t); fmpz_clear(b0);
    fmpz_clear(sixteen);
    bkf_clear(C0);
    return ret;
}

// vp(I) > 3 and vp(J) > 5, with the sharper condition at 3. Local solubility is
// assumed, which is legitimate since we only minimise ELS descendants.
static int is_nonmin(const fmpz_t p, slong vpi, slong vpj, slong vpd)
{
    if (fmpz_cmp_ui(p, 3) == 0)
    {
        return ((vpi > 4) && (vpj > 8)) || ((vpi == 4) && (vpj == 6) && (vpd > 14));
    }
    return (vpi > 3) && (vpj > 5);
}

slong bkf_minimise(bkf_t C, fmpz_mat_t M)
{
    fmpz_t I, J, g, tmp, disc;
    fmpz_factor_t fac;
    slong i, steps = 0;

    fmpz_mat_one(M);

    fmpz_init(I); fmpz_init(J); fmpz_init(g);
    fmpz_init(tmp); fmpz_init(disc);
    fmpz_factor_init(fac);

    bkf_invariants(I, J, C);

    if (fmpz_is_zero(I) && fmpz_is_zero(J))
    {
        goto cleanup;
    }

    // p^4 | I and p^6 | J force p | gcd(I, J), so only those primes can occur
    fmpz_gcd(g, I, J);
    if (fmpz_is_zero(g) || fmpz_is_pm1(g))
    {
        goto cleanup;
    }

    fmpz_factor(fac, g);

    for (i = 0; i < fac->num; i++)
    {
        const fmpz *p = fac->p + i;
        slong vpi, vpj, vpd = 0;

        vpi = fmpz_is_zero(I) ? 1000 : (slong) fmpz_remove(tmp, I, p);
        vpj = fmpz_is_zero(J) ? 1000 : (slong) fmpz_remove(tmp, J, p);

        if (fmpz_cmp_ui(p, 3) == 0)
        {
            // disc = 4 I^3 - J^2
            fmpz_pow_ui(disc, I, 3);
            fmpz_mul_ui(disc, disc, 4);
            fmpz_mul(tmp, J, J);
            fmpz_sub(disc, disc, tmp);
            vpd = fmpz_is_zero(disc) ? 1000 : (slong) fmpz_remove(tmp, disc, p);
        }

        while (is_nonmin(p, vpi, vpj, vpd))
        {
            if (!minim_p(C, p, M))
            {
                break;   // can only happen at p = 2
            }

            steps++;
            vpi -= 4;
            vpj -= 6;
            vpd -= 12;

            bkf_invariants(I, J, C);
        }
    }

cleanup:
    fmpz_clear(I); fmpz_clear(J); fmpz_clear(g);
    fmpz_clear(tmp); fmpz_clear(disc);
    fmpz_factor_clear(fac);
    return steps;
}

int bkf_equivalent(const bkf_t C1, const bkf_t C2, fmpz_mat_t M)
{
    bkf_t A, B, T;
    fmpz_mat_t M1, M2, U;
    fmpz_t I1, J1, I2, J2;
    slong p, r, B0 = BKF_EQUIV_SEARCH;
    int found = 0;

    bkf_init(A); bkf_init(B); bkf_init(T);
    fmpz_mat_init(M1, 2, 2);
    fmpz_mat_init(M2, 2, 2);
    fmpz_mat_init(U, 2, 2);
    fmpz_init(I1); fmpz_init(J1); fmpz_init(I2); fmpz_init(J2);

    bkf_set(A, C1->a, C1->b, C1->c, C1->d, C1->e);
    bkf_set(B, C2->a, C2->b, C2->c, C2->d, C2->e);

    // put both in minimal reduced shape first
    bkf_minimise(A, M1);
    bkf_reduce(A, M2);
    fmpz_mat_mul(M1, M1, M2);

    bkf_minimise(B, M2);
    {
        fmpz_mat_t M3;
        fmpz_mat_init(M3, 2, 2);
        bkf_reduce(B, M3);
        fmpz_mat_mul(M2, M2, M3);
        fmpz_mat_clear(M3);
    }

    bkf_invariants(I1, J1, A);
    bkf_invariants(I2, J2, B);

    // same covering requires the same minimal invariants
    if (!fmpz_equal(I1, I2) || !fmpz_equal(J1, J2))
    {
        goto cleanup;
    }

    // Both are minimal and reduced, so an equivalence between them is a unimodular
    // matrix with small entries. Rather than looping over all four entries, enumerate
    // the first column (p, r) with gcd(p, r) = 1, solve p s - q r = 1 for one (q, s)
    // with xgcd, and slide along the remaining one-parameter family (q, s) + t (p, r).
    // That is O(B^2 T) instead of O(B^4), so B can be taken much larger.
    {
        fmpz_t fp, fr, g, s0, q0, qq, ss;

        fmpz_init(fp); fmpz_init(fr); fmpz_init(g);
        fmpz_init(s0); fmpz_init(q0); fmpz_init(qq); fmpz_init(ss);

        for (p = -B0; p <= B0 && !found; p++)
        for (r = -B0; r <= B0 && !found; r++)
        {
            slong t;

            if (p == 0 && r == 0)
            {
                continue;
            }

            fmpz_set_si(fp, p);
            fmpz_set_si(fr, r);
            fmpz_xgcd(g, s0, q0, fp, fr);      // g = s0 p + q0 r
            if (!fmpz_is_one(g))
            {
                continue;
            }
            fmpz_neg(q0, q0);                  // p s0 - q0 r = 1

            for (t = -BKF_EQUIV_SLIDE; t <= BKF_EQUIV_SLIDE && !found; t++)
            {
                fmpz_set_si(qq, t);
                fmpz_mul(qq, qq, fp);
                fmpz_add(qq, qq, q0);

                fmpz_set_si(ss, t);
                fmpz_mul(ss, ss, fr);
                fmpz_add(ss, ss, s0);

                int flip;

                // both determinants: xgcd only produces det +1, but an equivalence
                // of coverings may just as well reverse orientation
                for (flip = 0; flip < 2 && !found; flip++)
                {
                    fmpz_set(fmpz_mat_entry(U, 0, 0), fp);
                    fmpz_set(fmpz_mat_entry(U, 0, 1), qq);
                    fmpz_set(fmpz_mat_entry(U, 1, 0), fr);
                    fmpz_set(fmpz_mat_entry(U, 1, 1), ss);

                    if (flip)
                    {
                        fmpz_neg(fmpz_mat_entry(U, 0, 1), fmpz_mat_entry(U, 0, 1));
                        fmpz_neg(fmpz_mat_entry(U, 1, 1), fmpz_mat_entry(U, 1, 1));
                    }

                    bkf_set(T, A->a, A->b, A->c, A->d, A->e);
                    bkf_apply_mat(T, U);

                    if (fmpz_equal(T->a, B->a) && fmpz_equal(T->b, B->b) && fmpz_equal(T->c, B->c)
                     && fmpz_equal(T->d, B->d) && fmpz_equal(T->e, B->e))
                    {
                        found = 1;
                        fmpz_mat_mul(M, M1, U);
                    }
                }
            }
        }

        fmpz_clear(fp); fmpz_clear(fr); fmpz_clear(g);
        fmpz_clear(s0); fmpz_clear(q0); fmpz_clear(qq); fmpz_clear(ss);
    }

cleanup:
    bkf_clear(A); bkf_clear(B); bkf_clear(T);
    fmpz_mat_clear(M1); fmpz_mat_clear(M2); fmpz_mat_clear(U);
    fmpz_clear(I1); fmpz_clear(J1); fmpz_clear(I2); fmpz_clear(J2);
    return found;
}
