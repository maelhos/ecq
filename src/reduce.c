#include "reduce.h"

#include <flint/fmpz_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>

#define REDUCE_PREC 256 // precision used in arb to solve over R, for now fixed, maybe a bad idea
#define REDUCE_MAX_STEPS 200

/* Substituting (x, y) -> (p x + q y, r x + s y) */
static void form_apply_mat(fmpz* out, const fmpz* in, slong n, const fmpz_mat_t M)
{
    fmpz_poly_t P1, P2, acc, term, pw1, pw2;
    slong i;

    fmpz_poly_init(P1);
    fmpz_poly_init(P2);
    fmpz_poly_init(acc);
    fmpz_poly_init(term);
    fmpz_poly_init(pw1);
    fmpz_poly_init(pw2);

    // P1 = p x + q, P2 = r x + s
    fmpz_poly_set_coeff_fmpz(P1, 1, fmpz_mat_entry(M, 0, 0));
    fmpz_poly_set_coeff_fmpz(P1, 0, fmpz_mat_entry(M, 0, 1));
    fmpz_poly_set_coeff_fmpz(P2, 1, fmpz_mat_entry(M, 1, 0));
    fmpz_poly_set_coeff_fmpz(P2, 0, fmpz_mat_entry(M, 1, 1));

    // in[i] is the coefficient of x^(n-i) y^i
    for (i = 0; i <= n; i++)
    {
        fmpz_poly_pow(pw1, P1, n - i);
        fmpz_poly_pow(pw2, P2, i);
        fmpz_poly_mul(term, pw1, pw2);
        fmpz_poly_scalar_mul_fmpz(term, term, in + i);
        fmpz_poly_add(acc, acc, term);
    }

    for (i = 0; i <= n; i++)
    {
        fmpz_poly_get_coeff_fmpz(out + i, acc, n - i);
    }

    fmpz_poly_clear(P1);
    fmpz_poly_clear(P2);
    fmpz_poly_clear(acc);
    fmpz_poly_clear(term);
    fmpz_poly_clear(pw1);
    fmpz_poly_clear(pw2);
}

void bkf_apply_mat(bkf_t C, const fmpz_mat_t M)
{
    fmpz *in, *out;

    in = _fmpz_vec_init(5);
    out = _fmpz_vec_init(5);

    fmpz_set(in + 0, C->a);
    fmpz_set(in + 1, C->b);
    fmpz_set(in + 2, C->c);
    fmpz_set(in + 3, C->d);
    fmpz_set(in + 4, C->e);

    form_apply_mat(out, in, 4, M);

    fmpz_set(C->a, out + 0);
    fmpz_set(C->b, out + 1);
    fmpz_set(C->c, out + 2);
    fmpz_set(C->d, out + 3);
    fmpz_set(C->e, out + 4);

    _fmpz_vec_clear(in, 5);
    _fmpz_vec_clear(out, 5);
}

void qfb_apply_mat(qfb_t q, const fmpz_mat_t M)
{
    fmpz *in, *out;

    in = _fmpz_vec_init(3);
    out = _fmpz_vec_init(3);

    fmpz_set(in + 0, q->a);
    fmpz_set(in + 1, q->b);
    fmpz_set(in + 2, q->c);

    form_apply_mat(out, in, 2, M);

    fmpz_set(q->a, out + 0);
    fmpz_set(q->b, out + 1);
    fmpz_set(q->c, out + 2);

    _fmpz_vec_clear(in, 3);
    _fmpz_vec_clear(out, 3);
}

/* Covariant quadratic for a quartic with three real resolvent roots (types 1 and 2).
   Normalised to (1, a1, a2). */
static void types12_covar(arb_t o1, arb_t o2, const bkf_t C, const arb_t xH, const arb_t phi, slong prec)
{
    arb_t a0, a1, a2, t;
    fmpz_t z, w;

    arb_init(a0); arb_init(a1); arb_init(a2); arb_init(t);
    fmpz_init(z);
    fmpz_init(w);

    // a0 = 3 (4 a phi - H)
    arb_mul_fmpz(a0, phi, C->a, prec);
    arb_mul_ui(a0, a0, 4, prec);
    arb_sub(a0, a0, xH, prec);
    arb_mul_ui(a0, a0, 3, prec);

    // a1 = 6 (b phi + (bc - 6ad))
    arb_mul_fmpz(a1, phi, C->b, prec);
    fmpz_mul(z, C->b, C->c);
    fmpz_mul(w, C->a, C->d);
    fmpz_submul_ui(z, w, 6);
    arb_add_fmpz(a1, a1, z, prec);
    arb_mul_ui(a1, a1, 6, prec);

    // a2 = 2 phi (c - phi) + (4c^2 - 9bd)
    arb_set_fmpz(t, C->c);
    arb_sub(t, t, phi, prec);
    arb_mul(a2, phi, t, prec);
    arb_mul_ui(a2, a2, 2, prec);
    fmpz_mul(z, C->c, C->c);
    fmpz_mul_ui(z, z, 4);
    fmpz_mul(w, C->b, C->d);
    fmpz_submul_ui(z, w, 9);
    arb_add_fmpz(a2, a2, z, prec);

    arb_div(o1, a1, a0, prec);
    arb_div(o2, a2, a0, prec);

    arb_clear(a0); arb_clear(a1); arb_clear(a2); arb_clear(t);
    fmpz_clear(z);
    fmpz_clear(w);
}

/* Julia covariant for a quartic with two real and two complex roots (type 3). */
static void type3_covar(arb_t o1, arb_t o2, const bkf_t C, const arb_t xH,
                        const arb_t rphi, const acb_t cphi, int Risneg, slong prec)
{
    arb_t a4, r3, rr1, ir1, rbeta, ibeta, alpha1, alpha2, p, lambda, q;
    arb_t ar, br, t1sq, t2sq, usq, a0, a1, a2, t;
    acb_t r1, ct, cr3;

    arb_init(a4); arb_init(r3); arb_init(rr1); arb_init(ir1);
    arb_init(rbeta); arb_init(ibeta); arb_init(alpha1); arb_init(alpha2);
    arb_init(p); arb_init(lambda); arb_init(q);
    arb_init(ar); arb_init(br); arb_init(t1sq); arb_init(t2sq); arb_init(usq);
    arb_init(a0); arb_init(a1); arb_init(a2); arb_init(t);
    acb_init(r1); acb_init(ct); acb_init(cr3);

    // a4 = 4a
    arb_set_fmpz(a4, C->a);
    arb_mul_ui(a4, a4, 4, prec);

    // r1 = sqrt((4a*cphi - H)/3)
    acb_mul_arb(ct, cphi, a4, prec);
    acb_sub_arb(ct, ct, xH, prec);
    acb_div_ui(ct, ct, 3, prec);
    acb_sqrt(r1, ct, prec);

    // r3 = +- sqrt(|4a*rphi - H|/3)
    arb_mul(t, a4, rphi, prec);
    arb_sub(t, t, xH, prec);
    arb_abs(t, t);
    arb_div_ui(t, t, 3, prec);
    arb_sqrt(r3, t, prec);
    if (Risneg)
    {
        arb_neg(r3, r3);
    }

    arb_abs(rr1, acb_realref(r1));
    arb_abs(ir1, acb_imagref(r1));

    // rbeta = (r3 - b)/4a, ibeta = 2 ir1 / 4a
    arb_set_fmpz(t, C->b);
    arb_sub(rbeta, r3, t, prec);
    arb_div(rbeta, rbeta, a4, prec);

    arb_mul_ui(ibeta, ir1, 2, prec);
    arb_div(ibeta, ibeta, a4, prec);

    // alpha1 = (2 rr1 - r3 - b)/4a, alpha2 = (-2 rr1 - r3 - b)/4a
    arb_mul_ui(alpha1, rr1, 2, prec);
    arb_sub(alpha1, alpha1, r3, prec);
    arb_sub(alpha1, alpha1, t, prec);
    arb_div(alpha1, alpha1, a4, prec);

    arb_mul_ui(alpha2, rr1, 2, prec);
    arb_neg(alpha2, alpha2);
    arb_sub(alpha2, alpha2, r3, prec);
    arb_sub(alpha2, alpha2, t, prec);
    arb_div(alpha2, alpha2, a4, prec);

    // p = -2 rbeta, lambda = 2 ibeta, q = (p^2 + lambda^2)/4
    arb_mul_si(p, rbeta, -2, prec);
    arb_mul_ui(lambda, ibeta, 2, prec);
    arb_mul(q, p, p, prec);
    arb_addmul(q, lambda, lambda, prec);
    arb_div_ui(q, q, 4, prec);

    // ar = |r1 + r3|, br = |r1 - r3|
    acb_set_arb(cr3, r3);
    acb_add(ct, r1, cr3, prec);
    acb_abs(ar, ct, prec);
    acb_sub(ct, r1, cr3, prec);
    acb_abs(br, ct, prec);

    // t1sq = ir1 ar^2, t2sq = ir1 br^2, usq = rr1 ar br
    arb_mul(t1sq, ar, ar, prec);
    arb_mul(t1sq, t1sq, ir1, prec);
    arb_mul(t2sq, br, br, prec);
    arb_mul(t2sq, t2sq, ir1, prec);
    arb_mul(usq, ar, br, prec);
    arb_mul(usq, usq, rr1, prec);

    // a0 = t1sq + t2sq + 2 usq
    arb_add(a0, t1sq, t2sq, prec);
    arb_addmul_ui(a0, usq, 2, prec);

    // a1 = -2 alpha1 t1sq - 2 alpha2 t2sq + 2 p usq
    arb_mul(a1, alpha1, t1sq, prec);
    arb_mul(t, alpha2, t2sq, prec);
    arb_add(a1, a1, t, prec);
    arb_mul_si(a1, a1, -2, prec);
    arb_mul(t, p, usq, prec);
    arb_addmul_ui(a1, t, 2, prec);

    // a2 = alpha1^2 t1sq + alpha2^2 t2sq + 2 q usq
    arb_mul(a2, alpha1, alpha1, prec);
    arb_mul(a2, a2, t1sq, prec);
    arb_mul(t, alpha2, alpha2, prec);
    arb_mul(t, t, t2sq, prec);
    arb_add(a2, a2, t, prec);
    arb_mul(t, q, usq, prec);
    arb_addmul_ui(a2, t, 2, prec);

    arb_div(o1, a1, a0, prec);
    arb_div(o2, a2, a0, prec);

    arb_clear(a4); arb_clear(r3); arb_clear(rr1); arb_clear(ir1);
    arb_clear(rbeta); arb_clear(ibeta); arb_clear(alpha1); arb_clear(alpha2);
    arb_clear(p); arb_clear(lambda); arb_clear(q);
    arb_clear(ar); arb_clear(br); arb_clear(t1sq); arb_clear(t2sq); arb_clear(usq);
    arb_clear(a0); arb_clear(a1); arb_clear(a2); arb_clear(t);
    acb_clear(r1); acb_clear(ct); acb_clear(cr3);
}

int bkf_reduce(bkf_t C, fmpz_mat_t M)
{
    fmpz_t I, J, H, R, disc, tmp, n;
    fmpz_poly_t cub;
    acb_ptr roots;
    arb_t xH, phi, cov1, cov2, zr, zi, t, u;
    fmpz_mat_t S;
    slong prec, step, i;
    int ret, has_real;

    fmpz_mat_one(M);

    fmpz_init(I); fmpz_init(J); fmpz_init(H); fmpz_init(R);
    fmpz_init(disc); fmpz_init(tmp); fmpz_init(n);

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

    // H = 8ac - 3b^2, R = b^3 + 8 a^2 d - 4abc
    fmpz_mul(H, C->a, C->c);
    fmpz_mul_ui(H, H, 8);
    fmpz_mul(tmp, C->b, C->b);
    fmpz_submul_ui(H, tmp, 3);

    fmpz_pow_ui(R, C->b, 3);
    fmpz_mul(tmp, C->a, C->a);
    fmpz_mul(tmp, tmp, C->d);
    fmpz_addmul_ui(R, tmp, 8);
    fmpz_mul(tmp, C->a, C->b);
    fmpz_mul(tmp, tmp, C->c);
    fmpz_submul_ui(R, tmp, 4);

    // disc = 4 I^3 - J^2
    fmpz_pow_ui(disc, I, 3);
    fmpz_mul_ui(disc, disc, 4);
    fmpz_mul(tmp, J, J);
    fmpz_sub(disc, disc, tmp);

    if (fmpz_is_zero(disc))
    {
        fmpz_clear(I); fmpz_clear(J); fmpz_clear(H); fmpz_clear(R);
        fmpz_clear(disc); fmpz_clear(tmp); fmpz_clear(n);
        return 0;
    }

    // resolvent cubic phi^3 - 3I phi + J
    fmpz_poly_init(cub);
    fmpz_poly_set_coeff_ui(cub, 3, 1);
    fmpz_mul_si(tmp, I, -3);
    fmpz_poly_set_coeff_fmpz(cub, 1, tmp);
    fmpz_poly_set_coeff_fmpz(cub, 0, J);

    arb_init(xH); arb_init(phi); arb_init(cov1); arb_init(cov2);
    arb_init(zr); arb_init(zi); arb_init(t); arb_init(u);
    fmpz_mat_init(S, 2, 2);

    ret = 0;

    for (prec = REDUCE_PREC; prec <= 16 * REDUCE_PREC; prec *= 2)
    {
        roots = _acb_vec_init(3);
        arb_fmpz_poly_complex_roots(roots, cub, 0, prec);

        arb_set_fmpz(xH, H);

        if (fmpz_sgn(disc) > 0)
        {
            // three real roots; order so that a*phi decreases
            arb_ptr rp[3];
            arb_t key[3];

            for (i = 0; i < 3; i++)
            {
                rp[i] = acb_realref(roots + i);
                arb_init(key[i]);
                arb_mul_fmpz(key[i], rp[i], C->a, prec);
            }

            // simple sort of 3 by descending key
            for (i = 0; i < 2; i++)
            {
                slong k;
                for (k = i + 1; k < 3; k++)
                {
                    if (arb_lt(key[i], key[k]))
                    {
                        arb_swap(key[i], key[k]);
                        arb_swap(rp[i], rp[k]);
                    }
                }
            }

            // type 2 iff H < 0 and H^2 > 16 a^2 I, else type 1
            fmpz_mul(tmp, H, H);
            {
                fmpz_t rhs;
                fmpz_init(rhs);
                fmpz_mul(rhs, C->a, C->a);
                fmpz_mul(rhs, rhs, I);
                fmpz_mul_ui(rhs, rhs, 16);

                if ((fmpz_sgn(H) < 0) && (fmpz_cmp(tmp, rhs) > 0))
                {
                    arb_set(phi, rp[1]);   // type 2
                }
                else
                {
                    arb_set(phi, rp[2]);   // type 1
                }
                fmpz_clear(rhs);
            }

            types12_covar(cov1, cov2, C, xH, phi, prec);

            for (i = 0; i < 3; i++)
            {
                arb_clear(key[i]);
            }
        }
        else
        {
            // one real root and a complex conjugate pair
            slong ri = -1, ci = -1;

            for (i = 0; i < 3; i++)
            {
                if (arb_contains_zero(acb_imagref(roots + i)))
                {
                    if (ri < 0) ri = i;
                }
                else if (ci < 0 || arb_is_positive(acb_imagref(roots + i)))
                {
                    if (ci < 0) ci = i;
                }
            }

            if (ri < 0 || ci < 0)
            {
                _acb_vec_clear(roots, 3);
                continue; // retry with more precision
            }

            type3_covar(cov1, cov2, C, xH, acb_realref(roots + ri), roots + ci,
                        fmpz_sgn(R) < 0, prec);
        }

        // covariant quadratic is x^2 + cov1 x + cov2, its root in the upper half
        // plane is z = (-cov1 + i sqrt(4 cov2 - cov1^2)) / 2
        arb_mul(t, cov1, cov1, prec);
        arb_mul_ui(u, cov2, 4, prec);
        arb_sub(u, u, t, prec);

        if (!arb_is_positive(u))
        {
            _acb_vec_clear(roots, 3);
            continue; // not positive definite at this precision, retry
        }

        arb_sqrt(zi, u, prec);
        arb_div_ui(zi, zi, 2, prec);
        arb_neg(zr, cov1);
        arb_div_ui(zr, zr, 2, prec);

        // reduce z to the fundamental domain, accumulating the substitution.
        // z -> z - n  corresponds to M *= [[1, n], [0, 1]]
        // z -> -1/z   corresponds to M *= [[0, 1], [-1, 0]]
        fmpz_mat_one(M);
        has_real = 1;

        for (step = 0; step < REDUCE_MAX_STEPS; step++)
        {
            // n = nearest integer to Re(z)
            if (!arb_get_unique_fmpz(n, zr))
            {
                arb_t half;
                arb_init(half);
                arb_set_d(half, 0.5);
                arb_add(half, zr, half, prec);
                if (!arb_get_unique_fmpz(n, half))
                {
                    arf_get_fmpz(n, arb_midref(half), ARF_RND_FLOOR);
                }
                arb_clear(half);
            }

            if (!fmpz_is_zero(n))
            {
                arb_sub_fmpz(zr, zr, n, prec);

                fmpz_mat_one(S);
                fmpz_set(fmpz_mat_entry(S, 0, 1), n);
                fmpz_mat_mul(M, M, S);
            }

            // |z|^2 = zr^2 + zi^2
            arb_mul(t, zr, zr, prec);
            arb_addmul(t, zi, zi, prec);

            arb_sub_ui(u, t, 1, prec);
            if (!arb_is_negative(u))
            {
                break; // |z| >= 1, reduced
            }

            // z -> -1/z = (-zr + i zi) / |z|^2
            arb_neg(zr, zr);
            arb_div(zr, zr, t, prec);
            arb_div(zi, zi, t, prec);

            fmpz_zero(fmpz_mat_entry(S, 0, 0));
            fmpz_one(fmpz_mat_entry(S, 0, 1));
            fmpz_set_si(fmpz_mat_entry(S, 1, 0), -1);
            fmpz_zero(fmpz_mat_entry(S, 1, 1));
            fmpz_mat_mul(M, M, S);
        }

        _acb_vec_clear(roots, 3);

        if (has_real)
        {
            ret = 1;
            break;
        }
    }

    if (ret)
    {
        bkf_apply_mat(C, M);
    }
    else
    {
        fmpz_mat_one(M);
    }

    fmpz_mat_clear(S);
    arb_clear(xH); arb_clear(phi); arb_clear(cov1); arb_clear(cov2);
    arb_clear(zr); arb_clear(zi); arb_clear(t); arb_clear(u);
    fmpz_poly_clear(cub);
    fmpz_clear(I); fmpz_clear(J); fmpz_clear(H); fmpz_clear(R);
    fmpz_clear(disc); fmpz_clear(tmp); fmpz_clear(n);

    return ret;
}
