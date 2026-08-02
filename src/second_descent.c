#include "second_descent.h"
#include "minimise.h"

// Everywhere local solubility of a descendant quartic: real, plus Q_p for every
// prime that can be bad (those dividing the curve data, collected in p_all).
static int descendant_els(const bkf_t G, const fmpz_factor_t p_all)
{
    slong i;

    if (!bkf_real_solubility(G))
    {
        return 0;
    }

    for (i = 0; i < p_all->num; i++)
    {
        if (!bkf_Qp_solubility(G, p_all->p + i))
        {
            return 0;
        }
    }

    return 1;
}

/* Solve the conic S^2 = D1 a^2 + C1 a b + D2 b^2 and parametrize its solutions.
   Mirrors the first descent path in main.c: build the diagonal ternary form
   (-1, 4 D1, delta), reduce, solve, shift back to the semi-diagonal solution.
   Returns 1 on success. */
static int solve_semi_diag_conic(qfb_t P1, qfb_t P2, qfb_t P3,
                                 const fmpz_t D1, const fmpz_t C1, const fmpz_t D2,
                                 const fmpz_t delta)
{
    fmpz_tqf_t tqf, tqf_red;
    fmpz_factor_t fa, fb, fc;
    fmpz_t T[3], z1, z2, z3, tmp, tmp2;
    int ret = 0;

    if (fmpz_is_zero(D1) || fmpz_is_zero(delta))
    {
        return 0;
    }

    // The ternary form is (-1, 4 D1, delta). With all three coefficients of the same
    // sign it has no real point at all, so it is insoluble and must not be handed to
    // fmpz_tqf_solve_reduced (which aborts when its certificate disagrees).
    if (fmpz_sgn(D1) < 0 && fmpz_sgn(delta) < 0)
    {
        return 0;
    }

    fmpz_tqf_init(tqf);
    fmpz_tqf_init(tqf_red);
    fmpz_factor_init(fa);
    fmpz_factor_init(fb);
    fmpz_factor_init(fc);
    fmpz_init(T[0]); fmpz_init(T[1]); fmpz_init(T[2]);
    fmpz_init(z1); fmpz_init(z2); fmpz_init(z3);
    fmpz_init(tmp); fmpz_init(tmp2);

    // coefficient 0 is -1
    fmpz_factor_set_one(fa);
    fa->sign = -1;

    // coefficient 1 is 4 D1
    fmpz_mul_ui(tmp, D1, 4);
    fmpz_factor(fb, tmp);

    // coefficient 2 is delta = C1^2 - 4 D1 D2
    fmpz_factor(fc, delta);

    fmpz_tqf_set(tqf, fa, fb, fc);
    fmpz_tqf_reduce(tqf_red, tqf, T);

    if (fmpz_tqf_solve_reduced(tqf_red, z1, z2, z3) != 1)
    {
        goto cleanup;
    }

    // undo the reduction shift
    fmpz_mul(z1, z1, T[0]);
    fmpz_mul(z2, z2, T[1]);
    fmpz_mul(z3, z3, T[2]);

    // diagonal -> semi diagonal
    fmpz_submul(z1, C1, z3);
    fmpz_mul(z2, z2, D1);
    fmpz_mul_2exp(z2, z2, 1);
    fmpz_mul(z3, z3, D1);
    fmpz_mul_2exp(z3, z3, 1);

    fmpz_gcd3(tmp, z1, z2, z3);
    if (fmpz_is_zero(tmp))
    {
        goto cleanup;
    }
    fmpz_divexact(z1, z1, tmp);
    fmpz_divexact(z2, z2, tmp);
    fmpz_divexact(z3, z3, tmp);

    // z2^2 == D1 z1^2 + C1 z1 z3 + D2 z3^2 ?
    fmpz_mul(tmp, z1, z1);
    fmpz_mul(tmp, tmp, D1);
    fmpz_mul(tmp2, z1, z3);
    fmpz_mul(tmp2, tmp2, C1);
    fmpz_add(tmp, tmp, tmp2);
    fmpz_mul(tmp2, z3, z3);
    fmpz_mul(tmp2, tmp2, D2);
    fmpz_add(tmp, tmp, tmp2);
    fmpz_mul(tmp2, z2, z2);

    if (!fmpz_equal(tmp, tmp2))
    {
        goto cleanup;
    }

    fmpz_tqf_parametrize_semi_diag(P1, P2, P3, D1, C1, D2, delta, z1, z2, z3);
    ret = 1;

cleanup:
    fmpz_tqf_clear(tqf);
    fmpz_tqf_clear(tqf_red);
    fmpz_factor_clear(fa);
    fmpz_factor_clear(fb);
    fmpz_factor_clear(fc);
    fmpz_clear(T[0]); fmpz_clear(T[1]); fmpz_clear(T[2]);
    fmpz_clear(z1); fmpz_clear(z2); fmpz_clear(z3);
    fmpz_clear(tmp); fmpz_clear(tmp2);
    return ret;
}

slong second_descent_d3(two_cover_t cover,
                        const fmpz_t d1, const fmpz_t c, const fmpz_t d2, const fmpz_t x0,
                        const qfb_t q1, const qfb_t q2, const qfb_t q3,
                        const fmpz_factor_t resf,
                        const fmpz_factor_t p_all,
                        morph_type mtype,
                        slong max_covers)
{
    fmpz_t d3, D1, C1, D2, delta, tmp;
    fmpz* primes;
    qfb_t P1, P2, P3;
    bkf_t Q1, Q2, Q3, G;
    fmpz_mat_t M;
    slong nprimes, i, k, added;
    ulong j;
    int sgn, has2;
    bkf_t *best_G, *best_Q1, *best_Q2, *best_Q3;
    fmpz *best_d3;
    slong *best_size;

    fmpz_init(d3); fmpz_init(D1); fmpz_init(C1); fmpz_init(D2);
    fmpz_init(delta); fmpz_init(tmp);
    qfb_init(P1); qfb_init(P2); qfb_init(P3);
    bkf_init(Q1); bkf_init(Q2); bkf_init(Q3); bkf_init(G);
    fmpz_mat_init(M, 2, 2);

    // d3 runs over squarefree divisors of 2 * res, with both signs
    has2 = 0;
    for (i = 0; i < resf->num; i++)
    {
        if (fmpz_cmp_ui(resf->p + i, 2) == 0) has2 = 1;
    }

    nprimes = resf->num + (has2 ? 0 : 1);
    primes = _fmpz_vec_init(nprimes);
    for (i = 0; i < resf->num; i++)
    {
        fmpz_set(primes + i, resf->p + i);
    }
    if (!has2)
    {
        fmpz_set_ui(primes + resf->num, 2);
    }

    best_G  = flint_malloc(max_covers * sizeof(bkf_t));
    best_Q1 = flint_malloc(max_covers * sizeof(bkf_t));
    best_Q2 = flint_malloc(max_covers * sizeof(bkf_t));
    best_Q3 = flint_malloc(max_covers * sizeof(bkf_t));
    best_size = flint_malloc(max_covers * sizeof(slong));
    best_d3 = _fmpz_vec_init(max_covers);
    for (i = 0; i < max_covers; i++)
    {
        bkf_init(best_G[i]); bkf_init(best_Q1[i]);
        bkf_init(best_Q2[i]); bkf_init(best_Q3[i]);
        best_size[i] = 0;
    }

    added = 0;

    /* All d3 giving an ELS descendant lie in one coset of H1 (Cremona d2 3.2.3), so
       they all represent the SAME class in S^2(E) and the coverings are equivalent.
       We therefore keep the max_covers smallest of them rather than the first ones
       found: same Selmer classes, but the sieve has a much easier time on a model
       with small coefficients. */
    for (sgn = 0; sgn < 2; sgn++)
    for (j = 0; j < (UWORD(1) << nprimes); j++)
    {
        fmpz_set_ui(d3, 1);
        for (k = 0; k < nprimes; k++)
        {
            if (j & (UWORD(1) << k))
            {
                fmpz_mul(d3, d3, primes + k);
            }
        }
        if (sgn)
        {
            fmpz_neg(d3, d3);
        }

        // conic (9): q1(a,b) = d3 s^2, cleared of denominators as
        // (d3 s)^2 = d3 q1(a, b), i.e. S^2 = D1 a^2 + C1 a b + D2 b^2
        fmpz_mul(D1, d3, q1->a);
        fmpz_mul(C1, d3, q1->b);
        fmpz_mul(D2, d3, q1->c);

        fmpz_mul(delta, C1, C1);
        fmpz_mul(tmp, D1, D2);
        fmpz_submul_ui(delta, tmp, 4);

        if (!solve_semi_diag_conic(P1, P2, P3, D1, C1, D2, delta))
        {
            continue;
        }

        // (a, b) = (P1, P3); compose the conic forms with it
        bkf_compose_qfb(Q1, q1, P1, P3);
        bkf_compose_qfb(Q2, q2, P1, P3);
        bkf_compose_qfb(Q3, q3, P1, P3);

        // descendant (10): q3(a,b) = d3 t^2, searched as d3 * Q3 = (d3 t)^2
        bkf_set(G, Q3->a, Q3->b, Q3->c, Q3->d, Q3->e);
        bkf_scalar_mul(G, d3);

        if (fmpz_is_zero(G->a) && fmpz_is_zero(G->b) && fmpz_is_zero(G->c)
            && fmpz_is_zero(G->d) && fmpz_is_zero(G->e))
        {
            continue;
        }

        if (!descendant_els(G, p_all))
        {
            continue;
        }

        // Simplify the descendant: divide out the square part of the content, minimise
        // (this is what actually lowers I and J, by the w^4, w^6 the second descent
        // introduces) and reduce. Every change of variable has to go through Q1, Q2, Q3
        // so that the morphism back to E stays consistent. bkf_minimise may transform G
        // even when it performs no minimisation step, so M is always applied.
        bkf_reduce_naive(G);

        bkf_minimise(G, M);
        bkf_apply_mat(Q1, M);
        bkf_apply_mat(Q2, M);
        bkf_apply_mat(Q3, M);

        if (bkf_reduce(G, M))
        {
            bkf_apply_mat(Q1, M);
            bkf_apply_mat(Q2, M);
            bkf_apply_mat(Q3, M);
        }

        // skip a descendant we already have (different d3 can reduce to the same quartic)
        {
            slong t;
            int dup = 0;
            for (t = 0; t < (slong) cover->size; t++)
            {
                bkf* e = (cover->cover + t)->eq;
                if (fmpz_equal(e->a, G->a) && fmpz_equal(e->b, G->b) && fmpz_equal(e->c, G->c)
                 && fmpz_equal(e->d, G->d) && fmpz_equal(e->e, G->e))
                {
                    dup = 1;
                    break;
                }
            }
            if (dup) continue;
        }

        // size of this candidate, as the bit length of its largest coefficient
        {
            slong sz = fmpz_bits(G->a);
            if (fmpz_bits(G->b) > sz) sz = fmpz_bits(G->b);
            if (fmpz_bits(G->c) > sz) sz = fmpz_bits(G->c);
            if (fmpz_bits(G->d) > sz) sz = fmpz_bits(G->d);
            if (fmpz_bits(G->e) > sz) sz = fmpz_bits(G->e);

            if (added < max_covers)
            {
                best_size[added] = sz;
                bkf_set(best_G[added], G->a, G->b, G->c, G->d, G->e);
                bkf_set(best_Q1[added], Q1->a, Q1->b, Q1->c, Q1->d, Q1->e);
                bkf_set(best_Q2[added], Q2->a, Q2->b, Q2->c, Q2->d, Q2->e);
                bkf_set(best_Q3[added], Q3->a, Q3->b, Q3->c, Q3->d, Q3->e);
                fmpz_set(best_d3 + added, d3);
                added++;
            }
            else
            {
                slong worst = 0, t2;
                for (t2 = 1; t2 < max_covers; t2++)
                {
                    if (best_size[t2] > best_size[worst]) worst = t2;
                }
                if (sz < best_size[worst])
                {
                    best_size[worst] = sz;
                    bkf_set(best_G[worst], G->a, G->b, G->c, G->d, G->e);
                    bkf_set(best_Q1[worst], Q1->a, Q1->b, Q1->c, Q1->d, Q1->e);
                    bkf_set(best_Q2[worst], Q2->a, Q2->b, Q2->c, Q2->d, Q2->e);
                    bkf_set(best_Q3[worst], Q3->a, Q3->b, Q3->c, Q3->d, Q3->e);
                    fmpz_set(best_d3 + worst, d3);
                }
            }
        }
    }

    for (i = 0; i < added; i++)
    {
        flint_printf("       d3 = ");
        fmpz_print(best_d3 + i);
        flint_printf("  descendant ");
        bkf_print_magma(best_G[i]);
        flint_printf("\n");

        add_cover_second_descent(cover, d1, best_d3 + i, x0, best_G[i],
                                 best_Q1[i], best_Q2[i], best_Q3[i], mtype);
    }

    for (i = 0; i < max_covers; i++)
    {
        bkf_clear(best_G[i]); bkf_clear(best_Q1[i]);
        bkf_clear(best_Q2[i]); bkf_clear(best_Q3[i]);
    }
    flint_free(best_G); flint_free(best_Q1); flint_free(best_Q2); flint_free(best_Q3);
    flint_free(best_size);
    _fmpz_vec_clear(best_d3, max_covers);

    _fmpz_vec_clear(primes, nprimes);
    fmpz_clear(d3); fmpz_clear(D1); fmpz_clear(C1); fmpz_clear(D2);
    fmpz_clear(delta); fmpz_clear(tmp);
    qfb_clear(P1); qfb_clear(P2); qfb_clear(P3);
    bkf_clear(Q1); bkf_clear(Q2); bkf_clear(Q3); bkf_clear(G);
    fmpz_mat_clear(M);

    return added;
}
