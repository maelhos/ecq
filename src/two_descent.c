#include "two_descent.h"

void descent_eq_init(descent_eq_t eq)
{
    qfb_init(eq->q1);
    qfb_init(eq->q2);
    qfb_init(eq->q3);

    fmpz_init(eq->b1);
    fmpz_init(eq->b2);
    fmpz_init(eq->b3);

    bkf_init(eq->sieve_eq);
}

void descent_eq_clear(descent_eq_t eq)
{
    qfb_clear(eq->q1);
    qfb_clear(eq->q2);
    qfb_clear(eq->q3);

    fmpz_clear(eq->b1);
    fmpz_clear(eq->b2);
    fmpz_clear(eq->b3);

    bkf_clear(eq->sieve_eq);
}

static inline void compute_sieving_eq(bkf_t F, const qfb_t q1, const qfb_t q2, 
                                      const fmpz_t b1, const fmpz_t b2, const fmpz_t r21)
{
    // = b1*(a*u^2 + b*u*v + c*v^2)^2 - b2*(ap*u^2 + bp*u*v + cp*v^2)
    // = (a^2*b1 - ap^2*b2)*u^4 + 
    //   (2*a*b*b1 - 2*ap*bp*b2)*u^3*v + 
    //   (b^2*b1 + 2*a*c*b1 - bp^2*b2 - 2*ap*cp*b2)*u^2*v^2 + 
    //   (2*b*c*b1 - 2*bp*cp*b2)*u*v^3 + 
    //   (c^2*b1 - cp^2*b2)*v^4
    
    fmpz_t a_b1, ap_b2, c_b1, cp_b2, tmp;
    
    fmpz_init(a_b1);
    fmpz_init(ap_b2);
    fmpz_init(c_b1);
    fmpz_init(cp_b2);
    fmpz_init(tmp);

    fmpz_mul(a_b1, q1->a, b1);
    fmpz_mul(ap_b2, q2->a, b2);
    fmpz_mul(c_b1, q1->c, b1);
    fmpz_mul(cp_b2, q2->c, b2);

    // (a^2*b1 - ap^2*b2)*u^4
    fmpz_mul(F->a, a_b1, q1->a);
    fmpz_mul(tmp, ap_b2, q2->a);
    fmpz_sub(F->a, F->a, tmp);

    // 2(a*b*b1 - ap*bp*b2)*u^3*v
    fmpz_mul(F->b, a_b1, q1->b);
    fmpz_mul(tmp, ap_b2, q2->b);
    fmpz_sub(F->b, F->b, tmp);
    fmpz_add(F->b, F->b, F->b);

    // (b^2*b1 + 2*a*c*b1 - bp^2*b2 - 2*ap*cp*b2)*u^2*v^2
    fmpz_mul(F->c, q1->b, q1->b);
    fmpz_mul(F->c, F->c, b1);
    fmpz_mul(tmp, a_b1, q1->c);
    fmpz_add(tmp, tmp, tmp);
    fmpz_add(F->c, F->c, tmp);
    fmpz_mul(tmp, q2->b, q2->b);
    fmpz_mul(tmp, tmp, b2);
    fmpz_sub(F->c, F->c, tmp);
    fmpz_mul(tmp, ap_b2, q2->c);
    fmpz_add(tmp, tmp, tmp);
    fmpz_sub(F->c, F->c, tmp);

    // 2(b*c*b1 - bp*cp*b2)*u*v^3
    fmpz_mul(F->d, c_b1, q1->b);
    fmpz_mul(tmp, cp_b2, q2->b);
    fmpz_sub(F->d, F->d, tmp);
    fmpz_add(F->d, F->d, F->d);

    // (c^2*b1 - cp^2*b2)*v^4
    fmpz_mul(F->e, c_b1, q1->c);
    fmpz_mul(tmp, cp_b2, q2->c);
    fmpz_sub(F->e, F->e, tmp);

    // divide all by r12
    fmpz_divexact(F->a, F->a, r21);
    fmpz_divexact(F->b, F->b, r21);
    fmpz_divexact(F->c, F->c, r21);
    fmpz_divexact(F->d, F->d, r21);
    fmpz_divexact(F->e, F->e, r21);

    fmpz_clear(a_b1);
    fmpz_clear(ap_b2);
    fmpz_clear(c_b1);
    fmpz_clear(cp_b2);
    fmpz_clear(tmp);
}


int two_descent(const fmpz_t e1, const fmpz_t e2, const fmpz_t e3)
{
    slong i, j, k;

    // set up full 2-torsion curve
    // y^2 = (x - e1)(x - e2)(x - e3);

    // compute divisors
    fmpz_t r21, r13, r32;
    fmpz_init(r21);
    fmpz_init(r13);
    fmpz_init(r32);

    fmpz_sub(r21, e2, e1);
    fmpz_sub(r13, e1, e3);
    fmpz_sub(r32, e3, e2);

    // factor divisors
    fmpz_factor_t r21_f, r13_f, r32_f;
    fmpz_factor_init(r21_f);
    fmpz_factor_init(r13_f);
    fmpz_factor_init(r32_f);

    fmpz_factor(r21_f, r21); fmpz_factor_sort(r21_f);
    fmpz_factor(r13_f, r13); fmpz_factor_sort(r13_f);
    fmpz_factor(r32_f, r32); fmpz_factor_sort(r32_f);

    flint_printf("e2 - e1 = "); fmpz_factor_print(r21_f); flint_printf("\n");
    flint_printf("e1 - e3 = "); fmpz_factor_print(r13_f); flint_printf("\n");
    flint_printf("e3 - e2 = "); fmpz_factor_print(r32_f); flint_printf("\n");

    // compute "reduced" discriminant factorization (no square)
    fmpz_factor_t delta_f, tmp_f;
    fmpz_t tmp, delta;
    fmpz_factor_init(delta_f);
    fmpz_factor_init(tmp_f);
    fmpz_init(tmp);
    fmpz_init(delta);

    _fmpz_factor_mul(tmp_f, r21_f, r13_f);
    _fmpz_factor_mul(delta_f, tmp_f, r32_f);
    delta_f->sign = 1;
    fmpz_factor_expand(delta, delta_f);

    flint_printf("delta = "); fmpz_factor_print(delta_f); flint_printf("\n");

    // compute square-free divisors
    ulong sqf_div_num = 1ULL << (delta_f->num + 1); 
    fmpz_factor_struct* sqf_divisors_delta;
    sqf_divisors_delta = flint_malloc(sqf_div_num * sizeof(fmpz_factor_struct));

    flint_printf("S = {");
    for (i = 0; i < sqf_div_num; i++)
    {
        fmpz_factor_init(sqf_divisors_delta + i);
        _fmpz_factor_fit_length(sqf_divisors_delta + i, __builtin_popcountl(i >> 1));
        sqf_divisors_delta[i].sign = (i & 1) ? -1 : 1;
        sqf_divisors_delta[i].num = 0;
        for (j = 1; j < delta_f->num + 1; j++)
        {
            if (i & (1ULL << j))
            {
                fmpz_set(sqf_divisors_delta[i].p + sqf_divisors_delta[i].num, delta_f->p + j - 1);
                sqf_divisors_delta[i].exp[sqf_divisors_delta[i].num] = 1;
                sqf_divisors_delta[i].num++;
            }
        }
        fmpz_factor_expand(tmp, sqf_divisors_delta + i);
        fmpz_print(tmp);
        flint_printf(" ");
    }
    flint_printf("}\n Starting main descent ...\n");

    // main descent
    fmpz_factor_struct *b1_f, *b2_f;
    fmpz_t b1, b2, b3;
    fmpz_factor_t b3_f, a_f, b_f, c_f, ap_f, bp_f, cp_f;
    fmpz_t tr[3];
    fmpz_tqf_t F, Fr;
    fmpz_t z1, z2, z3, x0p, y0p, z0p, test;
    descent_eq *cover, *curr_cover;
    slong cover_size, cover_alloc;

    cover = flint_malloc(sizeof(descent_eq) * 1);

    curr_cover = cover;
    descent_eq_init(curr_cover);

    cover_alloc = 1;
    cover_size = 0;

    fmpz_factor_init(b3_f);
    fmpz_factor_init(a_f);
    fmpz_factor_init(b_f);
    fmpz_factor_init(c_f);

    fmpz_factor_init(ap_f);
    fmpz_factor_init(bp_f);
    fmpz_factor_init(cp_f);

    fmpz_init(tr[0]); 
    fmpz_init(tr[1]); 
    fmpz_init(tr[2]);

    fmpz_tqf_init(F);
    fmpz_tqf_init(Fr);

    fmpz_init(b1);
    fmpz_init(b2);
    fmpz_init(b3);

    fmpz_init(z1);
    fmpz_init(z2);
    fmpz_init(z3);
    fmpz_init(x0p);
    fmpz_init(y0p);
    fmpz_init(z0p);

    fmpz_init(test);

    int err;

    for (i = 0; i < sqf_div_num; i++)
    {
        for (j = 0; j < sqf_div_num; j++)
        {
            b1_f = sqf_divisors_delta + i;
            b2_f = sqf_divisors_delta + j;
            _fmpz_factor_mul_square_free(b3_f, b1_f, b2_f);
            b3_f->sign = b1_f->sign * b2_f->sign;

            _fmpz_factor_mul(a_f, b1_f, r32_f);
            _fmpz_factor_mul(b_f, b2_f, r13_f);
            _fmpz_factor_mul(c_f, b3_f, r21_f);

            // check R-solubility
            if ((a_f->sign > 0 && b_f->sign > 0 && c_f->sign > 0) ||
                (a_f->sign < 0 && b_f->sign < 0 && c_f->sign < 0))
            {
                continue;
            }

            // set form
            fmpz_tqf_set(F, a_f, b_f, c_f);

            // reduce form
            fmpz_tqf_reduce(Fr, F, tr);

            // try to solve
            err = fmpz_tqf_solve_reduced(Fr, x0p, y0p, z0p);

            if (err != 1 || fmpz_is_zero(x0p) || fmpz_is_zero(y0p) || fmpz_is_zero(z0p))
            {
                continue;
            }

            fmpz_factor_expand(b1, b1_f);
            fmpz_factor_expand(b2, b2_f);
            fmpz_factor_expand(b3, b3_f);

            // compute the actual base solution
            fmpz_mul(z1, x0p, tr[0]);
            fmpz_mul(z2, y0p, tr[1]);
            fmpz_mul(z3, z0p, tr[2]);

            // compute the first base test
            fmpz_mul(test, z1,   z1);
            fmpz_mul(test, test, b1);
            fmpz_mul(tmp,  z2,   z2);
            fmpz_mul(tmp,  tmp,  b2);
            fmpz_sub(test, test, tmp);
            fmpz_divexact(test, test, r21);

            // do the actual test
            if (fmpz_is_square(test) && !fmpz_is_zero(test))
            {
                goto DESCENT_SUCCESS;
            }

            // parametrize
            fmpz_tqf_parametrize(curr_cover->q1, curr_cover->q2, curr_cover->q3,
                                 Fr, x0p, y0p, z0p);

            // scale to get unreduced solution space
            fmpz_mul(curr_cover->q1->a, curr_cover->q1->a, tr[0]);
            fmpz_mul(curr_cover->q1->b, curr_cover->q1->b, tr[0]);
            fmpz_mul(curr_cover->q1->c, curr_cover->q1->c, tr[0]);

            fmpz_mul(curr_cover->q2->a, curr_cover->q2->a, tr[1]);
            fmpz_mul(curr_cover->q2->b, curr_cover->q2->b, tr[1]);
            fmpz_mul(curr_cover->q2->c, curr_cover->q2->c, tr[1]);

            fmpz_mul(curr_cover->q3->a, curr_cover->q3->a, tr[2]);
            fmpz_mul(curr_cover->q3->b, curr_cover->q3->b, tr[2]);
            fmpz_mul(curr_cover->q3->c, curr_cover->q3->c, tr[2]);

            // copy and store bs
            fmpz_set(curr_cover->b1, b1);
            fmpz_set(curr_cover->b2, b2);
            fmpz_set(curr_cover->b3, b3);

            // compute sieving equation
            compute_sieving_eq(curr_cover->sieve_eq, curr_cover->q1, curr_cover->q2, b1, b2, r21);
            
            // check solubility and increase size if ok
            if (bkf_local_solubility(curr_cover->sieve_eq, delta, delta_f))
            {
                // reduce it 
                bkf_reduce_naive(curr_cover->sieve_eq);

                // check if we already have it (TODO: a better equivalence of quartics...)
                for (k = 0; k < cover_size; k++)
                {
                    if (bkf_eq(curr_cover->sieve_eq, (cover + k)->sieve_eq))
                    {
                        goto DUPLICATE_CURVE;
                    }
                }

                cover_size++;
                // compute the parametrization and store the binary quartic
                if (cover_size == cover_alloc)
                {
                    cover_alloc *= 2;
                    cover = flint_realloc(cover, sizeof(descent_eq) * cover_alloc);
                }

                // shift current equation and initialize it
                curr_cover = cover + cover_size;
                descent_eq_init(curr_cover);
            }
            DUPLICATE_CURVE:;
        }
    }

    flint_printf("Get %ld potential covering curves\n", cover_size);

    // check potential easy solvable
    // y^2 = ax^4 + bx^3z + cx^2z^2 + dxz^3 + ez^4
    // if a is a square than (1, 0) works (x = 0 gives a torsion point so we dont care about e being a square)
    for (i = 0; i < cover_size; i++)
    {
        curr_cover = cover + i;

        if (fmpz_is_square(curr_cover->sieve_eq->a) && !fmpz_is_zero(curr_cover->sieve_eq->a))
        {
            fmpz_set(z1, curr_cover->q1->a);
            fmpz_set(z2, curr_cover->q2->a);
            fmpz_set(z3, curr_cover->q3->a);

            fmpz_set(b1, curr_cover->b1);
            fmpz_set(b2, curr_cover->b2);
            fmpz_set(b3, curr_cover->b3);

            goto DESCENT_SUCCESS;
        }
    }

    // final descent thing
    fmpz_t B, U, V, U4, V4, U2V2, U1V3, U3V1, T_P, T_I;
    int ret;

    fmpz_init(B);
    fmpz_init(U);
    fmpz_init(V);

    fmpz_init(U4);
    fmpz_init(V4);
    fmpz_init(U2V2);
    fmpz_init(U1V3);
    fmpz_init(U3V1);

    fmpz_init(T_P);
    fmpz_init(T_I);

    fmpz_set_ui(B, 1);

    while (1)
    {
        // maybe remove u = 0
        for (fmpz_set_ui(U, 0); fmpz_cmp(U, B) < 0; fmpz_add_ui(U, U, 1))
        {
            fmpz_sub(V, B, U);
            fmpz_gcd(test, U, V);

            if (!fmpz_is_one(test))
            {
                continue;
            }

            // precompute powers of U and V
            bkf_sieve_precompute(U4, U3V1, U2V2, U1V3, V4, U, V);

            for (i = 0; i < cover_size; i++)
            {
                curr_cover = cover + i;

                // evaluate the bkf at U, V
                ret = bkf_sieve_bounded(test, curr_cover->sieve_eq, U4, U3V1, U2V2, U1V3, V4);
                if (ret == 1)
                {
                    qbf_eval(z1, curr_cover->q1, U, V);
                    qbf_eval(z2, curr_cover->q2, U, V);
                    qbf_eval(z3, curr_cover->q3, U, V);

                    fmpz_set(b1, curr_cover->b1);
                    fmpz_set(b2, curr_cover->b2);
                    fmpz_set(b3, curr_cover->b3);

                    goto DESCENT_SUCCESS;
                }
                
                if (ret == -1)
                {
                    fmpz_neg(V, V);

                    qbf_eval(z1, curr_cover->q1, U, V);
                    qbf_eval(z2, curr_cover->q2, U, V);
                    qbf_eval(z3, curr_cover->q3, U, V);

                    fmpz_set(b1, curr_cover->b1);
                    fmpz_set(b2, curr_cover->b2);
                    fmpz_set(b3, curr_cover->b3);

                    goto DESCENT_SUCCESS;
                }
            }
        }

        fmpz_add_ui(B, B, 1);
    }


    DESCENT_SUCCESS:;

    fmpz_t X_n, Y_n, Z;

    fmpz_init(X_n);
    fmpz_init(Y_n);
    fmpz_init(Z);

    fmpz_sqrt(Z, test);

    fmpz_mul(X_n, z1, z1);
    fmpz_mul(X_n, X_n, b1);
    fmpz_mul(tmp, e1, test); // test holds Z**2
    fmpz_add(X_n, X_n, tmp);

    fmpz_mul(Y_n, b1, b2);
    fmpz_mul(Y_n, Y_n, b3);
    fmpz_sqrt(Y_n, Y_n);
    fmpz_mul(Y_n, Y_n, z1);
    fmpz_mul(Y_n, Y_n, z2);
    fmpz_mul(Y_n, Y_n, z3);

    fmpq_t x, y, tmp_q, test_q;

    fmpq_init(x);
    fmpq_init(y);
    fmpq_init(tmp_q);
    fmpq_init(test_q);

    fmpz_mul(tmp, test, Z);
    fmpq_set_fmpz_frac(x, X_n, test);
    fmpq_set_fmpz_frac(y, Y_n, tmp);
    fmpq_canonicalise(x);
    fmpq_canonicalise(y);

    // print it
    flint_printf("("); fmpq_print(x); flint_printf(" : "); fmpq_print(y); flint_printf(")\n");

    // test solution
    fmpq_sub_fmpz(test_q, x, e1);
    fmpq_sub_fmpz(tmp_q, x, e2);
    fmpq_mul(test_q, test_q, tmp_q);
    fmpq_sub_fmpz(tmp_q, x, e3);
    fmpq_mul(test_q, test_q, tmp_q);

    fmpq_mul(tmp_q, y, y);

    fmpq_sub(test_q, test_q, tmp_q);

    if (!fmpq_is_zero(test_q))
    {
        flint_printf("Elliptic curve point is wrong, aborting\n");
        flint_abort();
    }

    return 1;
}