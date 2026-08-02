#include "two_descent.h"

void two_cover_init(two_cover_t eq)
{
    eq->alloc = 0;
    eq->size = 0;
    eq->cover = 0;
}

void two_cover_clear(two_cover_t eq)
{
    descent_eq* curr;

    for (slong i = 0; i < eq->size; i++)
    {
        curr = eq->cover + i;
        bkf_clear(curr->eq);
        // morphism is an inline member of descent_eq, not a separate allocation:
        // flint_free on it would be an invalid free of an interior pointer.
        morph_data_clear(curr->morphism, curr->mtype);
    }

    flint_free(eq->cover);
    eq->cover = 0;
    eq->size = 0;
    eq->alloc = 0;
}

// Drop covers whose quartic is literally repeated. Reduction is canonical enough that
// equal coverings coming from the same branch land on the same reduced form; genuinely
// equivalent covers reduced through a different covariant branch are NOT caught here
// (see bkf_equivalent, which is sound but incomplete).
void two_cover_dedup(two_cover_t cover)
{
    slong i, j, k = 0;

    for (i = 0; i < (slong) cover->size; i++)
    {
        bkf* a = (cover->cover + i)->eq;
        int dup = 0;

        for (j = 0; j < k; j++)
        {
            bkf* b = (cover->cover + j)->eq;
            if (fmpz_equal(a->a, b->a) && fmpz_equal(a->b, b->b) && fmpz_equal(a->c, b->c)
             && fmpz_equal(a->d, b->d) && fmpz_equal(a->e, b->e))
            {
                dup = 1;
                break;
            }
        }

        if (dup)
        {
            bkf_clear(a);
            morph_data_clear((cover->cover + i)->morphism, (cover->cover + i)->mtype);
        }
        else
        {
            if (k != i)
            {
                cover->cover[k] = cover->cover[i];
            }
            k++;
        }
    }

    cover->size = k;
}

void two_cover_print(const two_cover_t cover)
{
    descent_eq* curr;
    slong i;

    flint_printf("\n================= 2-coverings (%wu) =================\n", cover->size);
    flint_printf("each is  y^2 = a x^4 + b x^3 + c x^2 + d x + e  given as [a, b, c, d, e]\n\n");

    for (i = 0; i < (slong) cover->size; i++)
    {
        curr = cover->cover + i;

        int second = (curr->mtype == MORPHISM_SECOND_DESCENT
                   || curr->mtype == MORPHISM_ISOGENOUS_SECOND_DESCENT);

        // flint_printf does not accept width modifiers on %{fmpz}
        flint_printf("  %2wd  %-2s  %-6s  d1 = ", i + 1,
            (curr->mtype == MORPHISM_2ISOGENY || curr->mtype == MORPHISM_SECOND_DESCENT) ? "E" : "E'",
            second ? "2nd" : "1st");
        fmpz_print(curr->morphism->d1);
        flint_printf("  d3 = ");
        fmpz_print(curr->morphism->d3);
        flint_printf("  ");
        bkf_print_magma(curr->eq);
        flint_printf("\n");
    }

    flint_printf("\nMagma:\n  covers := [\n");
    for (i = 0; i < (slong) cover->size; i++)
    {
        flint_printf("    ");
        bkf_print_magma((cover->cover + i)->eq);
        flint_printf("%s\n", (i + 1 < (slong) cover->size) ? "," : "");
    }
    flint_printf("  ];\n\n");
}


void two_descent(fmpq_t X, const two_cover_t cover)
{
    // final descent thing
    fmpz_t test, T_P, T_I;
    fmpz* pows;
    fmpq_t a, b;
    fmpz* pol;
    ulong B, U, V, nb_primes;
    descent_eq* curr_cover;
    slong i;
    int ret;

    pol = _fmpz_vec_init(5);

    fmpq_init(a);
    fmpq_init(b);

    fmpz_init(test);
    pows = _fmpz_vec_init(5);

    fmpz_init(T_P);
    fmpz_init(T_I);

    // prepare sieve assisted search
    nb_primes = 16;
    B = 1ULL << (nb_primes / 2); 
    while (1)
    {
        for (i = 0; i < cover->size; i++)
        {
            curr_cover = cover->cover + i;

            // set the thing
            fmpz_set(pol + 4, curr_cover->eq->a);
            fmpz_set(pol + 3, curr_cover->eq->b);
            fmpz_set(pol + 2, curr_cover->eq->c);
            fmpz_set(pol + 1, curr_cover->eq->d);
            fmpz_set(pol + 0, curr_cover->eq->e);

            // we avoid p = 2 :(
            if (sieve_pol(&U, &V, pol, 5, B, n_primes_arr_readonly(nb_primes + 1) + 1, nb_primes))
            {
                bkf_power_precompute(pows, U, V);
                ret = bkf_power_bounded(test, curr_cover->eq, pows);

                if (ret != 1)
                {
                    flint_printf("Sieve failed to find a square, aborting\n");
                    flint_abort();
                }
                // (a, b) = (U/V, y/V^2) on the homogeneous space 

                fmpz_sqrt(test, test);
                flint_printf("Found solution for U = %ld, V = %ld, y = %{fmpz} on ", U, V, test);
                bkf_print(curr_cover->eq);
                flint_printf("\n");

                fmpq_set_ui(a, U, V);
                fmpq_set_fmpz(b, test);
                
                fmpz_set_ui(test, V);
                fmpz_mul(test, test, test);

                fmpq_div_fmpz(b, b, test);
    
                morphx(X, curr_cover->morphism, a, b, curr_cover->mtype);

                _fmpz_vec_clear(pol, 5);
                fmpq_clear(a); fmpq_clear(b);
                fmpz_clear(test);
                _fmpz_vec_clear(pows, 5);
                fmpz_clear(T_P); fmpz_clear(T_I);
                return;
            }
        }

        B *= 2;
        nb_primes += 2; // approximation since nb_prime = 2*log2(H)
    }

}

void add_cover_2isogeny(two_cover_t cover, const fmpz_t d1, const fmpz_t c, const fmpz_t d2, const fmpz_t x0, morph_type type)
{
    if (cover->alloc == cover->size)
    {
        cover->alloc++; cover->alloc *= 2;
        cover->cover = flint_realloc(cover->cover, cover->alloc * sizeof(descent_eq));
    }

    descent_eq* curr = cover->cover + cover->size;
    bkf_init_set_2iso(curr->eq, d1, c, d2);
    curr->mtype = type;
    
    if (type == MORPHISM_ISOGENOUS_2ISOGENY)
    {
        fmpz_init_set(curr->morphism->x0, x0);
        fmpz_init_set(curr->morphism->d1, d1);
    }
    else if (type == MORPHISM_2ISOGENY)
    {
        fmpz_init_set(curr->morphism->x0, x0);
        // morphx reads d1 for this type too, so it must be initialised here
        fmpz_init_set(curr->morphism->d1, d1);
    }
    else
    {
        flint_printf("Wrong morphism type for 2-isogeny quartic\n");
        flint_abort();
    }

    fmpz_init_set_ui(curr->morphism->d3, 1);

    cover->size++;
}

void add_cover_second_descent(two_cover_t cover, const fmpz_t d1, const fmpz_t d3, const fmpz_t x0, const bkf_t q,
                        const bkf_t Q1, const bkf_t Q2, const bkf_t Q3, morph_type type)
{
    if (cover->alloc == cover->size)
    {
        cover->alloc++; cover->alloc *= 2;
        cover->cover = flint_realloc(cover->cover, cover->alloc * sizeof(descent_eq));
    }

    descent_eq* curr = cover->cover + cover->size;
    bkf_init_set(curr->eq, q->a, q->b, q->c, q->d, q->e);
    curr->mtype = type;

    bkf_init(curr->morphism->Q1);
    bkf_init(curr->morphism->Q2);
    bkf_init(curr->morphism->Q3);

    bkf_set(curr->morphism->Q1, Q1->a, Q1->b, Q1->c, Q1->d, Q1->e);
    bkf_set(curr->morphism->Q2, Q2->a, Q2->b, Q2->c, Q2->d, Q2->e);
    bkf_set(curr->morphism->Q3, Q3->a, Q3->b, Q3->c, Q3->d, Q3->e);

    if (type != MORPHISM_SECOND_DESCENT && type != MORPHISM_ISOGENOUS_SECOND_DESCENT)
    {
        flint_printf("Wrong morphism type for second descent quartic\n");
        flint_abort();
    }

    fmpz_init_set(curr->morphism->x0, x0);
    fmpz_init_set(curr->morphism->d1, d1);
    fmpz_init_set(curr->morphism->d3, d3);

    cover->size++;
}
