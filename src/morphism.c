#include "morphism.h"

// C(t, 1) = a t^4 + b t^3 + c t^2 + d t + e
static void bkf_eval_fmpq(fmpq_t res, const bkf_t C, const fmpq_t t)
{
    fmpq_set_fmpz(res, C->a);
    fmpq_mul(res, res, t);
    fmpq_add_fmpz(res, res, C->b);
    fmpq_mul(res, res, t);
    fmpq_add_fmpz(res, res, C->c);
    fmpq_mul(res, res, t);
    fmpq_add_fmpz(res, res, C->d);
    fmpq_mul(res, res, t);
    fmpq_add_fmpz(res, res, C->e);
}

void morphx(fmpq_t x, const morph_data_t morph, const fmpq_t u, const fmpq_t v, morph_type mtype)
{
    switch (mtype)
    {
    case MORPHISM_2ISOGENY:
        fmpq_mul_fmpz(x, u, morph->d1);
        fmpq_mul(x, x, u);
        fmpq_add_fmpz(x, x, morph->x0);
        break;

    case MORPHISM_ISOGENOUS_2ISOGENY:
        fmpq_div(x, v, u);
        fmpq_div_2exp(x, x, 1);
        fmpq_mul(x, x, x);
        fmpq_add_fmpz(x, x, morph->x0);
        break;

    case MORPHISM_SECOND_DESCENT:
    {
        // x = d1 U^2 + x0 = d1 Q1/Q3 + x0
        fmpq_t r1, r3;

        fmpq_init(r1);
        fmpq_init(r3);

        bkf_eval_fmpq(r1, morph->Q1, u);
        bkf_eval_fmpq(r3, morph->Q3, u);

        if (fmpq_is_zero(r3))
        {
            flint_printf("morphx: Q3 vanishes at the point, degenerate descendant\n");
            flint_abort();
        }

        fmpq_div(x, r1, r3);
        fmpq_mul_fmpz(x, x, morph->d1);
        fmpq_add_fmpz(x, x, morph->x0);

        fmpq_clear(r1);
        fmpq_clear(r3);
        break;
    }

    case MORPHISM_ISOGENOUS_SECOND_DESCENT:
    {
        // x = V^2/(4 U^2) + x0 = Q2^2 / (4 Q1 Q3) + x0
        fmpq_t r1, r2, r3;

        fmpq_init(r1);
        fmpq_init(r2);
        fmpq_init(r3);

        bkf_eval_fmpq(r1, morph->Q1, u);
        bkf_eval_fmpq(r2, morph->Q2, u);
        bkf_eval_fmpq(r3, morph->Q3, u);

        fmpq_mul(r1, r1, r3);

        if (fmpq_is_zero(r1))
        {
            flint_printf("morphx: Q1 Q3 vanishes at the point, degenerate descendant\n");
            flint_abort();
        }

        fmpq_mul(x, r2, r2);
        fmpq_div(x, x, r1);
        fmpq_div_2exp(x, x, 2);
        fmpq_add_fmpz(x, x, morph->x0);

        fmpq_clear(r1);
        fmpq_clear(r2);
        fmpq_clear(r3);
        break;
    }

    default:
        flint_printf("morphx: unknown morphism type %d\n", (int) mtype);
        flint_abort();
        break;
    }
}

void morph_data_clear(morph_data_t morph, morph_type mtype)
{
    fmpz_clear(morph->x0);
    fmpz_clear(morph->d1);
    fmpz_clear(morph->d3);

    if (mtype == MORPHISM_SECOND_DESCENT || mtype == MORPHISM_ISOGENOUS_SECOND_DESCENT)
    {
        bkf_clear(morph->Q1);
        bkf_clear(morph->Q2);
        bkf_clear(morph->Q3);
    }
}
