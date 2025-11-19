#include "bkf.h"

void fmpz_bkf_init(fmpz_bkf_t C)
{
    fmpz_init(C->a);
    fmpz_init(C->b);
    fmpz_init(C->c);
    fmpz_init(C->d);
    fmpz_init(C->e);
}

void fmpz_bkf_set(fmpz_bkf_t C, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, fmpz_t e)
{
    fmpz_set(C->a, a);
    fmpz_set(C->b, b);
    fmpz_set(C->c, c);
    fmpz_set(C->d, d);
    fmpz_set(C->e, e);
}

void fmpz_bkf_init_set(fmpz_bkf_t C, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, fmpz_t e)
{
    fmpz_bkf_init(C);
    fmpz_bkf_set(C, a, b, c, d, e);
}

void fmpz_bkf_clear(fmpz_bkf_t C)
{
    fmpz_clear(C->a);
    fmpz_clear(C->b);
    fmpz_clear(C->c);
    fmpz_clear(C->d);
    fmpz_clear(C->e);
}

void fmpz_bkf_print(const fmpz_bkf_t C)
{
    fmpz_print(C->a);
    flint_printf("u^4 + ");
    fmpz_print(C->b);
    flint_printf("u^3v + ");
    fmpz_print(C->c);
    flint_printf("u^2v^2 + ");
    fmpz_print(C->d);
    flint_printf("uv^3 + ");
    fmpz_print(C->e);
    flint_printf("v^4");
}

void qbf_eval(fmpz_t r, const qfb_t C, const fmpz_t u, const fmpz_t v)
{
    // aU^2 + bUV + cV^2
    // U(aU + bV) + cV^2;

    fmpz_t tmp;

    fmpz_init(tmp);

    fmpz_mul(r, C->a, u);
    fmpz_mul(tmp, C->b, v);
    fmpz_add(r, r, tmp);
    fmpz_mul(r, r, u);
    fmpz_mul(tmp, v, v);
    fmpz_mul(tmp, tmp, C->c);
    fmpz_add(r, r, tmp);

    fmpz_clear(tmp);
}