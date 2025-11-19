#pragma once
#include <flint/fmpz.h>
#include <flint/qfb.h>

/* Represent a binary quartic form over ZZ
au^4 + bu^3v + cu^2v^2 + duv^3 + ev^4 */
typedef struct
{
    fmpz_t a, b, c, d, e;
} fmpz_bkf;

typedef fmpz_bkf fmpz_bkf_t[1];

void qbf_eval(fmpz_t r, const qfb_t C, const fmpz_t u, const fmpz_t v);

void fmpz_bkf_init(fmpz_bkf_t C);
void fmpz_bkf_set(fmpz_bkf_t C, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, fmpz_t e);
void fmpz_bkf_init_set(fmpz_bkf_t C, fmpz_t a, fmpz_t b, fmpz_t c, fmpz_t d, fmpz_t e);
void fmpz_bkf_clear(fmpz_bkf_t C);
void fmpz_bkf_print(const fmpz_bkf_t C);