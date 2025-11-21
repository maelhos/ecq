#pragma once
#include <flint/fmpz.h>
#include <flint/qfb.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/acb.h>

#include "factor_addition.h"

/* Represent a binary quartic form over ZZ
au^4 + bu^3v + cu^2v^2 + duv^3 + ev^4 
we dont require the discriminant to be set !
we only compute it and book-keep it when necessary */
typedef struct
{
    fmpz_t a, b, c, d, e;
    fmpz_factor_t disc;
} bkf;

typedef bkf bkf_t[1];

static inline void bkf_init(bkf_t C)
{
    fmpz_init(C->a);
    fmpz_init(C->b);
    fmpz_init(C->c);
    fmpz_init(C->d);
    fmpz_init(C->e);
}

static inline void bkf_set(bkf_t C, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_t e)
{
    fmpz_set(C->a, a);
    fmpz_set(C->b, b);
    fmpz_set(C->c, c);
    fmpz_set(C->d, d);
    fmpz_set(C->e, e);
}

static inline void bkf_init_set(bkf_t C, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d, const fmpz_t e)
{
    bkf_init(C);
    bkf_set(C, a, b, c, d, e);
}

static inline void bkf_clear(bkf_t C)
{
    fmpz_clear(C->a);
    fmpz_clear(C->b);
    fmpz_clear(C->c);
    fmpz_clear(C->d);
    fmpz_clear(C->e);
}

static inline void bkf_print(const bkf_t C)
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

static inline int bkf_eq(const bkf_t C1, const bkf_t C2)
{
    if (fmpz_cmp(C1->a, C2->a) != 0) return 0;
    if (fmpz_cmp(C1->b, C2->b) != 0) return 0;
    if (fmpz_cmp(C1->c, C2->c) != 0) return 0;
    if (fmpz_cmp(C1->d, C2->d) != 0) return 0;
    if (fmpz_cmp(C1->e, C2->e) != 0) return 0;

    return 1;
}

static inline void qbf_eval(fmpz_t r, const qfb_t C, const fmpz_t u, const fmpz_t v)
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

void bkf_sieve_precompute(fmpz_t U4, fmpz_t U3V1, fmpz_t U2V2, fmpz_t U1V3, fmpz_t V4, const fmpz_t U, const fmpz_t V);
int bkf_sieve_bounded(fmpz_t test, const bkf_t C, const fmpz_t U4, const fmpz_t U3V1, const fmpz_t U2V2, const fmpz_t U1V3, const fmpz_t V4);

int bkf_local_solubility(bkf_t C, const fmpz_t delta_curve, const fmpz_factor_t delta_curve_f);
int bkf_real_solubility(const bkf_t C);
int bkf_Qp_solubility(const bkf_t C, const fmpz_t p);
int bkf_Zp_solubility(const bkf_t C, const fmpz_t p, const fmpz_t x, ulong k);
void bkf_disc(fmpz_t disc, const bkf_t C);

void bkf_reduce_naive(bkf_t C);