#pragma once


#include <flint/fmpz.h>
#include <flint/qfb.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/acb.h>
#include <math.h>

// Maximum refinement depth for the p-adic solubility recursion (bkf_Zp_solubility).
#define BKF_ZP_MAX_DEPTH 10000

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


static inline void bkf_set_2iso(bkf_t C, const fmpz_t d1, const fmpz_t c, const fmpz_t d2)
{
    fmpz_set(C->a, d1);
    fmpz_zero(C->b);
    fmpz_set(C->c, c);
    fmpz_zero(C->d);
    fmpz_set(C->e, d2);
}

static inline void bkf_init_set_2iso(bkf_t C, const fmpz_t d1, const fmpz_t c, const fmpz_t d2)
{
    bkf_init(C);
    bkf_set_2iso(C, d1, c, d2);
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

static inline void qbf_eval(fmpz_t r, const qfb_t C, long u, long v)
{
    // aU^2 + bUV + cV^2
    // U(aU + bV) + cV^2;

    fmpz_t tmp;

    fmpz_init(tmp);

    fmpz_mul_si(r, C->a, u);
    fmpz_mul_si(tmp, C->b, v);
    fmpz_add(r, r, tmp);
    fmpz_mul_si(r, r, u);
    fmpz_ui_mul_ui(tmp, v, v);
    fmpz_mul(tmp, tmp, C->c);
    fmpz_add(r, r, tmp);

    fmpz_clear(tmp);
}

void bkf_power_precompute(fmpz* pows, long U, long V);
int bkf_power_bounded(fmpz_t test, const bkf *C, const fmpz* pows);
int bkf_local_solubility(bkf_t C, const fmpz_t delta_curve, const fmpz_factor_t delta_curve_f);
int bkf_everywhere_local_solubility_2isogeny(const bkf_t C, const fmpz_t c, const fmpz_t d, const fmpz_t dp, const fmpz_t d1, const fmpz_factor_t p_list);

int bkf_real_solubility(const bkf_t C);
int bkf_Qp_solubility(const bkf_t C, const fmpz_t p);
int bkf_Zp_solubility(const bkf_t C, const fmpz_t p, const fmpz_t x, ulong k);
void bkf_disc(fmpz_t disc, const bkf_t C);

void bkf_reduce_naive(bkf_t C);
void bkf_mul2qbf(bkf_t C, const qfb_t q1, const qfb_t q2);

// C(u, v) = q(P1(u, v), P3(u, v)), the composition of a binary quadratic with two
// binary quadratics, which is a binary quartic.
void bkf_compose_qfb(bkf_t C, const qfb_t q, const qfb_t P1, const qfb_t P3);

// Magma / Sage friendly form of the quartic y^2 = a x^4 + b x^3 + c x^2 + d x + e
static inline void bkf_print_magma(const bkf_t C)
{
    flint_printf("[%{fmpz}, %{fmpz}, %{fmpz}, %{fmpz}, %{fmpz}]", C->a, C->b, C->c, C->d, C->e);
}

// C <- k * C
static inline void bkf_scalar_mul(bkf_t C, const fmpz_t k)
{
    fmpz_mul(C->a, C->a, k);
    fmpz_mul(C->b, C->b, k);
    fmpz_mul(C->c, C->c, k);
    fmpz_mul(C->d, C->d, k);
    fmpz_mul(C->e, C->e, k);
}