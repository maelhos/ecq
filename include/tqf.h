#pragma once
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <flint/qfb.h>

#include "factor_addition.h"

#define TQL_TEST_MIN_BIT_LEN 3
#define TQL_TEST_MAX_BIT_LEN 3
#define TQL_TEST_NB_PER_SAMPLE 100

#define FMPZ_DEBUG(g) flint_printf(#g " = "); fmpz_print(g); flint_printf("\n"); 

/* Represent a ternary quadratic form over ZZ with known factorization
ax^2 + by^2 + cz^2 = 0 */
typedef struct
{
    fmpz_factor_t coeffs[3];
} fmpz_tqf;

typedef fmpz_tqf fmpz_tqf_t[1];

void fmpz_tqf_init(fmpz_tqf_t C);
void fmpz_tqf_set(fmpz_tqf_t C, fmpz_factor_t a, fmpz_factor_t b, fmpz_factor_t c);
void fmpz_tqf_init_set(fmpz_tqf_t C, fmpz_factor_t a, fmpz_factor_t b, fmpz_factor_t c);
void fmpz_tqf_clear(fmpz_tqf_t C);
void fmpz_tqf_print(const fmpz_tqf_t C);

void fmpz_tqf_reduce(fmpz_tqf_t A, const fmpz_tqf_t C, fmpz_t t[3]);
int _fmpz_tqf_test_sol(const fmpz* a, const fmpz* b, const fmpz* c, const fmpz* S);
int fmpz_tqf_certif(fmpz_t k, const fmpz_t a, const fmpz_t b, const fmpz_factor_t c);
int fmpz_tqf_solve_reduced(const fmpz_tqf_t R, fmpz_t z1, fmpz_t z2, fmpz_t z3);
void fmpz_tqf_parametrize(qfb_t q1, qfb_t q2, qfb_t q3, const fmpz_tqf_t R, 
    const fmpz_t x0, const fmpz_t y0, const fmpz_t z0);
void fmpz_tqf_parametrize_semi_diag(qfb_t q1, qfb_t q2, qfb_t q3, 
    const fmpz_t d1, const fmpz_t c, const fmpz_t d2, 
    const fmpz_t delta,
    const fmpz_t x0, const fmpz_t y0, const fmpz_t z0);

void test_tqf();