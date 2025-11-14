#pragma once
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include "factor_addition.h"

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
int fmpz_tqf_solve_reduced(fmpz_tqf_t R, fmpz_t z1, fmpz_t z2, fmpz_t z3);