#pragma once

#include <flint/fmpz_factor.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz.h>

void fmpz_factor_set(fmpz_factor_t z, const fmpz_factor_t x);
void _fmpz_factor_mul_square_free(fmpz_factor_t factor, const fmpz_factor_t a, const fmpz_factor_t b);