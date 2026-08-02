#pragma once
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include "factor_addition.h"

// returns the minimal model of a given curve
void minimal_model(fmpz_t a1, fmpz_t a2, fmpz_t a3, fmpz_t a4, fmpz_t a6, fmpz_t u, const fmpz_t c4, const fmpz_t c6, const fmpz_factor_t disc);

