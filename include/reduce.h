#pragma once

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/qfb.h>

#include "bkf.h"

/* Reduction of integral binary quartic forms, following
J. E. Cremona, "Reduction of binary cubic and quartic forms" */

// C(x, y) <- C(M00 x + M01 y, M10 x + M11 y)
void bkf_apply_mat(bkf_t C, const fmpz_mat_t M);

// q(x, y) <- q(M00 x + M01 y, M10 x + M11 y)
void qfb_apply_mat(qfb_t q, const fmpz_mat_t M);

// Reduce C in place, returning the substitution used in M.
// Returns 1 on success, 0 if the form is degenerate (disc == 0) and was left alone.
int bkf_reduce(bkf_t C, fmpz_mat_t M);
