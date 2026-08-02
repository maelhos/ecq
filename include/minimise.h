#pragma once

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "bkf.h"

/* Minimisation of integral binary quartic forms, following
J. E. Cremona, "Higher descents on elliptic curves" section 3.2.5, as implemented
in eclib's minim.cc.

Reduction (see reduce.h) only applies SL2(Z) and so cannot change the invariants
I and J. A descendant coming out of the second descent is usually also NON-minimal:
its invariants are I = w^4 I0 and J = w^6 J0 for some w > 1, and dividing those out
shrinks the coefficients much further than reduction alone can.

For each prime p with p^4 | I and p^6 | J the form is moved so that its multiple
root mod p sits at 0, then rescaled; this is not unimodular (the substitution has
determinant a power of p and the form is divided by a power of p), which is exactly
why it can lower I and J.

The returned M is the accumulated substitution only:
    C_min(x, y)  =  C(M00 x + M01 y, M10 x + M11 y) / u^2
for the power of p tracked internally in u. Companion forms that must follow the
same change of variable (the Q1, Q2, Q3 of the second descent) should be pushed
through M with bkf_apply_mat; the u factor does not affect them since the morphism
only uses ratios of the Qi. */

// Minimise C in place, accumulating the substitution in M. Returns the number of
// times the invariants were divided down (0 if C was already minimal).
slong bkf_minimise(bkf_t C, fmpz_mat_t M);

// I = 12ae - 3bd + c^2 and J = (72ae + 9bd - 2c^2)c - 27(ad^2 + b^2 e)
void bkf_invariants(fmpz_t I, fmpz_t J, const bkf_t C);

/* Test whether C1 and C2 are GL2(Q)-equivalent, i.e. define the same covering.
Necessary condition: their invariants agree up to (lambda^4, lambda^6). After
minimising and reducing both, an actual equivalence is a unimodular matrix with
small entries, which is what is searched for (bounded by BKF_EQUIV_SEARCH). */
#define BKF_EQUIV_SEARCH 60
#define BKF_EQUIV_SLIDE 12

int bkf_equivalent(const bkf_t C1, const bkf_t C2, fmpz_mat_t M);
