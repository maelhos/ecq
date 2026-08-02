#pragma once

#include "two_descent.h"
#include "tqf.h"
#include "reduce.h"

/* Second descent for one first descent cover

The first descent cover C_d1 : v^2 = d1 u^4 + c u^2 + d2 has its conic
Y^2 = d1 X^2 + c XZ + d2 Z^2 parametrized by (X:Y:Z) = (q1 : q2 : q3). A point of
C_d1 needs X/Z = q1/q3 to be a square. Requiring only that the PRODUCT q1*q3 be a
square gives back a curve isomorphic to C_d1 (same I, J).
The genuine 2-covering comes from splitting the condition: for a squarefree d3,

    q1(a, b) = d3 s^2      (9)
    q3(a, b) = d3 t^2      (10)

Solving (9) as a conic and parametrizing (a, b, s) = (P1, P3, P2)(u, v) turns (10)
into the descendant quartic

    D : d3 * q3(P1(u,v), P3(u,v)) = (d3 t)^2

which is a 2-covering of E, so its points have about 1/4 of the height.

d3 runs over the squarefree divisors of the resultant Res(q1, q3) = d' (in this
code's normalisation, where disc(q1) = 4 d2 and disc(q3) = 4 d1), together with 2
and both signs.

Returns the number of everywhere locally soluble descendants added to `cover`.*/
slong second_descent_d3(two_cover_t cover,
                        const fmpz_t d1, const fmpz_t c, const fmpz_t d2, const fmpz_t x0,
                        const qfb_t q1, const qfb_t q2, const qfb_t q3,
                        const fmpz_factor_t resf,
                        const fmpz_factor_t p_all,
                        morph_type mtype,
                        slong max_covers);
