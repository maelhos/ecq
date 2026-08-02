#pragma once

#include "two_descent.h"
#include "tqf.h"
#include "second_descent.h"

/* One side of the descent via 2-isogeny.

E  : y^2 = x(x^2 + c x + d)         with the "normal" homogeneous spaces
E' : y^2 = x(x^2 + c'x + d')        with c' = -2c, d' = c^2 - 4d

Both sides run exactly the same algorithm, only with (c, d, d') replaced by
(c', d', d''), so they are described by this struct and share one implementation. */

typedef struct
{
    const fmpz *c;                      // c   or c'
    const fmpz *d;                      // d   or d'   ; d1 * d2 = d
    const fmpz *delta;                  // d'  or d''  ; = c^2 - 4d
    const fmpz *x0;                     // the 2-torsion x-coordinate

    const fmpz_factor_struct *df;       // factorization of d, drives the d1 loop
    const fmpz_factor_struct *deltaf;   // factorization of delta -> TQF coefficient 2
    const fmpz_factor_struct *p_all;    // primes where local solubility is tested

    morph_type mtype_first;             // morphism for a first descent cover
    morph_type mtype_second;            // morphism for a second descent cover

    const char *label;                  // "E" / "E'", for printing
    int add_trivial;                    // does d1 = 1 contribute a cover?
} descent_side;

// Runs one side, adding covers to `cover`. Returns the number of ELS d1 found,
// which must be a power of 2.
ulong descent_side_run(two_cover_t cover, const descent_side *S,
                       int do_second_descent, slong max_descendants);
