#pragma once
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/qfb.h>

#include "bkf.h"

typedef struct
{
    fmpz_t x0;
    fmpz_t d1;
    fmpz_t d3;
    bkf_t Q1;
    bkf_t Q2;
    bkf_t Q3;
} morph_data;
typedef morph_data morph_data_t[1];

typedef enum
{
    MORPHISM_2ISOGENY,
    MORPHISM_ISOGENOUS_2ISOGENY,
    MORPHISM_SECOND_DESCENT,
    MORPHISM_ISOGENOUS_SECOND_DESCENT
} morph_type;

// Morphism memory convention :

/*
MORPHISM_2ISOGENY :
psi(u, v) = d1 u^2

MORPHISM_ISOGENOUS_2ISOGENY :
psi_(u, v) = v^2 / (4 u^2)

Second descent. For a squarefree d3 dividing the
resultant of q1 and q3 we split the "q1/q3 is a square" condition as

    q1(a, b) = d3 s^2    and    q3(a, b) = d3 t^2,

solve the first conic and parametrize (a, b, s) = (P1, P3, P2)(u, v). Then
Q1 = q1(P1,P3), Q2 = q2(P1,P3), Q3 = q3(P1,P3) are binary quartics and the
descendant searched is  d3 * Q3(u, v) = (d3 t)^2. On it,

    U^2 = q1/q3 = Q1/Q3       and       V = q2/q3 = Q2/Q3

recovers a point of the first descent cover, so the maps below are the same as for
the first descent, only evaluated on the composed quartics.

MORPHISM_SECOND_DESCENT :
psi(u, v) = d1 U^2 = d1 Q1(u,1) / Q3(u,1)

MORPHISM_ISOGENOUS_SECOND_DESCENT :
psi_(u, v) = V^2 / (4 U^2) = Q2(u,1)^2 / (4 Q1(u,1) Q3(u,1))
*/

void morphx(fmpq_t x, const morph_data_t morph, const fmpq_t u, const fmpq_t v, morph_type mtype);
void morph_data_clear(morph_data_t morph, morph_type mtype);
