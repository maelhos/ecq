#pragma once
#include <flint/fmpz_factor.h>
#include <flint/ulong_extras.h>
#include <flint/fmpq.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq_poly.h>
#include <flint/qfb.h>

#include "morphism.h"
#include "factor_addition.h"
#include "sieve.h"
#include "bkf.h"

typedef struct
{
    // morphism from the homogeneous space back to the curve
    morph_data_t morphism;
    morph_type mtype;

    bkf_t eq;
} descent_eq;
typedef descent_eq descent_eq_t[1];

typedef struct
{
    descent_eq* cover;
    ulong size;
    ulong alloc;
} two_cover;
typedef two_cover two_cover_t[1];

void two_descent(fmpq_t X, const two_cover_t cover);
void add_cover_2isogeny(two_cover_t cover, const fmpz_t d1, const fmpz_t c, const fmpz_t d2, const fmpz_t x0, morph_type type);
void add_cover_second_descent(two_cover_t cover, const fmpz_t d1, const fmpz_t d3, const fmpz_t x0, const bkf_t q,
    const bkf_t Q1, const bkf_t Q2, const bkf_t Q3, morph_type type);
void two_cover_print(const two_cover_t cover);
void two_cover_dedup(two_cover_t cover);
void two_cover_init(two_cover_t eq);
void two_cover_clear(two_cover_t eq);
