#pragma once
#include <flint/fmpz_factor.h>
#include <flint/fmpq.h>
#include <flint/qfb.h>

#include "factor_addition.h"
#include "bkf.h"
#include "tqf.h"

typedef struct descent_eq
{
    qfb_t q1, q2, q3;
    fmpz_t b1, b2, b3;
    bkf_t sieve_eq;
} descent_eq;

typedef descent_eq descent_eq_t[1];

void descent_eq_init(descent_eq_t eq);
void descent_eq_clear(descent_eq_t eq);

int two_descent(const fmpz_t e1, const fmpz_t e2, const fmpz_t e3);