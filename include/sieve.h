#pragma once
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz_vec.h>

#include "bitvec.h"

/* Sieves a given polynomial for a square value using quadratic
square, TODO: use LLL based sieving */
/* U is signed: the sieve searches x in [-H, H), so a perfectly good solution can have
   a negative numerator. V = Z is always >= 1. */
int sieve_pol(slong* U, ulong* V, const fmpz* pol, slong len, ulong H, const ulong* primes, slong nb_primes);
