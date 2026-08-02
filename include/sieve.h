#pragma once
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz_vec.h>

/* Sieves a given polynomial for a square value using quadratic
square, TODO: use LLL based sieving */
int sieve_pol(ulong* U, ulong* V, const fmpz* pol, slong len, ulong H, const ulong* primes, slong nb_primes);
