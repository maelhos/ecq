#pragma once
#include <flint/fmpz.h>

/* Represent an Elliptic Curve in Weirstrass Form over ZZ
y^2 = x^3 + a_2 x^2 + a_4 x + a_6 */
typedef struct
{
    fmpz_t a2, a4, a6;
} fmpz_EC;

typedef fmpz_EC fmpz_EC_t[1];

void fmpz_EC_init(fmpz_EC_t C, const fmpz_t a2, const fmpz_t a4, const fmpz_t a6);
void fmpz_EC_init_from_torsion(fmpz_EC_t C, const fmpz_t e1, const fmpz_t e2, const fmpz_t e3);
void fmpz_EC_print(const fmpz_EC_t C);
int  fmpz_EC_is_on_curve(const fmpz_EC_t C, const fmpz_t x, const fmpz_t y);