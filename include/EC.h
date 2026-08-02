#pragma once
#include <flint/fmpz.h>
#include <flint/fmpq.h>

/* Represent an Elliptic Curve in Weirstrass Form over ZZ
y^2 + a1 xy + a3 y = x^3 + a_2 x^2 + a_4 x + a_6 */
typedef struct
{
    fmpz_t a1, a2, a3, a4, a6;
} fmpz_EC;

typedef fmpz_EC fmpz_EC_t[1];

void fmpz_EC_init(fmpz_EC_t C, const fmpz_t a1, const fmpz_t a2, const fmpz_t a3, const fmpz_t a4, const fmpz_t a6);
void fmpz_EC_clear(fmpz_EC_t C);
int fmpz_EC_lift_x(fmpq_t y, const fmpz_EC_t C, const fmpq_t x);
void fmpz_EC_init_from_torsion(fmpz_EC_t C, const fmpz_t e1, const fmpz_t e2, const fmpz_t e3);
void fmpz_EC_print(const fmpz_EC_t C);

int  fmpz_EC_is_on_curve(const fmpz_EC_t C, const fmpq_t x, const fmpq_t y);
void fmpz_EC_iso(fmpq_t xp, fmpq_t yp, const fmpz_EC_t E, const fmpz_EC_t Ep, const fmpq_t x, const fmpq_t y, const fmpz_t r, const fmpz_t s, const fmpz_t t, const fmpz_t u);
void fmpz_EC_iso_inv(fmpq_t x, fmpq_t y, const fmpz_EC_t E, const fmpz_EC_t Ep, const fmpq_t xp, const fmpq_t yp, const fmpz_t r, const fmpz_t s, const fmpz_t t, const fmpz_t u);
