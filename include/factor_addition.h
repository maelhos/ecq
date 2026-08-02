#pragma once

#include <flint/fmpz_factor.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz.h>

// returned by fmpz_p_adic_order for a == 0, where the order is infinite
#define FMPZ_P_ADIC_ORDER_INF UWORD_MAX

ulong fmpz_p_adic_order(const fmpz_t a, const fmpz_t p);
void _fmpz_factor_gcd_fmpz(fmpz_factor_t factor, const fmpz_factor_t a, const fmpz_t b);
void fmpz_factor_set(fmpz_factor_t z, const fmpz_factor_t x);
void _fmpz_factor_mul_square_free(fmpz_factor_t factor, const fmpz_factor_t a, const fmpz_factor_t b);
int fmpz_factor_over_base(fmpz_factor_t z, const fmpz_t a, const fmpz_factor_t base); // returns 1 if ok

static inline void fmpz_factor_set_one(fmpz_factor_t factor)
{
    factor->sign = 1;
    factor->num = 0;
}

// assume sorted unique input and returns sorted unique
static inline void _fmpz_factor_mul_2_exp(fmpz_factor_t factor, ulong exp)
{
    if (factor->num == 0)
    {
        _fmpz_factor_append_ui(factor, 2, exp);
    }
    else
    {
        // input is sorted, so 2 can only ever be the first prime
        if (fmpz_cmp_ui(factor->p, 2) == 0)
        {
            factor->exp[0] += exp;
        }
        else
        {
            _fmpz_factor_append_ui(factor, 2, exp);
            fmpz_factor_sort(factor);
        }
    }

}