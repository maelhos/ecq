#include "EC.h"

void fmpz_EC_init(fmpz_EC_t C, const fmpz_t a2, const fmpz_t a4, const fmpz_t a6)
{
    fmpz_init_set(C->a2, a2);
    fmpz_init_set(C->a4, a4);
    fmpz_init_set(C->a6, a6);
}

void fmpz_EC_init_from_torsion(fmpz_EC_t C, const fmpz_t e1, const fmpz_t e2, const fmpz_t e3)
{
    /*
    y^2 = x^3 + a_2 x^2 + a_4 x + a_6
        = (x - e1)(x - e2)(x - e3)
        = x^3 - (e1 + e2 + e3) x^2 + (e1e2 + e1e3 + e2e3) x - e1e2e3
    */
    fmpz_init_set(C->a2, e1);
    fmpz_add(C->a2, C->a2, e2);
    fmpz_add(C->a2, C->a2, e3);
    fmpz_neg(C->a2, C->a2);

    fmpz_init_set(C->a4, e1);
    fmpz_mul(C->a4, C->a4, e2);

    fmpz_init_set(C->a6, C->a4);
    fmpz_mul(C->a6, C->a6, e3);
    fmpz_neg(C->a6, C->a6);

    fmpz_t tmp;
    fmpz_init_set(tmp, e1);
    fmpz_mul(tmp, tmp, e3);
    fmpz_add(C->a4, C->a4, tmp);

    fmpz_mul(tmp, e2, e3);
    fmpz_add(C->a4, C->a4, tmp);
    fmpz_clear(tmp);
}

inline char chr_sign(int sg)
{
    if (sg < 0) return '-';
    return '+';
} 

void fmpz_EC_print(const fmpz_EC_t C) 
{
    fmpz_t tmp;
    fmpz_init(tmp);

    flint_printf("y^2 = x^3");
    int sa2 = fmpz_sgn(C->a2);
    if (sa2)
    {
        fmpz_abs(tmp, C->a2);
        flint_printf(" %c ", chr_sign(sa2));
        fmpz_print(tmp);
        flint_printf("x^2");
    }

    int sa4 = fmpz_sgn(C->a4);
    if (sa4)
    {
        fmpz_abs(tmp, C->a4);
        flint_printf(" %c ", chr_sign(sa4));
        fmpz_print(tmp);
        flint_printf("x");
    }

    int sa6 = fmpz_sgn(C->a6);
    if (sa6)
    {
        fmpz_abs(tmp, C->a6);
        flint_printf(" %c ", chr_sign(sa6));
        fmpz_print(tmp);
    }
    flint_printf("\n");
    fmpz_clear(tmp);
}

int fmpz_EC_is_on_curve(const fmpz_EC_t C, const fmpz_t x, const fmpz_t y)
{
    fmpz_t y2, px;
    fmpz_init(y2);
    fmpz_init(px);

    fmpz_mul(y2, y, y);
    
    fmpz_add(px, x, C->a2);
    fmpz_mul(px, px, x);
    fmpz_add(px, x, C->a4);
    fmpz_mul(px, px, x);
    fmpz_add(px, x, C->a6);
    return fmpz_equal(y2, px);
}

