#include "EC.h"

void fmpz_EC_init(fmpz_EC_t C, const fmpz_t a1, const fmpz_t a2, const fmpz_t a3, const fmpz_t a4, const fmpz_t a6)
{
    fmpz_init_set(C->a1, a1);
    fmpz_init_set(C->a2, a2);
    fmpz_init_set(C->a3, a3);
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
    fmpz_init_set_ui(C->a1, 0);
    fmpz_init_set_ui(C->a3, 0);

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

    flint_printf("y^2");

    int sa1 = fmpz_sgn(C->a1);
    if (sa1)
    {
        fmpz_abs(tmp, C->a1);
        flint_printf(" %c ", chr_sign(sa1));
        fmpz_print(tmp);
        flint_printf("xy");
    }

    int sa3 = fmpz_sgn(C->a3);
    if (sa3)
    {
        fmpz_abs(tmp, C->a3);
        flint_printf(" %c ", chr_sign(sa3));
        fmpz_print(tmp);
        flint_printf("y");
    }

    flint_printf(" = x^3");
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

int fmpz_EC_is_on_curve(const fmpz_EC_t C, const fmpq_t x, const fmpq_t y)
{
    fmpq_t py, px;
    fmpq_init(py);
    fmpq_init(px);
    int ret;

    // y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6
    fmpq_mul_fmpz(py, x, C->a1);
    fmpq_add(py, py, y);
    fmpq_add_fmpz(py, py, C->a3);
    fmpq_mul(py, py, y);

    fmpq_set(px, x);
    fmpq_add_fmpz(px, px, C->a2);
    fmpq_mul(px, px, x);
    fmpq_add_fmpz(px, px, C->a4);
    fmpq_mul(px, px, x);
    fmpq_add_fmpz(px, px, C->a6);

    ret = fmpq_equal(px, py);

    fmpq_clear(py);
    fmpq_clear(px);

    return ret;
}

int fmpz_EC_lift_x(fmpq_t y, const fmpz_EC_t C, const fmpq_t x)
{
    fmpq_t delta, b, c, tmp;
    fmpz_t tmpz;
    int ret = 1;

    fmpz_init(tmpz);
    fmpq_init(delta);
    fmpq_init(b);
    fmpq_init(c);
    fmpq_init(tmp);

    fmpq_set(c, x);
    fmpq_add_fmpz(c, c, C->a2);
    fmpq_mul(c, c, x);
    fmpq_add_fmpz(c, c, C->a4);
    fmpq_mul(c, c, x);
    fmpq_add_fmpz(c, c, C->a6);

    // y^2 + (a1 x + a3) y - c = 0
    fmpq_mul_fmpz(b, x, C->a1);
    fmpq_add_fmpz(b, b, C->a3);

    // delta = b**2 - 4*c
    fmpq_mul(delta, b, b);
    fmpq_mul_ui(tmp, c, 4);
    fmpq_add(delta, delta, tmp);

    fmpq_canonicalise(delta);

    // check numerator
    fmpz_sqrt(tmpz, &delta->num);
    fmpz_set(&tmp->num, tmpz);
    fmpz_mul(tmpz, tmpz, tmpz);
    if (fmpz_cmp(tmpz, &delta->num) != 0)
    {
        ret = 0;
        goto CLEANUP;
    }

    // check denominator
    fmpz_sqrt(tmpz, &delta->den);
    fmpz_set(&tmp->den, tmpz);
    fmpz_mul(tmpz, tmpz, tmpz);
    if (fmpz_cmp(tmpz, &delta->den) != 0)
    {
        ret = 0;
        goto CLEANUP;
    }

    fmpq_canonicalise(tmp);

    // tmp is now the sqrt(delta)
    fmpq_neg(y, b);
    fmpq_sub(y, y, tmp);
    fmpq_div_2exp(y, y, 1);


    CLEANUP:;
    fmpz_clear(tmpz);
    fmpq_clear(delta);
    fmpq_clear(b);
    fmpq_clear(c);
    fmpq_clear(tmp);


    return ret;
}

void fmpz_EC_iso(fmpq_t xp, fmpq_t yp, const fmpz_EC_t E, const fmpz_EC_t Ep, const fmpq_t x, const fmpq_t y, const fmpz_t r, const fmpz_t s, const fmpz_t t, const fmpz_t u)
{
    if (!fmpz_EC_is_on_curve(E, x, y))
    {
        flint_printf("(Isomorphism computation) Input, point is not on curve, aborting");
        flint_abort();
    }

    // x' = (x - r) / u^2 
    fmpq_sub_fmpz(xp, x, r);
    fmpq_div_fmpz(xp, xp, u);
    fmpq_div_fmpz(xp, xp, u);

    // y' = (y - t - su^2 x') / u^3
    fmpq_mul_fmpz(yp, xp, u);
    fmpq_mul_fmpz(yp, yp, u);
    fmpq_mul_fmpz(yp, yp, s);
    fmpq_neg(yp, yp);

    fmpq_sub_fmpz(yp, yp, t);
    fmpq_add(yp, yp, y);
    fmpq_div_fmpz(yp, yp, u);
    fmpq_div_fmpz(yp, yp, u);
    fmpq_div_fmpz(yp, yp, u);


    if (!fmpz_EC_is_on_curve(Ep, xp, yp))
    {
        flint_printf("(Isomorphism computation) Ouput, point is not on curve, aborting");
        flint_abort();
    }
}

void fmpz_EC_iso_inv(fmpq_t x, fmpq_t y, const fmpz_EC_t E, const fmpz_EC_t Ep, const fmpq_t xp, const fmpq_t yp, const fmpz_t r, const fmpz_t s, const fmpz_t t, const fmpz_t u)
{
    if (!fmpz_EC_is_on_curve(E, xp, yp))
    {
        flint_printf("(Isomorphism computation) Input, point is not on curve, aborting");
        flint_abort();
    }
    fmpq_t tmp;

    fmpq_init(tmp);

    // x = u^2x' + r
    fmpq_mul_fmpz(x, xp, u);
    fmpq_mul_fmpz(x, x, u);
    fmpq_add_fmpz(x, x, r);

    // y = u^3y' + su^2x' + t
    fmpq_mul_fmpz(y, yp, u);
    fmpq_mul_fmpz(tmp, xp, s);
    fmpq_add(y, y, tmp);
    fmpq_mul_fmpz(y, y, u);
    fmpq_mul_fmpz(y, y, u);
    fmpq_add_fmpz(y, y, t);

    fmpq_clear(tmp);
    
    if (!fmpz_EC_is_on_curve(Ep, x, y))
    {
        flint_printf("(Isomorphism computation) Ouput, point is not on curve, aborting");
        flint_abort();
    }
}
void fmpz_EC_clear(fmpz_EC_t C)
{
    fmpz_clear(C->a1);
    fmpz_clear(C->a2);
    fmpz_clear(C->a3);
    fmpz_clear(C->a4);
    fmpz_clear(C->a6);
}
