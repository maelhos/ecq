#include "factor_addition.h"

void fmpz_factor_set(fmpz_factor_t z, const fmpz_factor_t x)
{
    if (z != x)
    {
        slong i, len;
        len = x->num;
        z->sign = x->sign;
        _fmpz_factor_fit_length(z, len);
        _fmpz_factor_set_length(z, len);
        _fmpz_vec_set(z->p, x->p, len);
        for (i = 0; i < len; i++)
        {
            z->exp[i] = x->exp[i];
        }
    }
}

void _fmpz_factor_mul_square_free(fmpz_factor_t factor, const fmpz_factor_t a, const fmpz_factor_t b)
{
    slong i, j;
    int cmp;

    _fmpz_factor_fit_length(factor, a->num + b->num);
    factor->sign = a->sign * b->sign;
    
    i = 0;
    j = 0;
    factor->num = 0;

    while ((i < a->num) && (j < b->num))
    {
        cmp = fmpz_cmp(a->p + i, b->p + j);
        if (cmp < 0) // a[i] < b[i]
        {
            if (a->exp[i] & 1)
            {
                fmpz_set(factor->p + factor->num, a->p + i);
                factor->exp[factor->num] = a->exp[i] & 1;
                factor->num++; 
            }
            i++;
        }
        else if (cmp > 0) // b[i] < a[i]
        {
            if (b->exp[j] & 1)
            {
                fmpz_set(factor->p + factor->num, b->p + j);
                factor->exp[factor->num] = b->exp[j] & 1;
                factor->num++; 
            }
            j++;
        }
        else // b[i] = a[i]
        {
            if ((a->exp[i] + b->exp[j]) & 1)
            {
                fmpz_set(factor->p + factor->num, a->p + i);
                factor->exp[factor->num] = (a->exp[i] + b->exp[j]) & 1;
                factor->num++; 
            }
            i++;
            j++;
        }
    }

    while (i < a->num)
    {
        if (a->exp[i] & 1)
        {
            fmpz_set(factor->p + factor->num, a->p + i);
            factor->exp[factor->num] = a->exp[i] & 1;
            factor->num++; 
        }
        i++;
    }

    while (j < b->num)
    {
        if (b->exp[j] & 1)
        {
            fmpz_set(factor->p + factor->num, b->p + j);
            factor->exp[factor->num] = b->exp[j] & 1;
            factor->num++; 
        }
        j++;
    }
}

int fmpz_factor_over_base(fmpz_factor_t z, const fmpz_t a, const fmpz_factor_t base)
{
    fmpz_t tmp;
    slong i, exp;
    int ret;

    fmpz_init_set(tmp, a);

    z->num = 0;
    z->sign = fmpz_sgn(a);


    for (i = 0; i < base->num; i++)
    {
        exp = 0;

        while (fmpz_divisible(tmp, base->p + i))
        {
            exp++;
            fmpz_divexact(tmp, tmp, base->p + i);
        }

        if (exp > 0)
        {
            _fmpz_factor_append(z, base->p + i, exp);
        }

        // the sign is carried separately in z->sign, so a fully factored input
        // leaves tmp = +-1 : testing against +1 only would reject every negative a
        if (fmpz_is_pm1(tmp))
        {
            fmpz_clear(tmp);
            return 1;
        }
    }

    ret = (fmpz_is_pm1(tmp));
    fmpz_clear(tmp);
    return ret;
}

ulong fmpz_p_adic_order(const fmpz_t a, const fmpz_t p)
{
    fmpz_t b;
    ulong ret = 0;

    // ord_p(0) is infinite : without this guard the loop below never terminates,
    // since p divides 0 and 0/p == 0. Curves with c6 = 0 (j = 1728) hit this.
    if (fmpz_is_zero(a))
    {
        return FMPZ_P_ADIC_ORDER_INF;
    }

    fmpz_init_set(b, a);

    while (fmpz_divisible(b, p))
    {
        fmpz_divexact(b, b, p);
        ret++;
    }

    fmpz_clear(b);
    return ret;
}


void _fmpz_factor_gcd_fmpz(fmpz_factor_t factor, const fmpz_factor_t a, const fmpz_t b)
{
    slong i;
    ulong o;

    _fmpz_factor_fit_length(factor, a->num);
    factor->sign = 1;
    factor->num = 0;

    for (i = 0; i < a->num; i++)
    {
        o = fmpz_p_adic_order(b, a->p + i);

        // gcd(a, 0) = a : ord_p(0) is infinite
        if (o == FMPZ_P_ADIC_ORDER_INF)
        {
            o = a->exp[i];
        }
        else if (o > (ulong) a->exp[i])
        {
            o = a->exp[i];
        }

        if (o > 0)
        {
            fmpz_set(factor->p + factor->num, a->p + i);
            factor->exp[factor->num] = o;
            factor->num++; 
        }
    }
    
}