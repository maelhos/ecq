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
