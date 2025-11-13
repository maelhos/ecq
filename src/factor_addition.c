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