#include <stdint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_vec.h>

#include "EC.h"
#include "two_descent.h"
#include "tqf.h"

int main_descent(int argc, char* argv[])
{
    ulong i, j;

    // set up full 2-torsion curve
    // y^2 = (x - e1)(x - e2)(x - e3)
    fmpz_t e1, e2, e3;
    fmpz_init_set_si(e1, 7265);
    fmpz_init_set_si(e2, 649);
    fmpz_init_set_si(e3, -7557);

    fmpz_EC_t curve;
    fmpz_EC_init_from_torsion(curve, e1, e2, e3);
    fmpz_EC_print(curve);

    // compute divisors
    fmpz_t r21, r31, r32;
    fmpz_init(r21);
    fmpz_init(r31);
    fmpz_init(r32);

    fmpz_sub(r21, e2, e1);
    fmpz_sub(r31, e3, e1);
    fmpz_sub(r32, e3, e2);

    // factor divisors
    fmpz_factor_t r21_f, r31_f, r32_f;
    fmpz_factor_init(r21_f);
    fmpz_factor_init(r31_f);
    fmpz_factor_init(r32_f);

    fmpz_factor(r21_f, r21); fmpz_factor_sort(r21_f);
    fmpz_factor(r31_f, r31); fmpz_factor_sort(r31_f);
    fmpz_factor(r32_f, r32); fmpz_factor_sort(r32_f);

    flint_printf("e2 - e1 = "); fmpz_factor_print(r21_f); flint_printf("\n");
    flint_printf("e3 - e1 = "); fmpz_factor_print(r31_f); flint_printf("\n");
    flint_printf("e3 - e2 = "); fmpz_factor_print(r32_f); flint_printf("\n");

    // compute "reduced" discriminant factorization (no square)
    fmpz_factor_t delta_f, tmp_f;
    fmpz_t tmp;
    fmpz_factor_init(delta_f);
    fmpz_factor_init(tmp_f);
    fmpz_init(tmp);

    _fmpz_factor_mul(tmp_f, r21_f, r31_f);
    _fmpz_factor_mul(delta_f, tmp_f, r32_f);
    delta_f->sign = 1;

    flint_printf("delta = "); fmpz_factor_print(delta_f); flint_printf("\n");

    // compute square-free divisors
    ulong sqf_div_num = 1ULL << (delta_f->num + 1); 
    fmpz_factor_struct* sqf_divisors_delta;
    sqf_divisors_delta = flint_malloc(sqf_div_num * sizeof(fmpz_factor_struct));

    flint_printf("S = {");
    for (i = 0; i < sqf_div_num; i++)
    {
        fmpz_factor_init(sqf_divisors_delta + i);
        _fmpz_factor_fit_length(sqf_divisors_delta + i, __builtin_popcountl(i >> 1));
        sqf_divisors_delta[i].sign = (i & 1) ? 1 : -1;
        sqf_divisors_delta[i].num = 0;
        for (j = 1; j < delta_f->num + 1; j++)
        {
            if (i & (1ULL << j))
            {
                fmpz_set(sqf_divisors_delta[i].p + sqf_divisors_delta[i].num, delta_f->p + j - 1);
                sqf_divisors_delta[i].exp[sqf_divisors_delta[i].num] = 1;
                sqf_divisors_delta[i].num++;
            }
        }
        fmpz_factor_expand(tmp, sqf_divisors_delta + i);
        fmpz_print(tmp);
        flint_printf(" ");

    }
    flint_printf("}\n");

    return 0;
}

int main_tqf(int argc, char* argv[])
{
    fmpz_t a, b, c;
    fmpz_t x, y, z;
    fmpz_factor_t af, bf, cf;

    fmpz_tqf_t TBF, TBF_r;
    fmpz_t tr[3];

    int err;

    // parameters for ax^2 + by^2 + cz^2 = 0
    fmpz_init_set_si(a, 7561);
    fmpz_init_set_si(b, 5 * 1181);
    fmpz_init_set_si(c, -1 * 2 * 2477);

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);

    // factor coeffs
    fmpz_factor_init(af);
    fmpz_factor_init(bf);
    fmpz_factor_init(cf);

    fmpz_factor(af, a);
    fmpz_factor(bf, b);
    fmpz_factor(cf, c);

    // init qbf
    fmpz_tqf_init_set(TBF, af, bf, cf);
    fmpz_tqf_init(TBF_r);

    fmpz_init(tr[0]); fmpz_init(tr[1]); fmpz_init(tr[2]);
    fmpz_tqf_print(TBF); flint_printf("\n");

    // reduce qbf
    fmpz_tqf_reduce(TBF_r, TBF, tr);

    // print reduction and translation
    fmpz_tqf_print(TBF_r); flint_printf("\n");
    flint_printf("tr = "); 
    fmpz_print(tr[0]); flint_printf(" "); 
    fmpz_print(tr[1]); flint_printf(" ");
    fmpz_print(tr[2]); flint_printf("\n");
    
    // solving qbf
    err = fmpz_tqf_solve_reduced(TBF_r, x, y, z);
    
    // print solution
    if (err == 1)
    {
        flint_printf("Found solution : (");
        fmpz_print(x); flint_printf(", ");
        fmpz_print(y); flint_printf(", ");
        fmpz_print(z); flint_printf(")\n");
    }
    else
    {
        flint_printf("No solution exist\n");
    }

    return 0;
}


int main(int argc, char* argv[])
{
    return main_tqf(argc, argv);
}