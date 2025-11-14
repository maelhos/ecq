#include <stdint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_vec.h>

#include "EC.h"
#include "two_descent.h"
#include "tqf.h"

#define TQL_TEST_MIN_BIT_LEN 10
#define TQL_TEST_MAX_BIT_LEN 100
#define TQL_TEST_NB_PER_SAMPLE 70


int test_tqf()
{
    // non working example : 2 * 17791 * 2041792994005229x^2 + 47591875233233478767y^2 + -1 * 3 * 19191403664720039567z^2 = 0
    
    fmpz c[3];
    fmpz sol[3];
    fmpz_factor_t af, bf, cf;

    fmpz_tqf_t TBF, TBF_r;
    fmpz_t tr[3];

    int err;
    ulong l, nb;

    flint_rand_t rstate;

    flint_rand_init(rstate);

    fmpz_init(c + 0);
    fmpz_init(c + 1);
    fmpz_init(c + 2);

    fmpz_init(sol + 0);
    fmpz_init(sol + 1);
    fmpz_init(sol + 2);

    fmpz_init(tr[0]); 
    fmpz_init(tr[1]); 
    fmpz_init(tr[2]);

    // init qbf
    fmpz_tqf_init(TBF);
    fmpz_tqf_init(TBF_r);

    // init factors
    fmpz_factor_init(af);
    fmpz_factor_init(bf);
    fmpz_factor_init(cf);

    flint_printf("Testing from bitlength %ld to %ld with %ld run per sample\n", 
        TQL_TEST_MIN_BIT_LEN, TQL_TEST_MAX_BIT_LEN, TQL_TEST_NB_PER_SAMPLE);

    for (l = TQL_TEST_MIN_BIT_LEN; l <= TQL_TEST_MAX_BIT_LEN; l++)
    {
        flint_printf("Test with bitlength %ld\n", l);

        nb = 0;
        while (nb < TQL_TEST_NB_PER_SAMPLE)
        {
            // parameters for c_0 x^2 + c_1 y^2 + c_2 z^2 = 0
            do 
            {
                fmpz_randbits(c + 0, rstate, l);
                fmpz_randbits(c + 1, rstate, l);
                fmpz_randbits(c + 2, rstate, l);
            } while (fmpz_is_zero(c + 0) || fmpz_is_zero(c + 1) || fmpz_is_zero(c + 2) ||
                     ((fmpz_sgn(c + 0) == fmpz_sgn(c + 1)) && (fmpz_sgn(c + 1) == fmpz_sgn(c + 2))));

            // factor coeffs
            fmpz_factor(af, c + 0);
            fmpz_factor(bf, c + 1);
            fmpz_factor(cf, c + 2);

            // set tqf
            fmpz_tqf_set(TBF, af, bf, cf);

            // reduce qbf
            fmpz_tqf_reduce(TBF_r, TBF, tr);
            
            // solving qbf
            err = fmpz_tqf_solve_reduced(TBF_r, sol + 0, sol + 1, sol + 2);
            
            // translating
            fmpz_mul(sol + 0, sol + 0, tr[0]);
            fmpz_mul(sol + 1, sol + 1, tr[1]);
            fmpz_mul(sol + 2, sol + 2, tr[2]);

            // print solution
            if (err == 1)
            {   
                flint_printf("\nTBF : "); fmpz_tqf_print(TBF); flint_printf("\n");
                flint_printf("TBF reduced : "); fmpz_tqf_print(TBF_r); flint_printf("\n");
                flint_printf("tr : "); 
                fmpz_print(tr[0]); flint_printf(" "); 
                fmpz_print(tr[1]); flint_printf(" ");
                fmpz_print(tr[2]); flint_printf("\n");

                flint_printf("Found solution : ("); 
                fmpz_print(sol + 0); flint_printf(", ");
                fmpz_print(sol + 1); flint_printf(", ");
                fmpz_print(sol + 2); flint_printf(")\n");

                if (!_fmpz_tqf_test_sol(c + 0, c + 1, c + 2, sol))
                {
                    flint_printf("Solution FAIL !!\n");
                    flint_abort();
                }
                nb++;
            }

        }
    }

    return 0;
}


int main(int argc, char* argv[])
{
    return test_tqf();
}