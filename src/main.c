#include <stdint.h>
#include <stdlib.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_vec.h>

#include "EC.h"
#include "two_descent.h"
#include "tqf.h"

int main(int argc, char* argv[])
{
    fmpz_t e1, e2, e3;
    fmpz_init_set_si(e1, 7265);
    fmpz_init_set_si(e2, 649);
    fmpz_init_set_si(e3, -7557);

    two_descent(e1, e2, e3);
}
