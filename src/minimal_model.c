#include "minimal_model.h"

void minimal_model(fmpz_t a1, fmpz_t a2, fmpz_t a3, fmpz_t a4, fmpz_t a6, fmpz_t u, const fmpz_t c4, const fmpz_t c6, const fmpz_factor_t disc)
{
    fmpz_t tmp, three, tau, c4p, c6p, b2p, b4p, b6p, a, b;
    fmpz_factor_t g;
    slong i, m;
    ulong d, val; 

    fmpz_set_ui(u, 1);
    fmpz_init_set_ui(three, 3);
    fmpz_init(tau);
    fmpz_init(tmp);
    fmpz_init(c4p);
    fmpz_init(c6p);
    fmpz_init(b2p);
    fmpz_init(b4p);
    fmpz_init(b6p);
    fmpz_init(a);
    fmpz_init(b);

    fmpz_factor_init(g);

    fmpz_mul(tmp, c6, c6);
    _fmpz_factor_gcd_fmpz(g, disc, tmp);

    // set tau default
    val = fmpz_val2(c4);
    if (val == 0 && fmpz_fdiv_ui(c6, 4) == 3)
    {
        fmpz_add_ui(tau, c6, 1);
        fmpz_fdiv_q_2exp(tau, tau, 2);
        fmpz_neg(tau, tau);
    }
    else if (val > 3 && ((fmpz_fdiv_ui(c6, 32) == 0) || fmpz_fdiv_ui(c6, 32) == 8))
    {
        // outside the loop there is no scaling yet, so this is c6 itself (b is still 0 here)
        fmpz_fdiv_q_2exp(tau, c6, 3);
    }
    else
    {
        fmpz_zero(tau);
    }

    for (i = 0; i < g->num; i++)
    {
        d = g->exp[i] / 12;
        if (fmpz_cmp_ui(g->p + i, 2) == 0)
        {
            fmpz_fdiv_q_2exp(a, c4, 4*d);
            fmpz_fdiv_q_2exp(b, c6, 6*d);

            val = fmpz_val2(a);
            if (val == 0 && fmpz_fdiv_ui(b, 4) == 3)
            {
                fmpz_add_ui(tau, b, 1);
                fmpz_fdiv_q_2exp(tau, tau, 2);
                fmpz_neg(tau, tau);
            }
            else if (val > 3 && ((fmpz_fdiv_ui(b, 32) == 0) || fmpz_fdiv_ui(b, 32) == 8))
            {
                fmpz_fdiv_q_2exp(tau, b, 3);
            }
            else
            {
                d--;
                fmpz_zero(tau);
            }
        }
        else if (fmpz_cmp_ui(g->p + i, 3) == 0)
        {
            if (fmpz_p_adic_order(c6, three) == 6*d + 2)
            {
                d--;
            }
        }

        fmpz_pow_ui(tmp, g->p + i, d);
        fmpz_mul(u, u, tmp);
    }

    // c4' = c4 / u**4
    fmpz_pow_ui(tmp, u, 4);
    fmpz_divexact(c4p, c4, tmp);

    // c6' = c6 / u**6
    fmpz_pow_ui(tmp, u, 6);
    fmpz_divexact(c6p, c6, tmp);

    // a1' = c6' mod 2
    fmpz_fdiv_r_2exp(a1, c4p, 1);
    
    // a2' = (-a1' - c6') mod 3
    fmpz_add(a2, a1, c6p);
    fmpz_neg(a2, a2);
    m = fmpz_fdiv_ui(a2, 3);
    if (m == 2) { m = -1;}
    fmpz_set_si(a2, m);

    // a3' = tau + a1'a2' mod 2
    fmpz_mul(a3, a1, a2);
    fmpz_add(a3, a3, tau);
    fmpz_fdiv_r_2exp(a3, a3, 1);

    // b2' = a1' + 4a2'
    fmpz_mul_ui(b2p, a2, 4);
    fmpz_add(b2p, b2p, a1);

    // a4' = (b2'^2 - c4' - 24a1'a3') / 48
    fmpz_mul(a4, a1, a3);
    fmpz_mul_si(a4, a4, -24);
    fmpz_sub(a4, a4, c4p);
    fmpz_mul(tmp, b2p, b2p);
    fmpz_add(a4, a4, tmp);
    fmpz_divexact_ui(a4, a4, 48);

    // a6' = (36b2'(a1'a3' + 2a4') - b2'^3 - c6' - 216a3') / 864
    fmpz_mul(a6, a1, a3);
    fmpz_addmul_ui(a6, a4, 2);
    fmpz_mul(a6, a6, b2p);
    fmpz_mul_ui(a6, a6, 36);
    fmpz_mul(tmp, tmp, b2p); // tmp still had b2'^2
    fmpz_sub(a6, a6, tmp);
    fmpz_sub(a6, a6, c6p);
    fmpz_submul_ui(a6, a3, 216);
    fmpz_divexact_ui(a6, a6, 864);

    fmpz_clear(three);
    fmpz_clear(tau);
    fmpz_clear(tmp);
    fmpz_clear(c4p);
    fmpz_clear(c6p);
    fmpz_clear(b2p);
    fmpz_clear(b4p);
    fmpz_clear(b6p);
    fmpz_clear(a);
    fmpz_clear(b);

    fmpz_factor_clear(g);
}