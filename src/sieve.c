#include "sieve.h"

int sieve_pol(ulong* U, ulong* V, const fmpz* pol, slong len, ulong H, const ulong* primes, slong nb_primes)
{
    fmpz** masks;
    fmpz* curr_mask, *ev;
    ulong p, nb_squares, Hp, idx, width, curr_len, Z; 
    slong wi, w, r, i, z_mod, X;
    nmod_t* mods;
    nn_ptr reduced, buff;
    int is_const_square, empty, ret;
    fmpz_t tmp, tmp2, base_candidate, candidates;

    fmpz_init(tmp);
    fmpz_init(tmp2);
    fmpz_init(base_candidate);
    fmpz_init(candidates);

    masks   = flint_malloc(sizeof(fmpz*) * nb_primes);
    mods    = flint_malloc(sizeof(nmod_t) * nb_primes);
    reduced = _nmod_vec_init(len);
    buff    = _nmod_vec_init(primes[nb_primes - 1] + 1); // for the biggest prime, at most a bit more than half are squares
    ev      = _fmpz_vec_init(len);

    width = 2 * H + 1;
    ret = 0;

    // precompute base candidate
    fmpz_setbit(base_candidate, width);
    fmpz_sub_ui(base_candidate, base_candidate, 1);

    for (i = 0; i < nb_primes; i++)
    {
        p = primes[i];
        masks[i] = _fmpz_vec_init(p);
        nmod_init(mods + i, p);

        // reduce H mod p (for sieving index)
        NMOD_RED(Hp, H, mods[i]);

        // reduce polynomial
        _fmpz_vec_get_nmod_vec(reduced, pol, len, mods[i]);

        // compute evals that are squares
        nb_squares = 0;
        for (w = 0; w < p; w++)
        {
            if (n_jacobi(_nmod_poly_evaluate_nmod(reduced, len, w, mods[i]), p) >= 0)
            {
                buff[nb_squares] = w;
                nb_squares++;
            }
        }

        // compute whether constant coeff is a square
        is_const_square = n_jacobi(reduced[len - 1], p) >= 0;

        // compute masks
        for (z_mod = 0; z_mod < p; z_mod++)
        {
            curr_mask = masks[i] + z_mod;
            if (z_mod == 0)
            {
                if (is_const_square)
                {
                    for (r = 1; r < p; r++)
                    {
                        idx = nmod_add(r, Hp, mods[i]);
                        fmpz_setbit(curr_mask, idx);
                    }
                    
                }
            }
            else
            {
                for (wi = 0; wi < nb_squares; wi++)
                {
                    w = buff[wi];
                    idx = nmod_add(nmod_mul(w, z_mod, mods[i]), Hp, mods[i]);
                    fmpz_setbit(curr_mask, idx);
                }
            }

            // shift it
            curr_len = p;
            while (curr_len < width)
            {
                fmpz_mul_2exp(tmp, curr_mask, curr_len);
                fmpz_or(curr_mask, curr_mask, tmp);

                curr_len *= 2;
            }

            // trim it
            fmpz_and(curr_mask, curr_mask, base_candidate);
        }
    }

    // start the actual sieving
    for (Z = 1; Z <= H; Z++)
    {
        fmpz_set(candidates, base_candidate);
        
        empty = 0;
        for (i = 0; i < nb_primes; i++)
        {
            NMOD_RED(idx, Z, mods[i]);
            fmpz_and(candidates, candidates, masks[i] + idx);
            if (fmpz_is_zero(candidates))
            {
                empty = 1;
                break;
            }
        }
        
        if (empty)
        {
            continue;
        }

        // sieve
        while (!fmpz_is_zero(candidates))
        {
            idx = fmpz_val2(candidates);
            fmpz_clrbit(candidates, idx);

            X = idx - H;
            if (Z != 1)
            {
                if (n_gcd(FLINT_ABS(X), FLINT_ABS(Z)) > 1)
                {
                    continue;
                }
            }

            // finally evaluate
            _fmpz_vec_set(ev, pol, len);
            fmpz_set_ui(tmp, 1);
            fmpz_set_ui(tmp2, 1);

            for (i = 0; i < len; i++)
            {
                fmpz_mul(ev + len - i - 1, ev + len - i - 1, tmp);
                fmpz_mul(ev + i, ev + i, tmp2);
                fmpz_mul_ui(tmp, tmp, Z);
                fmpz_mul_ui(tmp2, tmp2, X);
            }
            _fmpz_vec_sum(tmp, ev, len);

            if (fmpz_sgn(tmp) <= 0 || X == 0) // no torsion
            {
                continue;
            }

            if (fmpz_is_square(tmp))
            {
                *U = X;
                *V = Z;
                ret = 1;
                goto CLEAR;
            }
            
        }
    }
    
    CLEAR:;
    for (i = 0; i < nb_primes; i++)
    {
        _fmpz_vec_clear(masks[i], primes[i]);
    }

    fmpz_clear(tmp);
    fmpz_clear(tmp2);
    fmpz_clear(base_candidate);
    fmpz_clear(candidates);

    flint_free(masks);
    flint_free(mods);
    _nmod_vec_clear(reduced);
    _nmod_vec_clear(buff);
    _fmpz_vec_clear(ev, len);
    
    return ret;
}

void test_sieve()
{
    ulong U = 0, V = 0;
    fmpz pol[5];

    fmpz_init(pol + 0);
    fmpz_init(pol + 1);
    fmpz_init(pol + 2);
    fmpz_init(pol + 3);
    fmpz_init(pol + 4);

    fmpz_set_str(pol + 0, "34603733393418305512", 10);
    fmpz_set_str(pol + 1, "63884056232237304008", 10);
    fmpz_set_str(pol + 2, "20010515629693103504", 10);
    fmpz_set_str(pol + 3, "68575292223621408192", 10);
    fmpz_set_str(pol + 4, "195717489278972053008", 10);

    flint_printf("ret = %d\n", sieve_pol(&U, &V, pol, 5, 50, n_primes_arr_readonly(12 + 1) + 1, 12));
    flint_printf("U, V = %ld, %ld\n", U, V);
}