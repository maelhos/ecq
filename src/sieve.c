#include "sieve.h"

int sieve_pol(slong* U, ulong* V, const fmpz* pol, slong len, ulong H, const ulong* primes, slong nb_primes)
{
    bitvec** masks;
    bitvec* curr_mask;
    bitvec_t base_candidate, candidates;
    fmpz_t tmp, tmp2;
    fmpz* ev;
    ulong p, nb_squares, Hp, idx, width, curr_len, Z;
    slong wi, w, r, i, z_mod, X;
    size_t li;
    nmod_t* mods;
    nn_ptr reduced, buff;
    int is_const_square, empty, ret;

    fmpz_init(tmp);
    fmpz_init(tmp2);

    masks   = flint_malloc(sizeof(bitvec*) * nb_primes);
    mods    = flint_malloc(sizeof(nmod_t) * nb_primes);
    reduced = _nmod_vec_init(len);
    buff    = _nmod_vec_init(primes[nb_primes - 1] + 1); // for the biggest prime, at most a bit more than half are squares
    ev      = _fmpz_vec_init(len);

    width = 2 * H;
    ret = 0;

    // precompute base candidate
    bv_init(base_candidate, width, 0xff);
    bv_init(candidates, width, 0);

    for (i = 0; i < nb_primes; i++)
    {
        p = primes[i];
        masks[i] = flint_malloc(sizeof(bitvec) * p);
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
            bv_init(curr_mask, width, 0);

            if (z_mod == 0)
            {
                if (is_const_square)
                {
                    for (r = 1; r < p; r++)
                    {
                        idx = nmod_add(r, Hp, mods[i]);
                        bv_set_bit(curr_mask, idx);
                    }

                }
            }
            else
            {
                for (wi = 0; wi < nb_squares; wi++)
                {
                    w = buff[wi];
                    idx = nmod_add(nmod_mul(w, z_mod, mods[i]), Hp, mods[i]);
                    bv_set_bit(curr_mask, idx);
                }
            }

            /* Replicate the period p pattern over the whole width by repeated doubling.
               bv_or_shift_into fuses the shift and the or, and is safe with dest == src
               because it walks the limbs downwards, reading only limbs it has not
               written yet. It also trims to nbits, so no separate masking with
               base_candidate is needed. */
            curr_len = p;
            while (curr_len < width)
            {
                bv_or_shift_into(curr_mask, curr_mask, curr_len);
                curr_len *= 2;
            }
        }
    }

    // start the actual sieving
    for (Z = 1; Z <= H; Z++)
    {
        bv_set(candidates, base_candidate);

        empty = 0;
        for (i = 0; i < nb_primes; i++)
        {
            NMOD_RED(idx, Z, mods[i]);
            // the AND reports emptiness itself, so the vector is only walked once
            if (!bv_and_into(candidates, masks[i] + idx))
            {
                empty = 1;
                break;
            }
        }

        if (empty)
        {
            continue;
        }

        // sieve : walk the surviving bits limb by limb rather than rescanning from
        // the start for every candidate
        for (li = 0; li < candidates->len; li++)
        {
            bv_limb word = candidates->rep[li];

            while (word)
            {
                idx = li * BITS_PER_LIMB + (ulong) __builtin_ctzll(word);
                word &= word - 1;   // clear the lowest set bit

                X = (slong) idx - (slong) H;
                if (Z != 1)
                {
                    if (n_gcd(FLINT_ABS(X), FLINT_ABS((slong) Z)) > 1)
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
                    fmpz_mul_si(tmp2, tmp2, X);
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
    }

    CLEAR:;
    for (i = 0; i < nb_primes; i++)
    {
        for (z_mod = 0; z_mod < (slong) primes[i]; z_mod++)
        {
            bv_clear(masks[i] + z_mod);
        }
        flint_free(masks[i]);
    }

    bv_clear(base_candidate);
    bv_clear(candidates);

    fmpz_clear(tmp);
    fmpz_clear(tmp2);

    flint_free(masks);
    flint_free(mods);
    _nmod_vec_clear(reduced);
    _nmod_vec_clear(buff);
    _fmpz_vec_clear(ev, len);

    return ret;
}

void test_sieve()
{
    slong U = 0;
    ulong V = 0;
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
