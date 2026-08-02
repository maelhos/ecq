#include <stdint.h>
#include <stdlib.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_vec.h>

#include "EC.h"
#include "two_descent.h"
#include "minimal_model.h"
#include "tqf.h"
#include "reduce.h"
#include "second_descent.h"
#include "descent.h"

// how many ELS descendants to keep per first descent cover
/* How many ELS descendants to keep per first descent cover.

Cremona d2 3.2.3: the d3 giving ELS descendants of a given C_d1 form one complete
coset of H1 = S^(phi')(E/Q)/K in H0, so every one of them represents the SAME class
in S^2(E) and they are all equivalent as coverings. One representative per d1 is
therefore enough, and is what makes the printed list of 2-coverings match the output
of Magma's TwoDescent.

Keeping more only helps when the first descendant has no point small enough for the
sieve to find, in which case Cremona loops through the rest of the coset; raise this
with ECQ_DESCENDANTS=n if a search stalls. */
#define ECQ_MAX_DESCENDANTS_PER_COVER 1

int main(int argc, char* argv[])
{
    fmpz_t a1, a2, a3, a4, a6;
    fmpz_t a1m, a2m, a3m, a4m, a6m;
    fmpz_t b2, b4, b6, c4, c6, disc, tmp, tmp2;
    fmpz_t u, s, r, t;
    fmpz_factor_t disc_f;

    // second descent on by default, can be removed...
    const int do_second_descent = (getenv("ECQ_NO_SECOND_DESCENT") == NULL);
    const char *nd_env = getenv("ECQ_DESCENDANTS");
    const slong max_descendants = nd_env ? atol(nd_env) : ECQ_MAX_DESCENDANTS_PER_COVER;

    fmpz_init(a1);
    fmpz_init(a2);
    fmpz_init(a3);
    fmpz_init(a4);
    fmpz_init(a6);

    flint_printf("> "); fmpz_fread(stdin, a1);
    flint_printf("> "); fmpz_fread(stdin, a2);
    flint_printf("> "); fmpz_fread(stdin, a3);
    flint_printf("> "); fmpz_fread(stdin, a4);
    flint_printf("> "); fmpz_fread(stdin, a6);

    fmpz_init(a1m);
    fmpz_init(a2m);
    fmpz_init(a3m);
    fmpz_init(a4m);
    fmpz_init(a6m);

    fmpz_init(b2);
    fmpz_init(b4);
    fmpz_init(b6);
    fmpz_init(c4);
    fmpz_init(c6);
    fmpz_init(disc);
    fmpz_init(tmp);
    fmpz_init(tmp2);

    fmpz_factor_init(disc_f);

    fmpz_mul(b2, a1, a1);
    fmpz_addmul_ui(b2, a2, 4);

    fmpz_mul(b4, a1, a3);
    fmpz_addmul_ui(b4, a4, 2);

    fmpz_mul(b6, a3, a3);
    fmpz_addmul_ui(b6, a6, 4);

    fmpz_mul(c4, b2, b2);
    fmpz_submul_ui(c4, b4, 24);

    fmpz_mul(c6, b2, b2); // c6 = b2(-b2**2 + 36b4) - 216b6
    fmpz_neg(c6, c6);
    fmpz_addmul_ui(c6, b4, 36);
    fmpz_mul(c6, c6, b2);
    fmpz_submul_ui(c6, b6, 216);

    fmpz_mul(tmp, c6, c6);
    fmpz_mul(disc, c4, c4);
    fmpz_mul(disc, disc, c4);
    fmpz_sub(disc, disc, tmp);
    fmpz_divexact_ui(disc, disc, 1728);

    fmpz_factor(disc_f, disc);
    fmpz_factor_sort(disc_f);

    flint_printf("disc = %{fmpz}\n", disc);

    fmpz_init(u);
    minimal_model(a1m, a2m, a3m, a4m, a6m, u, c4, c6, disc_f);
    
    fmpz_EC_t E, Em;
    fmpz_EC_init(E, a1, a2, a3, a4, a6);
    fmpz_EC_init(Em, a1m, a2m, a3m, a4m, a6m);

    // print
    flint_printf("E : "); fmpz_EC_print(E);
    flint_printf("Em : "); fmpz_EC_print(Em);

    // get back isomorphism parameters
    fmpz_init(r);
    fmpz_init(t);
    fmpz_init(s);

    // s = (ua1' - a1) / 2
    fmpz_mul(s, u, a1m);
    fmpz_sub(s, s, a1);
    fmpz_divexact_ui(s, s, 2);

    // r = (u^2a2' - a2 + sa1 + s^2) / 3
    fmpz_mul(r, u, u);
    fmpz_mul(r, r, a2m);
    fmpz_sub(r, r, a2);
    fmpz_mul(tmp, s, a1);
    fmpz_add(r, r, tmp);
    fmpz_mul(tmp, s, s);
    fmpz_add(r, r, tmp);
    fmpz_divexact_ui(r, r, 3);

    // t = (u^3a3' - a3 - ra1) / 2
    fmpz_mul(t, u, u);
    fmpz_mul(t, t, u);
    fmpz_mul(t, t, a3m);
    fmpz_sub(t, t, a3);
    fmpz_mul(tmp, r, a1);
    fmpz_sub(t, t, tmp);
    fmpz_divexact_ui(t, t, 2);

    flint_printf("Using isomorphism : [%{fmpz}, %{fmpz}, %{fmpz}, %{fmpz}]\n", u, r, s, t);

    // test for 2-isogenies
    fmpz_poly_t P;
    fmpz_poly_factor_t Pf;
    slong i;
    ulong u_square;

    fmpz_poly_init2(P, 4);
    fmpz_poly_factor_init(Pf);

    if (fmpz_is_zero(a1m) && fmpz_is_zero(a3m))
    {
        // work on the minimal model throughout: X is an x-coordinate on Em,
        // mapped back to E at the very end via x = u^2 x' + r
        fmpz_poly_set_coeff_fmpz(P, 0, a6m);
        fmpz_poly_set_coeff_fmpz(P, 1, a4m);
        fmpz_poly_set_coeff_fmpz(P, 2, a2m);
        fmpz_poly_set_coeff_ui  (P, 3, 1);
        u_square = 1;
    }
    else
    {
        // set aux coeff of minimal model reusing same vars...
        fmpz_mul(b2, a1m, a1m);
        fmpz_addmul_ui(b2, a2m, 4);

        fmpz_mul(b4, a1m, a3m);
        fmpz_addmul_ui(b4, a4m, 2);

        fmpz_mul(b6, a3m, a3m);
        fmpz_addmul_ui(b6, a6m, 4);

        // P = x^3 + b2 x^2 + 8b4 x + 16b6
        fmpz_mul_ui(tmp, b6, 16);
        fmpz_poly_set_coeff_fmpz(P, 0, tmp);
        fmpz_mul_ui(tmp, b4, 8);
        fmpz_poly_set_coeff_fmpz(P, 1, tmp);
        fmpz_poly_set_coeff_fmpz(P, 2, b2);
        fmpz_poly_set_coeff_ui  (P, 3, 1);

        u_square = 4;
    }

    fmpz_poly_factor(Pf, P);

    if (Pf->num > 1)
    {
        fmpz_t x0, c, d, cp, dp, dpp;
        fmpz_t d1, d2;
        fmpz_factor_t df, dpf, dppf, p_all;
        two_cover_t cover;

        two_cover_init(cover);

        fmpz_init(d1);
        fmpz_init(d2);

        fmpz_init(x0);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(dp);
        fmpz_init(cp);
        fmpz_init(dpp);
        fmpz_factor_init(df);
        fmpz_factor_init(dpf);
        fmpz_factor_init(dppf);
        fmpz_factor_init(p_all);

        flint_printf("Descent via 2-isogeny (first descent + optional second descent using ternary quadratic forms)\n");

        
        // minimal curve resulting point
        fmpq_t X, Y;
        
        fmpq_init(X);
        fmpq_init(Y);

        for (i = 0; i < Pf->num; i++)
        {
            if (fmpz_poly_degree(Pf->p + i) == 1)
            {
                //Compute the isogeny coefficients, morphisms and decude homogeneous space
                // extract root in tmp
                fmpz_neg(x0, (Pf->p + i)->coeffs);
                flint_printf("\nUsing isogeny number %ld with x-torsion %{fmpz} on the minimal model\n", i + 1, x0);

                // c = 3x0 + s2
                fmpz_mul_ui(c, x0, 3);
                fmpz_add(c, c, P->coeffs + 2);

                // d = (c + s2)*x0 + s4
                fmpz_add(d, c,P->coeffs + 2);
                fmpz_mul(d, d, x0);
                fmpz_add(d, d, P->coeffs + 1);

                // c' = -2c
                fmpz_mul_si(cp, c, -2);

                // d' = c^2 - 4d
                fmpz_mul(dp, c, c);
                fmpz_submul_ui(dp, d, 4);

                //d'' = c'^2 - 4d' = 16d
                fmpz_mul_ui(dpp, d, 16);

                if (fmpz_is_zero(d) || fmpz_is_zero(dp))
                {
                    flint_printf("Curve was singular, aborting\n");
                    flint_abort();
                }
                
                // factor dd' and start finding quartics
                fmpz_factor(df, d);
                fmpz_factor(dpf, dp);
                
                // preallocate enough in factorization to not realloc later
            
                fmpz_factor_sort(df);
                fmpz_factor_sort(dpf);  
                _fmpz_factor_mul(p_all, df, dpf);


                // manually add 2 if not present
                if ((p_all->num > 0) && (fmpz_cmp_ui(p_all->p, 2) != 0))
                {
                    _fmpz_factor_append_ui(p_all, 2, 1);
                }

                // set dppf = 2**4 * dpf
                fmpz_factor_set(dppf, df);

                if (dppf->num > 0)
                {
                    if ((fmpz_cmp_ui(dppf->p, 2) != 0))
                    {
                        _fmpz_factor_append_ui(dppf, 2, 4);
                        fmpz_factor_sort(dppf);
                    }
                    else
                    {
                        dppf->exp[0] += 4;
                    }
                }
                else
                {
                    _fmpz_factor_append_ui(dppf, 2, 4);
                    fmpz_factor_sort(dppf);
                }

                {
                    descent_side side;
                    ulong nb_norm, nb_isog;

                    // normal homogeneous spaces over E
                    side.c = c; side.d = d; side.delta = dp; side.x0 = x0;
                    side.df = df; side.deltaf = dpf; side.p_all = p_all;
                    side.mtype_first  = MORPHISM_2ISOGENY;
                    side.mtype_second = MORPHISM_SECOND_DESCENT;
                    side.label = "E ";
                    side.add_trivial = 0;   // d1 = 1 only gives the 2-torsion point

                    nb_norm = descent_side_run(cover, &side, do_second_descent,
                                               max_descendants);

                    // isogenous homogeneous spaces over E'
                    side.c = cp; side.d = dp; side.delta = dpp; side.x0 = x0;
                    side.df = dpf; side.deltaf = dppf; side.p_all = p_all;
                    side.mtype_first  = MORPHISM_ISOGENOUS_2ISOGENY;
                    side.mtype_second = MORPHISM_ISOGENOUS_SECOND_DESCENT;
                    side.label = "E'";
                    side.add_trivial = 1;   // the trivial class maps through the isogeny

                    nb_isog = descent_side_run(cover, &side, do_second_descent,
                                               max_descendants);

                    if ((nb_norm == 0) || (nb_norm != (UWORD(1) << FLINT_FLOG2(nb_norm)))
                     || (nb_isog == 0) || (nb_isog != (UWORD(1) << FLINT_FLOG2(nb_isog))))
                    {
                        flint_printf("Selmer group order is not a power of 2, aborting...\n");
                        flint_abort();
                    }

                    flint_printf("\n  Selmer: #S^(phi)(E'/Q) = %wu, #S^(phi')(E/Q) = %wu\n",
                                 nb_norm, nb_isog);
                }
            }
        }
        
        if (cover->size == 0)
        {
            flint_printf("2-cover is empty... aborting\n");
            flint_abort();
        }

        two_cover_dedup(cover);
        two_cover_print(cover);

        flint_printf("Searching for points on the %wu coverings ...\n", cover->size);
        two_descent(X, cover);

        // scale by U : X was in the scaled coordinate xi = u_square * x', so this
        // recovers x', an x-coordinate on the minimal model Em
        fmpz_mul_ui(&X->den, &X->den, u_square);
        fmpq_canonicalise(X);

        flint_printf("Output coord on Em %{fmpq}\n", X);

        // map back Em -> E along the isomorphism [u, r, s, t] : x = u^2 x' + r
        fmpz_mul(tmp, u, u);
        fmpz_mul(&X->num, &X->num, tmp);
        fmpq_canonicalise(X);
        fmpz_addmul(&X->num, &X->den, r);
        fmpq_canonicalise(X);

        flint_printf("Output coord on E %{fmpq}\n", X);

        if (!fmpz_EC_lift_x(Y, E, X))
        {
            flint_printf("Lift failed, not point with such x-coord\n");
            flint_abort();
        }   

        fmpq_canonicalise(Y);

        flint_printf("(%{fmpq} : %{fmpq} : 1)\n", X, Y);

        two_cover_clear(cover);
        fmpq_clear(X);
        fmpq_clear(Y);
        fmpz_clear(x0); fmpz_clear(c); fmpz_clear(d);
        fmpz_clear(cp); fmpz_clear(dp); fmpz_clear(dpp);
        fmpz_clear(d1); fmpz_clear(d2);
        fmpz_factor_clear(df); fmpz_factor_clear(dpf);
        fmpz_factor_clear(dppf); fmpz_factor_clear(p_all);

        fmpz_EC_clear(E);
        fmpz_EC_clear(Em);
        fmpz_poly_clear(P);
        fmpz_poly_factor_clear(Pf);
        fmpz_factor_clear(disc_f);

        // bruh, don't ask
        fmpz_clear(a1); fmpz_clear(a2); fmpz_clear(a3); fmpz_clear(a4); fmpz_clear(a6);
        fmpz_clear(a1m); fmpz_clear(a2m); fmpz_clear(a3m); fmpz_clear(a4m); fmpz_clear(a6m);
        fmpz_clear(b2); fmpz_clear(b4); fmpz_clear(b6);
        fmpz_clear(c4); fmpz_clear(c6); fmpz_clear(disc);
        fmpz_clear(tmp); fmpz_clear(tmp2);
        fmpz_clear(u); fmpz_clear(r); fmpz_clear(s); fmpz_clear(t);

        return 0;
    }
    else
    {
        flint_printf("TODO: Implement syzygy sieve for general case\n");
        return 1;
    }
}