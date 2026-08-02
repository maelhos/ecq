#!/usr/bin/env sage
"""
Find elliptic curves with a rational 2-torsion point whose generator has a large
canonical height, to stress the point search.

Rather than looking for the generator (which is the expensive thing we are trying to
stress in the first place), the height is *predicted* from BSD. For analytic rank 1 and
trivial Sha,

    h(P) = Reg = L'(E, 1) * |E_tors|^2 / (Omega * prod_p c_p)

Checked against tests/curves/goal (n = 10239): predicted 87.23, actual 87.236.

When Sha is non trivial the real height is smaller by |Sha|, always a perfect square:
n = 10741 predicts 98.23 but the generator sits at 10.91, and 98.23/10.91 = 9.

The family below is y^2 = x^3 - n^2 x (the congruent number curves), which has full
2-torsion so the 2-isogeny descent applies. Conjecturally n is congruent exactly when
n = 5, 6, 7 mod 8, which is used as a cheap filter for odd analytic rank.
"""

from sage.all import EllipticCurve, ZZ, RR, is_squarefree


def predicted_height(E):
    if E.analytic_rank() != 1:
        return None
    Lp = E.lseries().deriv_at1()[0]
    Om = E.period_lattice().omega()
    return RR(Lp * E.torsion_order() ** 2 / (Om * E.tamagawa_product()))


def scan(lo, hi, threshold=90.0):
    out = []
    for n in range(lo, hi):
        if n % 8 not in (5, 6, 7) or not is_squarefree(n):
            continue
        try:
            E = EllipticCurve([0, 0, 0, -ZZ(n) ** 2, 0])
            h = predicted_height(E)
            if h is not None and h > threshold:
                out.append((float(h), n))
                print("n = %-8d predicted height %8.2f   curve [0,0,0,%d,0]" % (n, h, -n * n))
        except Exception:
            pass
    out.sort()
    return out


if __name__ == "__main__":
    scan(10000, 12000, threshold=90.0)
