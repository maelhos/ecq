from sage.all import *
from tqf import *
import itertools

PP = QQ["Z"]
Z = PP.gen()

PQ = QQ["x, y"]
x, y = PQ.gens()

UV = PolynomialRing(QQ, "U, V, W")
U, V, W = UV.gens()

e1, e2, e3 = 7265, 649, -7557 # hard
#e1, e2, e3 = -45171, 16881, -57427 # hard
#e1, e2, e3 = 765, 654, -757 # medium
#e1, e2, e3 = -745, 164, -572 # easy

f = y**2 - (x - e1)*(x - e2)*(x - e3)

E = EllipticCurve(f)
print(E)

print(f"rank = {E.rank()}")
print(f"disc = {E.discriminant()}")
bad_reduction = [p for p, _ in E.discriminant().factor()] + [1, -1]
print(f"{bad_reduction = }")
# We consider primes of bad potentially reduction (we could compute more precisly but this doesnt really matter...)

S = set(prod(elt) for l in range(1, len(bad_reduction)) for elt in itertools.combinations(bad_reduction, l))
SS = sorted(list(S), key=lambda x: abs(x))
print(f"{SS = }")
print(f"{len(SS) = }")


Pr = QQ["z1, z2, z3"]
z1, z2, z3 = Pr.gens()

r21 = e2 - e1
r13 = e1 - e3
r32 = e3 - e2
#from descent import reduceIsotropic3, solveIsotropicReduced3, iso3iter
print()

sieve_eqs = []

DEBUG_CTR = 0
for i in range(len(SS)):
    for j in range(len(SS)):
        b1, b2 = SS[i], SS[j]

        R = b1 * b2
        b3 = prod(f**(k % 2) for f, k in factor(R))
        if b1*b2*b3 < 0: b3 *= -1

        a = b1*r32
        b = b2*r13
        c = b3*r21

        if (a > 0 and b > 0 and c > 0) or (a < 0 and b < 0 and c < 0): continue

        # a * z1**2 + b * z2**2 + c*z3**2 = 0
        (ap, bp, cp), tr = reduceIsotropic3(a, b, c)
        S = solveIsotropicReduced3(ap, bp, cp)

        if S == None: continue
        if S[0] == 0 or S[1] == 0 or S[2] == 0: continue

        #if DEBUG_CTR == 4:
        #    exit()
        #else:
        #    DEBUG_CTR += 1
        
        # test base solution
        test = (b1*(tr[0] * S[0])**2 - b2*(tr[1] * S[1])**2) / r21
        #if test == 0: continue

        if test.is_square() and test != 0:
            z = test.sqrt()
            z1, z2, z3 = tr[0] * S[0], tr[1] * S[1], tr[2] * S[2]
            x = (b1*z1**2 + e1*z**2) / z**2
            y = (b1*b2*b3).sqrt() * z1*z2*z3 / z**3
            print(f"({x} : {y})")
            assert y**2 - (x - e1)*(x - e2)*(x - e3) == 0

        # parametrization

        it = iso3iter(ap, bp, cp, S[0], S[1], S[2])
        z1p_space, z2p_space, z3p_space = it.sol(U, V)
        z1_space, z2_space, z3_space = tr[0]*z1p_space, tr[1]*z2p_space, tr[2]*z3p_space

        # compute the new form
        sieve_eq = (b1*z1_space**2 - b2*z2_space**2) / r21
        sieve_eqs.append((z1_space, z2_space, z3_space, b1, b2, b3, sieve_eq))
print("no basic : start sieving")

B = 1

while True:
    for u in range(B):
        if gcd(u, B) != 1: continue
        v = B - u

        for z1_s, z2_s, z3_s, b1, b2, b3, eq in sieve_eqs:
            for o in (-1, 1):
                test = eq(U=u, V=v*o)

                if test.is_square() and test != 0:
                    z = test.sqrt()
                    z1, z2, z3 = z1_s(U=u, V=v*o), z2_s(U=u, V=v*o), z3_s(U=u, V=v*o)
                    x = (b1*z1**2 + e1*z**2) / z**2
                    y = (b1*b2*b3).sqrt() * z1*z2*z3 / z**3
                    print(f"({x} : {y}) using sieve {u, v} with {z1, z2, z3} and {test = }")
                    assert y**2 - (x - e1)*(x - e2)*(x - e3) == 0
    B += 1