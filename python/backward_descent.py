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
#e1, e2, e3 = 765, 654, -757 # medium
#e1, e2, e3 = -45171, 16881, -57427 # easy

f = y**2 - (x - e1)*(x - e2)*(x - e3)

E = EllipticCurve(f)
print(E)
G = E.gens()
print(G)

print(f"rank = {E.rank()}")
print(f"disc = {E.discriminant()}")
bad_reduction = [p for p, _ in E.discriminant().factor()] + [1, -1]
print(f"{bad_reduction = }")
# We consider primes of bad potentially reduction (we could compute more precisly but this doesnt really matter...)

S = set(prod(elt) for l in range(1, len(bad_reduction)) for elt in itertools.combinations(bad_reduction, l))
S = sorted(S)
print(f"{S = }")
print(f"{len(S) = }")

# Here we assume we have a generator to show that it's descent has smaller height coefficants
G = E.gens()[0]
print("G =", G)
Gy = G.xy()[1]
Gx = G.xy()[0]

Gxn, Gyn = Gx.numerator(), Gy.numerator()
d = Gx.denominator().sqrt()

Pr = QQ["z1, z2, z3"]
z1, z2, z3 = Pr.gens()

# Sort from smallest to biggest pairs to mitigate the over estimation of bad reduction
S_prod_sorted = sorted(itertools.product(S, repeat=2), key=lambda x: abs(x[0]) + abs(x[1]))
print("sorted...")
## Over Q we can do even better

r21 = e2 - e1
r13 = e1 - e3
r32 = e3 - e2
#from descent import reduceIsotropic3, solveIsotropicReduced3, iso3iter
print()

for b1, b2 in S_prod_sorted:
    R = b1 * b2
    b3 = prod(f**(k % 2) for f, k in factor(R))
    if b1*b2*b3 < 0: b3 *= -1

    a = b1*r32
    b = b2*r13
    c = b3*r21

    if (a > 0 and b > 0 and c > 0) or (a < 0 and b < 0 and c < 0): continue
    (ap, bp, cp), tr = reduceIsotropic3(a, b, c)
    S = solveIsotropicReduced3(ap, bp, cp)
    if S == None: continue

    I = [a * z1**2 + b * z2**2 + c*z3**2]
    Ip = Pr.ideal(I + [Gxn - (e1*d**2 + b1*z1**2), Gyn - int(sqrt(b1*b2*b3)) * z1*z2*z3])

    # Here we find the coefficient of the quadric (is it really a quadric ?) to show that indeed they are smaller
    if Ip.dimension() == 0 and len(Va := Ip.variety()) > 0:
        z1v, z2v, z3v = Va[0][z1], Va[0][z2], Va[0][z3]
        if not (z1v in ZZ and z2v in ZZ and z3v in ZZ): continue

        print(f"b1 = {factor(b1)}, b2 = {factor(b2)}, b3 = {factor(b3)}, solution : {Va[0]}")
        assert a * z1v**2 + b * z2v**2 + c*z3v**2 == 0
        print("Z =", d)

        # parametrization
        it = iso3iter(ap, bp, cp, S[0], S[1], S[2])
        z1p_space, z2p_space, z3p_space = it.sol(U, V)
        z1_space, z2_space, z3_space = tr[0]*z1p_space, tr[1]*z2p_space, tr[2]*z3p_space

        # this is purely backtest
        assert a * z1_space**2 + b * z2_space**2 + c*z3_space**2 == 0
        assert a * (tr[0] * S[0])**2 + b * (tr[1] * S[1])**2 + c*(tr[2] * S[2])**2 == 0

        eq1 = z2_space * z3v - z3_space * z2v
        eq2 = z3_space * z1v - z1_space * z3v

        L = eq1.gcd(eq2)
        c_U = L.coefficient(U) 
        c_V = L.coefficient(V)
        sol_U = -c_V
        sol_V = c_U

        v_z1, v_z2, v_z3 = it.sol(sol_U, sol_V)
        k = QQ(v_z1) / QQ(z1v)
        assert (z1v, z2v, z3v) == (v_z1 / k, v_z2 / k, v_z3 / k)
        print(f"{sol_U = }, {sol_V = }")
        print(f"projective scaling factor = {k}")
        print("z^2 =", factor((b1*v_z1**2 - b2*v_z2**2) / (e2 - e1)))
        
        # compute the new form
        sieve_eq = (b1*z1_space**2 - b2*z2_space**2) / (e2 - e1)
        print("The new equation :")
        print("H^2 =", sieve_eq)


        exit()