from sage.all import *
from tqf import *
from local_solubility import ELS, disc, invariants

PQ = QQ["x, y"]
x, y = PQ.gens()

Pu = QQ["X"]
X = Pu.gen()

#f = y**2 + x*y + y - (x**3 - 234*x + 1352) #(x - e)*(x**2 - 8*x + 3)
# = 42
#f = y**2  - (x - e)*(x**2 - 2*x + 2)
f = y**2 -  (x+8888)*(x**2 + 58888*x + 868198)
E = EllipticCurve(f)
#print(factor((x)*((x + e)**2 + 6*(x + e) + 7)))
#print(factor((x - e)*(x**2 + 6*x + 1)))
print(E.ainvs())

"""
R<x> := PolynomialRing(Rationals());
E := EllipticCurve([0, 67776, 0, 524264742, 7716543824]);

TwoDescent(E)

[
    Hyperelliptic Curve defined by y^2 = 5540*x^4 + 45816*x^3 + 105612*x^2 -
        203092*x + 139353 over Rational Field,
    Hyperelliptic Curve defined by y^2 = x^4 - 82224*x^2 + 3464323752 over
    Rational Field,
    Hyperelliptic Curve defined by y^2 = 46831*x^4 - 136636*x^3 + 7494*x^2 +
        104636*x + 9583 over Rational Field
]
"""

P = E.gens()
print(P)
a1, a2, a3, a4, a6 = E.ainvs()

if a1 == a3 == 0:
    U = 1
    s2 = a2
    s4 = a4
    s6 = a6
else:
    U = 2
    s2 = a1**2 + 4*a2
    s4 = 8*(a1*a3 + 2*a4)
    s6 = 16*(a3**2 + 4*a6)

rts = (X**3 + s2*X**2 + s4*X + s6).roots()
for ax0 in rts:
    x0 = ax0[0]
    if a1 == a3 == 0:
        morph = (-x0, 0, 0, U)
    else:
        morph = (-x0, 0, 0, U) #(-x0, -a1/2, (-a3 + x0*a1) / 2, u)
        

    print("x0 =", x0)
    c = 3*x0 + s2
    d = (c + s2)*x0 + s4
    cp = -2*c
    dp = c**2 - 4*d

    if d*dp == 0:
        print("Anomalous curve")
        exit()

    Fd = factor(d)
    Fdp = factor(dp)
    F16dp = Fdp * factor(2)
    F16dpp = Fd * factor(2)

    # first descent
    H_list = []
    for mask in range(1 << len(Fd)):
        d1 = 1
        for i in range(len(Fd)):
            if mask & (1 << i):
                d1 *= Fd[i][0]

        if d1 == 1: continue
        d2 = d / d1
        if ELS(d1, 0, c, 0, d2):
            H_list.append(("N", d1, c, d2))

    for mask in range(1 << len(Fdp)):
        d1p = 1
        for i in range(len(Fdp)):
            if mask & (1 << i):
                d1p *= Fdp[i][0]
        d2p = dp / d1p
        if ELS(d1p, 0, cp, 0, d2p):
            H_list.append(("I", d1p, cp, d2p))
    # second descent
    
    extract_d1_c_d2 = lambda p: (int(p.coefficient(x**2)), int(p.coefficient(x*y)), int(p.coefficient(y**2)))
    C_list = []
    for t, d1, c, d2 in H_list:
        s = fullSolveSemiDiag(x, y, d1, c, d2)
        if s == None:
            continue
        q1, q2, q3 = s
        
        Fcts = F16dp if t == "N" else F16dpp
        for mask in range(1 << len(Fcts)):
            d3 = 1
            for i in range(len(Fcts)):
                if mask & (1 << i):
                    d3 *= Fcts[i][0]
            # we want d3 * q1 = y^2 and d3 * q3 = y^2
            d11, c1, d21 = extract_d1_c_d2(d3*q1)
            s1 = fullSolveSemiDiag(x, y, d11, c1, d21)
            if s1 == None: continue
            q11, q21, q31 = s1

            d13, c3, d23 = extract_d1_c_d2(d3*q3)
            if not checkFullSolveSemiDiag(d13, c3, d23): continue

            new_quartic = q3(x=q11, y=q31)
            vls = (int(new_quartic.coefficient(x**4)), 
                   int(new_quartic.coefficient(x**3*y)), 
                   int(new_quartic.coefficient(x**2*y**2)), 
                   int(new_quartic.coefficient(x*y**3)),
                   int(new_quartic.coefficient(y**4)))
            if ELS(*vls):
                I, J = invariants(*vls)
                print(HyperellipticCurve(sum(X**i * vls[5 - i - 1] for i in range(5))).isomorphisms())
                print("=>", vls, (factor(I), factor(J)))

    exit()
    
    # find points
    for t, d1, c, d2 in H_list:
        print((d1, 0, c, 0, d2))
        H = HyperellipticCurve(d1*X**4 + c*X**2 + d2)
        l_rat = H.rational_points(bound=100)
        for P in l_rat:
            if P[0] != 0 and P[1] != 0:
                u, v, _ = P
                if t == "N":
                    print("  ==>", E.lift_x((d1 * u**2 + x0) / U**2) )
                elif t == "I":
                    print("  ==>", E.lift_x(((v**2 / (4*u**2) + x0)) / U**2))
