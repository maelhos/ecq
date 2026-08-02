from sage.all import *
import tqf

PQ = QQ["x, y"]
x, y = PQ.gens()

Pu = QQ["X"]
X = Pu.gen()

#f = y**2 + x*y + y - (x**3 - 234*x + 1352) #(x - e)*(x**2 - 8*x + 3)
# = 42
#f = y**2  - (x - e)*(x**2 - 2*x + 2)
#f = y**2 -  (x - 7146) * (x - 530) * (x + 7676)
f = y**2 -  (x - 71) * (x**2 + 74*x + 58512)
f = y**2 -  (x+8888)*(x**2 + 58888*x + 868198)
E = EllipticCurve(f)
#print(factor((x)*((x + e)**2 + 6*(x + e) + 7)))
#print(factor((x - e)*(x**2 + 6*x + 1)))
print(E)
print(E.ainvs())

Pts = E.gens()
Pts = [P.xy() for P in Pts if P != 0]
print(Pts)
print()
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

    c = 3*x0 + s2
    d = (c + s2)*x0 + s4
    cp = -2*c
    dp = c**2 - 4*d

    print("x0 =", x0)
    print("c =", c)
    print("cp =", cp)
    print("d =", factor(d))
    print("dp =", factor(dp))

    if d*dp == 0:
        print("Anomalous curve")
        exit()

    Fd = factor(d)
    Fdp = factor(dp)

    H_list = []
    for mask in range(1 << len(Fd)):
        d1 = 1
        for i in range(len(Fd)):
            if mask & (1 << i):
                d1 *= Fd[i][0]
        d2 = d / d1
        H_list.append((d1, d2, d1*x**4 + c*x**2 + d2))

    Hp_list = []
    for mask in range(1 << len(Fdp)):
        d1p = 1
        for i in range(len(Fdp)):
            if mask & (1 << i):
                d1p *= Fdp[i][0]
        d2p = dp / d1p
        Hp_list.append((d1p, d2p, d1p*x**4 + cp*x**2 + d2p))

    print("\nbase")
    for d1, d2, H in H_list:
        print("\n", H)
        ## morphism gives Px = d1 * x**2 + x0
        ## morphism gives Py = d1 * x * y - Py
        for Px, Py in Pts:
            I = PQ.ideal([d1 * x**2 + x0 - Px , d1 * x * y - Py ])
            Iv = I.variety()
            if len(Iv) > 0 and Iv != [{x: 0, y: 0}]:
                if {x: 0, y: 0} in Iv: Iv.remove({x: 0, y: 0})
                print((Px, Py), "->", Iv[0])
        """
        # d1*x**4 + c*x**2 + d2 -> d1*x**2 + c*xz + d2 z**2
        (add, bdd, cdd), _  = tqf.semiToDiag(d1, c, d2)
        (app, bpp, cpp), tr = tqf.reduceIsotropic3(add, bdd, cdd)
        
        S = tqf.solveIsotropicReduced3(app, bpp, cpp)
        if S == None: 
            continue
        
        xd, yd, zd = [tr[i] * S[i] for i in range(3)]
        x00, y00, z00 = tqf.diagSolToSemiSol(xd, yd, zd, d1, c)

        assert y00**2 == d1*x00**2 + c*x00*z00 + d2*z00**2

        print(f"tqf base zol of {H} : {(x00, y00, z00)}")
        Q1, Q2, Q3 = tqf.iterSquareNormalSemiDiag(x, y, x00, y00, z00, dp, d1, c, d2)
        nQ = Q1*Q3
        print("New quartic :", nQ)"""

    print("\nisogenous")
    for d1p, d2p, Hp in Hp_list:
        print("\n", Hp)
        ## morphism gives Px = y**2 / (4*x**2) + x0
        ## morphism gives Py = y*(d1p*x**4 - d2p) / (8*x**3)
        for Px, Py in Pts:
            I = PQ.ideal([y**2 + (x0 - Px) * (4*x**2), y*(d1p*x**4 - d2p) - Py * (8*x**3)])
            Iv = I.variety()
            if len(Iv) > 0 and Iv != [{x: 0, y: 0}]:
                if {x: 0, y: 0} in Iv: Iv.remove({x: 0, y: 0})
                print((Px, Py), "->", Iv[0])