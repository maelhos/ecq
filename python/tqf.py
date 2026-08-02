from sage.all import *
from random import randrange
from itertools import product

def ssqrt(a, n):
    r = int(GF(n)(a).sqrt())
    return min(r, n - r)

def solCertif(b1, b2, b3):
    fcts = factor(abs(b3))

    rts = []
    ffs = []
    for f, l in fcts:
        assert l == 1
        ff = GF(f)
        ffs.append(f)
        a = -ff(b1) * ff(b2)
        if a.is_square(): rts.append(ssqrt(a, f))
        else: return None

    return int(CRT(rts, ffs))

def reduceIsotropic3(c1, c2, c3):
    cffs = [c1, c2, c3]
    sgns = [1 if c > 0 else -1 for c in cffs]
    cffs = [abs(c) for c in cffs]

    g = GCD(cffs)
    cffs = [cffs[i] // g for i in range(3)]
    translation = [1, 1, 1]

    modif = True
    while modif:
        modif = False

        for i in range(3):
            r = prod(f**(k // 2) for f, k in factor(cffs[i]))
            cffs[i] = prod(f**(k % 2) for f, k in factor(cffs[i]))
            translation[(i+1) % 3] *= r
            translation[(i+2) % 3] *= r

        for i in range(3):
            g = gcd(cffs[(i+1) % 3], cffs[(i+2) % 3])
            if g > 1:
                modif = True
                cffs = [cffs[i]*g for i in range(3)]
                break
    cffs = [cffs[i]*sgns[i] for i in range(3)]
    gt = GCD(translation)
    translation = [translation[i] // gt for i in range(3)]
    return cffs, translation

def smallVects(basis, limit=1):
    for i, j, k in product(range(-limit, limit + 1), repeat=3):
        cb = i*basis[0] + j*basis[1] + k*basis[2]
        if (cb.norm() != 0): yield cb

def solveIsotropicReduced3(a, b, c):
    if (a == 1 and b == -1) or (a == -1 and b == 1): return (1, 1, 0)
    if (a == 1 and c == -1) or (a == -1 and c == 1): return (1, 0, 1)
    if (b == 1 and c == -1) or (b == -1 and c == 1): return (0, 1, 1)

    if (a > 0) == (b > 0) == (c > 0): return None

    g, u, v = xgcd(b, c)
    assert g == 1, "form not reduced"
    gp, ap, bp = xgcd(a, b*c)
    assert gp == 1, "form not reduced"

    k3 = solCertif(a, b, c)
    k2 = solCertif(c, a, b)
    k1 = solCertif(b, c, a)

    if k1 == None or k2 == None or k3 == None:
        return None
    
    # replacement just in case
    #F = DiagonalQuadraticForm(QQ, [a, b, c])
    #return F.solve()

    alpha = (bp*c*k1)   % abs(a)
    beta  = (u*ap*b*k3) % abs(b*c)
    gamma = (v*ap*c*k2) % abs(b*c)

    v1 = vector(ZZ, [b*c, 0, 0])
    v2 = vector(ZZ, [a*beta, a, 0])
    v3 = vector(ZZ, [alpha*beta + gamma, alpha, 1])

    vs = [v1, v2, v3]
    f = lambda x, y, z: (a*x**2 + b*y**2 + c*z**2)

    abc4 = (a*b*c) % 4

    def e(x, y, z):
        return (((f(x % 4, y % 4, z % 4) % 4) >> (((abc4) & 1) ^ 1)) & 1)
    
    assert e(b*c, 0, 0) == (b*c) % 2
    assert e(0, a*c, 0) == (a*c) % 2
    assert e(0, 0, a*b) == (a*b) % 2

    vi = None
    i = None
    for j, v in enumerate(vs):
        if e(*v) == 1: 
            vi = v
            i = j
            break
    assert vi != None
    
    w = [2*vj if i == j else (vj - vi if e(*vj) == 1 else vj) for j, vj in enumerate(vs)]
    
    A, B, C = RR(abs(a)), RR(abs(b)), RR(abs(c))
    prec = 1
    T = diagonal_matrix(ZZ, [int(RR(prec) * sqrt(A)), int(RR(prec) * sqrt(B)), int(RR(prec) * sqrt(C))])
    red = (matrix(ZZ, w) * T).LLL() / T
    
    for sol in smallVects(red):
        if a*sol[0]**2 + b*sol[1]**2 + c*sol[2]**2 == 0:
            g = GCD(sol)
            return [s // g for s in sol]
        
    print(a, b, c)
    print("FAILLL")
    exit(1)

class iso3iter:
    def __init__(self, a, b, c, x0, y0, z0):
        assert (x0, y0, z0) != (0, 0, 0)
        assert gcd(a, b) == 1 and gcd(a, c) == 1 and gcd(c, b) == 1
        assert all(k < 2 for _, k in factor(a*b*c))
        assert a*x0**2 + b*y0**2 + c*z0**2 == 0
        self.a, self.b, self.c = a, b, c

        Fz02 = Zmod(z0**2)
        e = int(Fz02(x0) / Fz02(y0))
        
        self.Q1 = lambda U, V: a*x0*U**2 + (2 * (e*a*x0 + b*y0) // z0) * U*V + ((a*x0*e**2 + 2*e*b*y0 - b*x0) // z0**2) * V**2
        self.Q2 = lambda U, V: -a*y0*U**2 + (2 * (a*x0 - a*y0*e) // z0) * U*V + ((2*a*x0*e - a*y0*e**2 + b*y0) // z0**2) * V**2
        self.Q3 = lambda U, V: a*z0*U**2 + 2*a*e*U*V + ((a*e**2 + b) // z0) * V**2

    def sol(self, u, v, safe=True):
        x1, y1, z1 = self.Q1(u, v), self.Q2(u, v), self.Q3(u, v)  
        
        if safe: assert self.a*x1**2 + self.b*y1**2 + self.c*z1**2 == 0
        return x1, y1, z1
    
    def b1b2b3(self):
        return self.b1, self.b2, self.b3

def iterSquareDiscSemiDiag(u, v, disc, d1, c, d2):
    assert ZZ(disc).is_square()
    de = ZZ(disc).sqrt()

    a1 = gcd(d1, (de + c) // 2)
    a2 = d1 // a1

    Q1 = a1*(de - c) / (2*d1) * u**2 + (de + c) / (2*a1) * v**2
    Q2 = de * u*v
    Q3 = a1 * u**2 - a2 * v**2

    assert Q2**2 == d1*Q1**2 + c*Q1*Q3 + d2*Q3**2

    return Q1, Q2, Q3

def iterNormalSemiDiag(u, v, x0, y0, z0, delta, d1, c, d2):
    g, s, t = xgcd(x0, z0)
    assert y0 != 0 and g == 1
    e = (t*(2*d1*x0 + c*z0) - s*(2*d2*z0 + c*x0)) % (4*y0**2)

    Q1 = x0*u**2 + ((x0*c + 2*z0*d2 + x0*e)/y0)*u*v  + ((e*( 2*(2*z0*d2 + x0*c) + x0*e) + x0*delta)/(4*y0**2))*v**2
    Q2 = y0*u**2 + e*u*v + ((e**2 - delta)/(4*y0))*v**2
    Q3 = z0*u**2 + ((-2*x0*d1 - z0*c + z0*e)/y0)*u*v + ((e*(-2*(2*x0*d1 + z0*c) + z0*e) + z0*delta)/(4*y0**2))*v**2

    assert Q2**2 == d1*Q1**2 + c*Q1*Q3 + d2*Q3**2
    return Q1, Q2, Q3

def iterSemiDiag(u, v, x0, y0, z0, delta, d1, c, d2):
    if y0 == 0: return iterSquareDiscSemiDiag(u, v, delta, d1, c, d2)
    return iterNormalSemiDiag(u, v, x0, y0, z0, delta, d1, c, d2)

def checkFullSolveSemiDiag(d1, c, d2):
    delta = c**2 - 4*d1*d2
    if ZZ(delta).is_square():
        return True
    A, B, C = semiToDiag(d1, c, d2)
    (na, nb, nc), _ = reduceIsotropic3(A, B, C)
    s0 = solveIsotropicReduced3(na, nb, nc)
    return s0 != None

def fullSolveSemiDiag(u, v, d1, c, d2):
    delta = c**2 - 4*d1*d2
    if ZZ(delta).is_square():
        return iterSquareDiscSemiDiag(u, v, delta, d1, c, d2)
    A, B, C = semiToDiag(d1, c, d2)
    (na, nb, nc), (t1, t2, t3) = reduceIsotropic3(A, B, C)
    s0 = solveIsotropicReduced3(na, nb, nc)
    if s0 == None:
        return None
    x0, y0, z0 = s0
    x0, y0, z0 = x0*t1, y0*t2, z0*t3
    assert A*x0**2 + B*y0**2 + C*z0**2 == 0
    x0, y0, z0 = diagSolToSemiSol(x0, y0, z0, d1, c)

    assert y0**2 == d1*x0**2 + c*x0*z0 + d2*z0**2
    return iterSemiDiag(u, v, x0, y0, z0, delta, d1, c, d2)

def semiToDiag(d1, c, d2):
    # from equation of the form y**2 = d1 x**2 + c xz + d2 z**2
    # to equation x'**2 + (4d1d2 - c**2)z**2 - 4d1 y**2 with x' = 2a x + b z
    return 1, -4*d1, 4*d1*d2 - c**2

def diagSolToSemiSol(x, y, z, d1, c):
    S = (x - c*z, 2*y*d1 , 2*z*d1)
    g = gcd(S)
    return (S[0] // g, S[1] // g, S[2] // g)

if __name__ == "__main__":
    nb_test = 1000

    # test semi-diagonal form
    length = 10
    nb = 0
    while nb < nb_test:
        d1, c, d2 = randrange(-2**length, 2**length), randrange(-2**length, 2**length), randrange(-2**length, 2**length)
        if (d1 == 0 or c == 0 or d2 == 0): continue

        ad, bd, cd = semiToDiag(d1, c, d2)
        (ap, bp, cp), tr = reduceIsotropic3(ad, bd, cd)
    
        S = solveIsotropicReduced3(ap, bp, cp)
        if S == None: 
            continue

        xd, yd, zd = [tr[i] * S[i] for i in range(3)]
        x0, y0, z0 = diagSolToSemiSol(xd, yd, zd, d1, c)

        print(x0, y0, z0)
        if y0**2 == d1*x0**2 + c*x0*z0 + d2*z0**2:
            print(" : PASS !")
            nb += 1
        else:
            print(" : FAIL ! :", (x0, y0, z0), y0**2, d1*x0**2 + c*x0*z0 + d2*z0**2)
            exit()
            
        g, s, t = xgcd(x0, z0)
        if g != 1:
            print(" : FAIL ! : g != 1 for iter")
            exit()
        assert s*x0 + t*z0 == 1
        delta = c**2 - 4*d1*d2
        Fuv = QQ["u, v"]
        u, v = Fuv.gens()
        
        if y0 == 0:
            assert ZZ(delta).is_square()
            de = ZZ(delta).sqrt()    
            a1 = gcd(d1, (de + c) // 2)
            a2 = d1 // a1

            Q1 = a1*(de - c) / (2*d1) * u**2 + (de + c) / (2*a1) * v**2
            Q2 = de * u*v
            Q3 = a1 * u**2 - a2 * v**2
            print("Y ======================= 0")
        else:
            e = (t*(2*d1*x0 + c*z0) - s*(2*d2*z0 + c*x0)) % (4*y0**2)

            Q1 = x0 * u**2 + (x0*e + c*x0 + 2*d2*z0) / y0 * u*v + (x0*e**2 + x0*delta + 2*c*x0*e + 4*d2*z0*e) / (4*y0**2) * v**2
            Q2 = y0 * u**2                            + e * u*v                                 + (e**2 - delta) / (4*y0) * v**2
            Q3 = z0 * u**2 + (z0*e - c*z0 - 2*d1*x0) / y0 * u*v + (z0*e**2 + z0*delta - 2*c*z0*e - 4*d1*x0*e) / (4*y0**2) * v**2
            
            Q3 = z0 * u**2 + (z0*e - c*z0 - 2*d1*x0) / y0 * u*v + (e*(z0*e - 2*(c*z0 + 2*d1*x0)) + z0*delta) / (4*y0**2) * v**2


        if Q2**2 == d1*Q1**2 + c*Q1*Q3 + d2*Q3**2:
            print(" : PASS !")
            nb += 1
        else:
            print(" : FAIL ! :", Q2**2 - (d1*Q1**2 + c*Q1*Q3 + d2*Q3**2))
            exit()

        if y0 == 0: exit()
            ## iterate
    exit()

    # test general ternary quadratic
    for length in range(10, 100):
        print(f"length test : {length}")

        nb = 0
        while nb < nb_test:
            a, b, c = randrange(-2**length, 2**length), randrange(-2**length, 2**length), randrange(-2**length, 2**length)
            if (a > 0 and b > 0 and c > 0) or (a < 0 and b < 0 and c < 0) or a == 0 or b == 0 or c == 0: continue

            (ap, bp, cp), tr = reduceIsotropic3(a, b, c)
    
            S = solveIsotropicReduced3(ap, bp, cp)
            if S == None: 
                continue
            
            print(f"\na = {factor(a)}, b = {factor(b)}, c = {factor(c)}")
            print(f"ap = {factor(ap)}, bp = {factor(bp)}, cp = {factor(cp)}")
            print("tr =", tr)

            x, y, z = [tr[i] * S[i] for i in range(3)]
            if a*x**2 + b*y**2 + c*z**2 == 0:
                print(" : PASS !")
                nb += 1
            else:
                print(" : FAIL ! :", (x, y, z))
                exit()
