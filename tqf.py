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

    g, u, v = xgcd(b, c)
    assert g == 1, "form not reduced"
    gp, ap, bp = xgcd(a, b*c)
    assert gp == 1, "form not reduced"

    k3 = solCertif(a, b, c)
    k2 = solCertif(c, a, b)
    k1 = solCertif(b, c, a)

    if k1 == None or k2 == None or k3 == None:
        return None
    
    alpha = (bp*c*k1)   % abs(a)
    beta  = (u*ap*b*k3) % abs(b*c)
    gamma = (v*ap*c*k2) % abs(b*c)

    v1 = vector(ZZ, [b*c, 0, 0])
    v2 = vector(ZZ, [a*beta, a, 0])
    v3 = vector(ZZ, [alpha*beta + gamma, alpha, 1])

    vs = [v1, v2, v3]
    f = lambda x, y, z: (a*x**2 + b*y**2 + c*z**2)

    def e(x, y, z):
        jj = f(x, y, z)
        dd = a*b*c

        return (QQ(jj) / QQ(dd)) % 2 
    
    #assert e(b*c, 0, 0) == (b*c) % 2
    #assert e(0, a*c, 0) == (a*c) % 2
    #assert e(0, 0, a*b) == (a*b) % 2

    vi = None
    i = None
    for j, v in enumerate(vs):
        if e(*v) == 1: 
            vi = v
            i = j
            break
    assert vi != None
    
    w = [2*vi if i == j else (vj - vi if e(*vj) == 1 else vj) for j, vj in enumerate(vs)]
    
    A, B, C = RR(abs(a)), RR(abs(b)), RR(abs(c))
    prec = 1 # kinda random ngl
    # Create the transformation matrix T
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
    
if __name__ == "__main__":
    nb_test = 10

    """
    a, b, c = -2543, -150934, 79
    (ap, bp, cp), tr = reduceIsotropic3(a, b, c)
    print(factor(a), ",", factor(b), ",", factor(c))
    print(factor(ap), ",", factor(bp), ",", factor(cp))
    print(tr)

    S = solveIsotropicReduced3(ap, bp, cp)
    print("solution =", S)
    """

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