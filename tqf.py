from sage.all import *
from random import randrange
from itertools import product

def solCertif(b1, b2, b3):
    fcts = factor(abs(b3))

    rts = []
    ffs = []
    for f, l in fcts:
        assert l == 1
        ff = GF(f)
        ffs.append(f)
        try:
            a = -ff(b1) * ff(b2)
            if a.is_square(): rts.append(int(a.sqrt()))
            else: return None
        except: return None
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
        if (cb.norm() != 0): yield (i, j, k), cb

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
    
    alpha = (bp*c*k1) % a
    beta  = (u*ap*b*k3) % (b*c)
    gamma = (v*ap*c*k2) % (b*c)
    
    v1 = vector(ZZ, [b*c, 0, 0])
    v2 = vector(ZZ, [a*beta, a, 0])
    v3 = vector(ZZ, [alpha*beta + gamma, alpha, 1])

    vs = [v1, v2, v3]
    f = lambda x, y, z: a*x**2 + b*y**2 + c*z**2
    e = lambda x, y, z: int((QQ(f(x, y, z)) / QQ(a*b*c)) % 2)

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
    
    w = [2*vi if i == j else (vj - vi if e(*vj) == 1 else vj) for j, vj in enumerate(vs)]
    A, B, C = RR(abs(a)), RR(abs(b)), RR(abs(c))
    prec = 1e10 # kinda random ngl
    # Create the transformation matrix T
    T = diagonal_matrix(ZZ, [int(RR(prec) * sqrt(A)), int(RR(prec) * sqrt(B)), int(RR(prec) * sqrt(C))])

    red = (matrix(ZZ, w) * T).LLL() / T
    for sol in smallVects(red):
        if a*sol[0]**2 + b*sol[1]**2 + c*sol[2]**2 == 0:
            return sol
        
    print(a, b, c)
    print("FAILLL")
    exit(1)

vector_stat = {}


def solveIsotropicReducedOPTI3(a, b, c):

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
    f = lambda x, y, z: (a*x**2 + b*y**2 + c*z**2) % 4
    e = lambda x, y, z: int((QQ(f(x, y, z)) / QQ(a*b*c % 4)) % 2)

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
    
    w = [2*vi if i == j else (vj - vi if e(*vj) == 1 else vj) for j, vj in enumerate(vs)]
    
    A, B, C = RR(abs(a)), RR(abs(b)), RR(abs(c))
    prec = 1 # kinda random ngl
    # Create the transformation matrix T
    T = diagonal_matrix(ZZ, [int(RR(prec) * sqrt(A)), int(RR(prec) * sqrt(B)), int(RR(prec) * sqrt(C))])

    red = (matrix(ZZ, w) * T).LLL() / T
    for cmb, sol in smallVects(red):
        if a*sol[0]**2 + b*sol[1]**2 + c*sol[2]**2 == 0:
            vector_stat.setdefault(cmb, 0);
            vector_stat[cmb] += 1
            return sol
        
    print(a, b, c)
    print("FAILLL")
    exit(1)

nb_test = 10

if __name__ == "__main__":
    a, b, c = 7561, 5 * 1181, -1 * 2 * 2477
    (ap, bp, cp), tr = reduceIsotropic3(a, b, c)
    print(factor(a), ",", factor(b), ",", factor(c))
    print(factor(ap), ",", factor(bp), ",", factor(cp))
    print(tr)

    S = solveIsotropicReducedOPTI3(ap, bp, cp)
    print("solution =", S)
    exit()
    for length in range(10, 100):
        print(f"length test : {length}")

        nb = 0
        while nb < nb_test:
            a, b, c = randrange(-2**length, 2**length), randrange(-2**length, 2**length), randrange(-2**length, 2**length)
            if (a > 0 and b > 0 and c > 0) or (a < 0 and b < 0 and c < 0) or a == 0 or b == 0 or c == 0: continue

            (ap, bp, cp), tr = reduceIsotropic3(a, b, c)
    
            S = solveIsotropicReducedOPTI3(ap, bp, cp)
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

print(vector_stat)