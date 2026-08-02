from sage.all import *


# (a, b, c, d, e)
eval4 = lambda a, b, c, d, e, x: e + x*(d + x*(c + x*(b + x*a)))
eval4diff = lambda a, b, c, d, x: d + x*(2*c + x*(3*b + 4*a*x))

def lemma6(a, b, c, d, e, x, p, n):
    tmp = eval4(a, b, c, d, e, x)

    if GF(p)(tmp).is_square():
        return 1
    
    l = ZZ(tmp).valuation(p)
    tmp = eval4diff(a, b, c, d, x)

    if tmp == 0:
        if l >= 2*n: return 0
        return -1
    
    m = ZZ(tmp).valuation(p)

    if (l >= m + n) and (n > m): return 1
    if (l >= 2*n) and (m >= n): return 0

    return -1;

def lemma7(a, b, c, d, e, x, n):
    tmp = eval4(a, b, c, d, e, x)

    l = ZZ(tmp).valuation(2) if tmp else 0
    tmp //= 2**l

    if (l % 2 == 0):
        if tmp % 8 == 1:
            return 1

    gx_odd = tmp % 4
    if (gx_odd == 0):
        gx_odd += 1

    tmp = eval4diff(a, b, c, d, x)
    if tmp == 0:
        if l >= 2*n: return 0
        if (l == 2*n-2) and (gx_odd == 1): return 0
        return -1

    m = ZZ(tmp).valuation(2) if tmp else 0

    if (l >= m + n) and (n > m): return 1
    if (n > m) and (l == m + n - 1) and (l % 2 == 0): return 1
    if (n > m) and (l == m + n - 2) and (gx_odd == 1) and (l % 2 == 0): return 1
    if (m >= n) and (l >= 2*n): return 0
    if (m >= n) and (l == 2*n - 2) and (gx_odd == 1): return 0

    return -1

def Zp_sol(a, b, c, d, e, p, x, k):
    code = ...
    if p == 2:
        code = lemma7(a, b, c, d, e, x, k)
    else:
        code = lemma6(a, b, c, d, e, x, p, k)

    if code == 1:
        return True
    if code == -1:
        return False
    
    for t in range(p):
        if Zp_sol(a, b, c, d, e, p, t*p**k + x, k + 1):
            return True
    
    return False

def Qp_sol(a, b, c, d, e, p):
    if Zp_sol(a, b, c, d, e, p, 0, 0):
        return True
    
    if Zp_sol(e, d, c, b, a, p, 0, 1):
        return True
    
    return False

def RR_Sol(a, b, c, d, e):
    if a > 0: return True
    _X = RR["X"].gen()
    ev = eval4(a, b, c, d, e, _X)
    return len(ev.roots()) > 0

def invariants(a, b, c, d, e):
    I = 12*a*e - 3*b*d + c**2
    J = 72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*e*b**2 - 2*c**3
    return I, J

def disc(a, b, c, d, e):
    I, J = invariants(a, b, c, d, e)
    return 4*I**3 - J**2

def ELS(a, b, c, d, e):
    a, b, c, d, e = int(a), int(b), int(c), int(d), int(e)
    if not RR_Sol(a, b, c, d, e): return False

    for p , _ in ZZ(2*disc(a, b, c, d, e)).factor():
        if not Qp_sol(a, b, c, d, e, p): return False
    
    return True
    

