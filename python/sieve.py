import itertools
import numpy as np
import gmpy2

def product(l: list[int]):
    p = l[0]
    for x in l[1:]:
        p *= x
    return p


def crt_basis(mod):
    M = product(mod)
    basis = []
    for i, mi in enumerate(mod):
        basis.append(M // mi * pow(M // mi, -1, mi))
    return basis

primes = [8, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
print(len(primes))
extraprimes = [47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139]
print(len(primes + extraprimes), "primes")

psquares = {p: set(x* x % p for x in range(p)) for p in (primes + extraprimes)}

ps = primes[len(primes) // 4 : len(primes) * 3 // 4]
qs = [p for p in primes if p not in ps]
print(product(ps))
print(product(qs))

def poly(x):
    return 1 + x*x + x*x*x*x

crtp = crt_basis(ps)
resp = []
for p in ps:
    squares = set(x* x % p for x in range(p))
    # x^2 = (1 - P1 m) iff m = (1 - x^2)/P1
    res = [i for i in range(p) if poly(i) % p in squares]
    resp.append(res)
print(product([len(r) for r in resp]))

crtq = crt_basis(qs)
resq = []
for q in qs:
    squares = set(x* x % q for x in range(q))
    # x^2 = (1 - P1 m) iff m = (1 - x^2)/P1
    res = [i for i in range(q) if poly(i) % q in squares]
    resq.append(res)
print(product([len(r) for r in resq]))

P, Q = product(ps), product(qs)
PQ = P * Q
A, B = crt_basis([P, Q])
print(A, B)

pp = np.array([0], dtype=np.uint64)
for ci, ri in zip(crtp, resp):
    cr = [ci * rii for rii in ri]
    pp = np.concatenate([pp + cri for cri in cr])
print("len(pp)", len(pp))
pp = list(pp)
for i in range(len(pp)):
    pp[i] = gmpy2.mpz(A * int(pp[i]) % PQ)
pp.sort()
print("pp ready")

qq = np.array([0], dtype=np.uint64)
for ci, ri in zip(crtq, resq):
    cr = [ci * rii for rii in ri]
    qq = np.concatenate([qq + cri for cri in cr])
print("len(qq)", len(qq))
qq = list(qq)
for i in range(len(qq)):
    qq[i] = gmpy2.mpz(PQ - (B * int(qq[i]) % PQ))
qq.sort()
print("qq ready")

# multiples of P => 2n
# multiples of Q => 2n+1
ip, iq = 0, 0
xp, xq = pp[ip], qq[iq]
H = 2**34
while ip < len(pp) and iq < len(qq):
    if abs(xp - xq) < H:
        x = xp - xq
        fx = poly(x)
        #assert all(fx % p in psquares[p] for p in primes)
        if all(fx % p in psquares[p] for p in extraprimes):
            print(x)
    if xp < xq:
        ip += 1
        if ip == len(pp):
            break
        xp = pp[ip]
    else:
        iq += 1
        if iq == len(qq):
            break
        xq = qq[iq]
    

