import sys
from sage.all import *
from tqdm import *
import secrets

def is_square(n):
    if n < 0: return False
    if n == 0: return True
    return isqrt(n) ** 2 == n

def solve_hyperelliptic_naive_sieve(coeffs, H, primes):
    mods = {}
    
    deg = len(coeffs) - 1
    
    for p in primes:
        mods[p] = {}
        # Precompute squares mod p
        squares = set([pow(i, 2, p) for i in range(p)])
        
        # Iterate all possible z mod p
        for z_res in range(p):
            allowed_x = []
            z_val = z_res
            
            # For this z_res, check all x_res
            for x_res in range(p):
                # Compute Homogeneous F(x, z) mod p
                # F(X, Z) = sum(a_i * X^i * Z^(4-i))
                val = 0
                for i in range(deg + 1):
                    term = (coeffs[i] * pow(x_res, i, p) * pow(z_val, 4 - i, p))
                    val = (val + term) % p
                
                if val in squares:
                    allowed_x.append(x_res)
            
            mods[p][z_res] = set(allowed_x)

    for Z in range(1, H + 1):
        width = 2 * H + 1
        candidates = [True] * width
        
        for p in primes:
            z_mod = Z % p
            allowed_x_residues = mods[p][z_mod]
            for r in range(p):
                if r not in allowed_x_residues:
                    start_x = -H
                    rem = start_x % p
                    first_idx = (r - rem + p) % p
                    
                    for i in range(first_idx, width, p):
                        candidates[i] = False
        for i in range(width):
            if candidates[i]:
                X = i - H
                if Z != 1: # optimization
                    a, b = abs(X), Z
                    while b: a, b = b, a % b
                    if a > 1: continue

                val = 0
                for k in range(deg, -1, -1): 
                    term = coeffs[k] * (X**k) * (Z**(4-k))
                    val += term
                
                if is_square(val):
                    return (X, Z)
    return None

def solve_hyperelliptic_optimized(coeffs, H, primes):

    f_poly = ZZ["x"](coeffs)
    allowed_roots = {}
    lead_is_square = {}
    
    # padd
    while len(coeffs) < 5: coeffs.append(0)
    a4 = coeffs[-1]

    for p in primes:
        Fp = GF(p)
        f_p = f_poly.change_ring(Fp)
        
        # Determine valid inputs w such that f(w) is a square
        valid_w = []
        for w in range(p):
            val = f_p(w)
            if val.is_square():
                valid_w.append(w)
        allowed_roots[p] = set(valid_w)
        
        lead_is_square[p] = Fp(a4).is_square()

    for Z in range(1, H + 1):
        possible_Z = True

        z_divisors = []
        z_non_divisors = []
        
        for p in primes:
            if Z % p == 0:
                if not lead_is_square[p]:
                    possible_Z = False
                    break
                z_divisors.append(p)
            else:
                z_non_divisors.append(p)
        
        if not possible_Z:
            continue

        width = 2 * H + 1
        candidates = [True] * width
        
        for p in z_non_divisors:
            z_mod = Z % p

            allowed_x_residues = set()
            for w in allowed_roots[p]:
                allowed_x_residues.add((w * z_mod) % p)
            for r in range(p):
                if r not in allowed_x_residues:
                    start_idx = (r + H) % p
                    for k in range(start_idx, width, p):
                        candidates[k] = False

        for p in z_divisors:
            r = 0 
            start_idx = (r + H) % p
            for k in range(start_idx, width, p):
                candidates[k] = False

        for i in range(width):
            if candidates[i]:
                X = i - H
                
                # GCD Primitiveness Check
                if Z != 1: 
                    if GCD(X, Z) != 1: continue
                val = 0
                for k, coeff in enumerate(coeffs):
                    val += coeff * (X**k) * (Z**(4-k))
                
                if val >= 0 and Integer(val).is_square():
                    return (X, Z)

    return None

def solve_hyperelliptic_bitset(coeffs, H, primes):
    f_poly = ZZ["x"](coeffs)
    while len(coeffs) < 5: coeffs.append(0)
    a4 = coeffs[4]
    
    width = 2 * H + 1
    masks = {} 
    
    for p in primes:
        masks[p] = {}
        Fp = GF(p)
        f_p = f_poly.change_ring(Fp)
        
        squares = set([pow(i, 2, p) for i in range(p)])
        valid_roots = [w for w in range(p) if f_p(w) in squares]
            
        a4_is_sq = Fp(a4).is_square()
        
        for z_mod in range(p):
            allowed_x_indices = 0
            
            if z_mod == 0:
                if a4_is_sq:
                    for r in range(1, p):
                        idx = (r + H) % p
                        allowed_x_indices |= (1 << idx)
            else:
                for w in valid_roots:
                    x_res = (w * z_mod) % p
                    idx = (x_res + H) % p
                    allowed_x_indices |= (1 << idx)
            
            full_mask = allowed_x_indices
            current_len = p
            while current_len < width:
                full_mask |= (full_mask << current_len)
                current_len *= 2
            
            masks[p][z_mod] = full_mask & ((1 << width) - 1)
            #print(masks[p][z_mod])
    for Z in range(1, H + 1):
        candidates = (1 << width) - 1
        
        for p in primes:
            candidates &= masks[p][Z % p]
            if not candidates:
                break
        
        if not candidates:
            continue

        while candidates:
            lowest_bit = candidates & -candidates
            candidates ^= lowest_bit
            
            idx = lowest_bit.bit_length() - 1
            X = idx - H

            if Z != 1:
                if GCD(X, Z) != 1:
                    continue
            #print((X, Z))
                  
            #
            val = 0
            for k, c in enumerate(coeffs):
                val += c * (X**k) * (Z**(4-k))
                
            if val >= 0 and Integer(val).is_square():
                continue
                return (X, Z)
                
    return None

def getTest(height, deg, nbit=16):
    num = secrets.randbelow(2*height) - height
    den = secrets.randbelow(height - 1) + 1
    r = QQ(num) / QQ(den)
    #print("r =", r)
    base = [secrets.randbits(nbit) for _ in range(deg + 1)]
    R = QQ["x"]
    ev = R(base)(x=r)
    m = ev.numerator() * ev.denominator()
    ret = [m*k for k in base]
    assert R(ret)(x=r).is_square()
    return ret

def solve_doubly_focused(coeffs, H):
    f_poly = ZZ["x"](coeffs)
    # 1. Optimal Moduli Selection (m1, m2 ~ sqrt(H))
    target = float(H)**0.5
    primes_iter = primes(3, stop=infinity)
    
    m1, m2 = 1, 1
    P1, P2 = [], []
    
    # Greedy balancing
    for p in primes_iter:
        if m1 * m2 > H * 4: break # Sufficiently large modulus
        if m1 < m2:
            m1 *= p; P1.append(p)
        else:
            m2 *= p; P2.append(p)
            
    # 2. Precompute Valid Roots (Normalized)
    # R[m] = { w in 0..m-1 | f(w) is square mod m }
    R_sets = {m1: [], m2: []}
    
    for m, P_list in [(m1, P1), (m2, P2)]:
        # Build allowed residues via CRT
        valid_residues = [0]
        curr_m = 1
        
        for p in P_list:
            Fp = GF(p)
            fp = f_poly.change_ring(Fp)
            sq = set(x**2 for x in Fp)
            good_p = [r for r in range(p) if fp(r) in sq]
            
            # Lift to new modulus
            new_valid = []
            for r in good_p:
                inv = pow(curr_m, -1, p)
                for val in valid_residues:
                    k = ((r - val) * inv) % p
                    new_valid.append(val + k * curr_m)
            
            valid_residues = new_valid
            curr_m *= p
            
        R_sets[m] = sorted(valid_residues)

    inv_m1_mod_m2 = pow(m1, -1, m2)
    inv_m2_mod_m1 = pow(m2, -1, m1)
    
    solutions = []
    
    # 3. Iterate Denominators
    for Z in range(1, H + 1):

        if GCD(Z, m1) != 1 or GCD(Z, m2) != 1:
            # Fallback or simple skip (for high speed search of coprime points)
            continue 

        S1 = [(r * Z) % m1 for r in R_sets[m1]]
        S2 = [(r * Z) % m2 for r in R_sets[m2]]
        
        valid_k = [ (s * inv_m1_mod_m2) % m2 for s in S2 ]
        
        
        list_a1 = sorted([ k * m1 for k in valid_k ])
        
        valid_j = [ ((-s) * inv_m2_mod_m1) % m1 for s in S1 ]
        list_a2 = sorted([ j * m2 for j in valid_j ])
        
        i, j = 0, 0
        len1, len2 = len(list_a1), len(list_a2)

        M = m1 * m2
        
        while i < len1 and j < len2:
            val1 = list_a1[i]
            val2 = list_a2[j]
            diff = val1 - val2
            
            if abs(diff) <= H:
                # Potential Candidate
                X = diff
                val = 0
                for k, c in enumerate(coeffs):
                    val += c * (X**k) * (Z**(4-k))
                if val >= 0 and Integer(val).is_square():
                    solutions.append((X, Z))

                if diff < 0: i += 1
                else: j += 1
            elif diff < -H:
                i += 1
            else: # diff > H
                j += 1
                
    return solutions

# test complexity
from time import time
"""
H = 10
for Hm in range(20):
    coeffs = getTest(H, 4)
    nb = get_optimal_primes(H)
    #primes = [p for p, _ in zip(Primes(), range(nb + 1)) if p != 2]

    t = time()
    pts = solve_doubly_focused(coeffs, H)
    print(H, ":", time() - t)

    H *= 2
"""

H = 2**17
coeffs = getTest(H, 4)
coeffs = [34603733393418305512, 63884056232237304008, 20010515629693103504, 68575292223621408192, 195717489278972053008]
coeffs = [1, 1, 1, 1, 1]
nb = 10

import time

def get_optimal_primes(H):
    return ceil(2 * log(H, 2))

#for nb in range(4, 30):
nb = get_optimal_primes(H)
print(f"Nb = {nb}")
primes = [p for p, _ in zip(Primes(), range(nb + 1)) if p != 2]

"""
start = time.time()
pts = solve_hyperelliptic_naive_sieve(coeffs, H, primes)
print(f"Rational points found basic (H={H}): {pts} in {time.time() - start}")

start = time.time()
pts = solve_hyperelliptic_optimized(coeffs, H, primes)
print(f"Rational points found opti (H={H}): {pts} in {time.time() - start}")
"""

start = time.time()
pts = solve_hyperelliptic_bitset(coeffs, H, primes)
print(f"Rational points found bitset (H={H}): {pts} in {time.time() - start}")
