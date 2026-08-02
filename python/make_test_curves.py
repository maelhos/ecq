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
l = 14
f = y**2 -  (x - 2**l - 9)*(x**2 + 39 * 2**l * x + 2**l + 37)

print(factor(y**2 - f))

def free_QQ_EC_gen(C):
    magma_shot = f"P<x>:=PolynomialRing(RationalField());E:=EllipticCurve([{str(C.ainvs())[1:-1]}]);SetClassGroupBounds(\"GRH\");print Generators(E);"
    ret = magma_free(magma_shot)
    print(ret)
    if not ":" in ret: return []
    
    gen_str = ret.replace(" ", "").split("\n")[0].strip()[1:-1]
    gens = [C([QQ(e) for e in k[1:-1].split(":")]) for k in gen_str.split(",")]
    return gens

n = 9458877894

E = EllipticCurve(QQ, [-n**2, 0])
#E = EllipticCurve(f)

#print(factor((x)*((x + e)**2 + 6*(x + e) + 7)))
#print(factor((x - e)*(x**2 + 6*x + 1)))
G = free_QQ_EC_gen(E)[-1]
print(G)

print(E.ainvs())
print(E.gens())


"""
R<x> := PolynomialRing(Rationals());
E := EllipticCurve([0, 67776, 0, 524264742, 7716543824]);

TwoDescent(E)

[
    Hyperelliptic Curve defined by y^2 = 8092*x^4 - 119952*x^3 + 245844*x^2 -
        31152*x - 9231 over Rational Field,
    Hyperelliptic Curve defined by y^2 = 12465*x^4 + 85344*x^3 + 182646*x^2 -
        40232*x + 31237 over Rational Field,
    Hyperelliptic Curve defined by y^2 = x^4 - 82224*x^2 + 3464323752 over
    Rational Field
]
[
    Mapping from: Hyperelliptic Curve defined by y^2 = 8092*x^4 - 119952*x^3 +
        245844*x^2 - 31152*x - 9231 over Rational Field to CrvEll: E
    with equations :
    567718536*$.1^4*$.2 - 2331415632*$.1^3*$.2*$.3 + 4644210192*$.1^2*$.2*$.3^2
        - 1191849480*$.1*$.2*$.3^3 - 22592*$.2^3 + 438883938*$.2*$.3^4
    -12308289634032*$.1^6 + 50994623445360*$.1^5*$.3 -
        10051402301640*$.1^4*$.3^2 - 43960293160560*$.1^3*$.3^3 +
        92707879787460*$.1^2*$.3^4 - 75383843568012*$.1*$.3^5 +
        6168518746614*$.3^6
    $.2^3,
    Mapping from: Hyperelliptic Curve defined by y^2 = 12465*x^4 + 85344*x^3 +
        182646*x^2 - 40232*x + 31237 over Rational Field to CrvEll: E
    with equations :
    75777831*$.1^4*$.2 + 1549724292*$.1^3*$.2*$.3 + 2819789214*$.1^2*$.2*$.3^2 -
        1945296420*$.1*$.2*$.3^3 - 22592*$.2^3 - 849722153*$.2*$.3^4
    -3212532522729*$.1^6 - 10654077415842*$.1^5*$.3 + 26499282723135*$.1^4*$.3^2
        + 64794293720460*$.1^3*$.3^3 + 151159281418545*$.1^2*$.3^4 +
        121640054017278*$.1*$.3^5 - 23737811482759*$.3^6
    $.2^3,
    Mapping from: Hyperelliptic Curve defined by y^2 = x^4 - 82224*x^2 +
        3464323752 over Rational Field to CrvEll: E
    with equations :
    13704*$.1^4*$.2 - 2900924904*$.1^2*$.2*$.3^2 - 22592*$.2^3 +
        47475092697408*$.2*$.3^4
    887063604*$.1^5*$.3 - 3073075512871922208*$.1*$.3^5
    $.2^3
]
"""

P = E.gens()
print(P)
exit()
