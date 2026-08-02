from sage.all import *


P = ZZ["x0, y0, z0, a, b, c, d, e"]
x0, y0, z0, a, b, c, d, e = P.gens()


Px = P["u, v"]
u, v = Px.gens()

Q1 = x0*u**2 + 2*(b*x0 + 2*c*z0)*u*v + x0*d*v**2
Q2 = y0*u**2 - y0*d*v**2
Q3 = z0*u**2 - 2*(b*z0 + 2*a*x0)*u*v + z0*d*v**2

print(Q1(u=u + e*v / (2*y0), v=v / (2*y0)))
print(Q2(u=u + e*v / (2*y0), v=v / (2*y0)))
print(Q3(u=u + e*v / (2*y0), v=v / (2*y0)))