from sage.all import *

# to equation y^2 = d1 x^2 + c xz + d2 z^2
R = QQ["d1, d2, d, c, x0, y0, z0"]
d1, d2, d, c, x0, y0, z0 = R.gens()
Q = R.quo(d1*d2 - d)
d1, d2, d, c, x0, y0, z0 = Q.gens()

# notation from the cremona paper
A, B, C, D = d1, c, d2, 1
delta = B**2 - 4*A*C

Fuv = Q["U, V"]
U, V = Fuv.gens()


