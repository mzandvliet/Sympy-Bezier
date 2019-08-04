import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *

def bezier_log():
    # Exploring relationships between bezier curves and logarithms
    # Todo: projective
    x = symbols('x')
    p0 = Matrix([0,     0])
    p1 = Matrix([256,   0])
    p2 = Matrix([256,   256])
    points = (p0, p1, p2)
    bases = bezier_bases(2, x)
    bez2_x = make_bezier(points, bases)(x)[0] / 32
    bez2_x = bez2_x.subs(x, x/256)

    pprint(bez2_x)

    log2_x = log(x, 2)

    graph = plot(log2_x, bez2_x, (x, 1, 256))
