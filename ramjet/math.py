from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.util import *

''' Polynomial helpers '''

def to_polynomial(expr, param):
    expr = simplify(expr)
    expr = collect(expr, param)  # collect in terms of a*t^0, b*t^1, c*t^2, ...
    poly = Poly(expr, param)
    return poly

# not really needed, but it was educational to write
def solve_quadratic(expr, t):
    poly = to_polynomial(expr, t)
    print("Got polynomial of degree: " + str(poly.degree()))

    coeffs = poly.coeffs()
    v1, v2, v3 = symbols('v1 v2 v3')
    substitutions = {
        v1: coeffs[0],
        v2: coeffs[1],
        v3: coeffs[2],
    }

    expr_simple = expr.subs(invert_dict(substitutions))

    # solve v1*t^2 + v2*t + v3 = 0
    eq = Eq(expr_simple, 0)
    roots = solve(eq, t)

    root_a = roots[0].subs(substitutions)
    root_b = roots[1].subs(substitutions)

    return (root_a, root_b)

''' Some Linear Algebra helpers'''

def symbolic_vector_2d(name):
    bases = ('x', 'y')
    bases = map(lambda c: name + "_" + c + ", ", bases)
    bases = reduce(lambda a, b: a + b, bases)
    x, y = symbols(bases)
    return Matrix([x, y])


def symbolic_vector_3d(name):
    bases = ('x', 'y', 'z')
    bases = map(lambda c: name + "_" + c + ", ", bases)
    bases = reduce(lambda a, b: a + b, bases)
    x, y, z = symbols(bases)
    return Matrix([x, y, z])


def dot_2d(p1, p2):
    return p1[0] * p2[0] + p1[1] * p2[1]

def cross_2d(p1, p2):
    return p1[0] * p2[1] - p1[1] * p2[0]


''' de Casteljau Bezier '''

def bernstein_basis(n, i, param):
    basis = binomial(n, i) * param**i * (1 - param)**(n-i)

    return basis


def bezier_bases(n, param):
    bases = []
    for i in range(0, n+1):
        bases.append(bernstein_basis(n, i, param))

    return bases


def make_bezier_expr(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i" % (
            len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    expr = reduce((lambda x, y: x + y), terms)
    return lambda t: expr

def get_curve_point_deltas(points, multiplier):
    deltas = []
    for i in range(0, len(points)-1):
        deltas.append(multiplier * (points[i+1] - points[i]))
    return deltas


def patch_pos_3d(patch, u, v):
    degree = len(patch[0])-1
    bases_u = bezier_bases(degree, u)
    bases_v = bezier_bases(degree, v)

    pos = Matrix([0, 0, 0])

    for i in range(0, degree+1):
        for j in range(0, degree+1):
            pos += patch[i][j] * bases_u[i] * bases_v[j]

    return pos


def quadratic_patch_pos_3d(patch, u, v):
    bases_u = bezier_bases(2, u)
    bases_v = bezier_bases(2, v)

    pos = Matrix([0, 0, 0])

    for i in range(0, 3):
        for j in range(0, 3):
            pos += patch[i][j] * bases_u[i] * bases_v[j]

    return pos


def quadratic_patch_normal_3d(patch_d, u, v):
    bases_u = bezier_bases(1, u)
    bases_v = bezier_bases(1, v)

    normal = Matrix([0, 0, 0])

    for i in range(0, 2):
        for j in range(0, 2):
            tangents = patch_d[i][j]
            normal += tangents[0].cross(tangents[1]) * \
                (bases_u[i] * bases_v[j])

    return normal
