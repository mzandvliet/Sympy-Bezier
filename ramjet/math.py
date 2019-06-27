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


def make_bezier(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i" % (
            len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    expr = reduce((lambda x, y: x + y), terms)
    return lambda t: expr

def differentiate_curve_points(points):
    points_dt = []
    input_degree = len(points[0])-1
    for i in range(0, input_degree):
        points_dt.append((input_degree) * (points[i+1] - points[i]))
    return points_dt

def differentiate_patch_points(patch):
    '''
    Todo: this calculates point deltas only,
    but could actually return the full patch
    derivative formulate, just zip with bersteins
    '''

    patch_du = []
    patch_dv = []

    input_degree = len(patch[0])-1
    for v in range(0, input_degree):
        du = []
        dv = []
        for u in range(0, input_degree):
            du.append((input_degree) * (patch[v][u+1] - patch[v][u]))
            dv.append((input_degree) * (patch[v+1][u] - patch[v][u]))

        patch_du.append(du)
        patch_dv.append(dv)

    return (patch_du, patch_dv)


def differentiate_patch_points_u(patch):
    input_degree = len(patch[0])-1
    
    patch_du = []
    for v in range(0, input_degree+1):
        du = []
        for u in range(0, input_degree):
            du.append((input_degree) * (patch[v][u+1] - patch[v][u]))

        patch_du.append(du)

    return patch_du


def differentiate_patch_points_v(patch):
    input_degree = len(patch)-1
    
    patch_dv = []
    for v in range(0, input_degree):
        dv = []
        for u in range(0, input_degree+1):
            dv.append((input_degree) * (patch[v+1][u] - patch[v][u]))

        patch_dv.append(dv)

    return patch_dv

def make_patch(patch, u, v):
    '''
    Given matrix of points and two parameters, constructs a function that
    samples a position along the given surface. Arbitrary degree.
    '''

    degree_u = len(patch[0])-1
    degree_v = len(patch)-1
    bases_u = bezier_bases(degree_u, u)
    bases_v = bezier_bases(degree_v, v)

    pos = Matrix([0, 0, 0])

    for vIdx in range(0, degree_v+1):
        for uIdx in range(0, degree_u+1):
            pos += patch[vIdx][uIdx] * bases_u[uIdx] * bases_v[vIdx]

    return pos
