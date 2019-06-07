from sympy import *
from sympy.vector import *
import matplotlib.pyplot as plt
import math
from functools import reduce


def invert_dict(dict):
    return {v: k for k, v in dict.items()}

# convert sympy expression to inline latex ready for showing in MatPlotLib
def show_expr_latex(expr):
    show_latex_str(latex(expr, mode='inline'))

    # https://docs.sympy.org/latest/modules/printing.html#sympy.printing.latex.latex

# Render a latex string using MatPlotLib
def show_latex_str(latex_str):
    plt.rc('text', usetex=True)
    plt.title(latex_str)
    plt.show()

    # https://matplotlib.org/users/usetex.html
    # Note: requires working install of LateX on PATH

def to_latex_docstring(expr):
    doc_srt = r"""
\documentclass{article}
\begin{document}
    %s
\end{document}""" % (
        latex(expr, mode='inline')
    )

    return doc_srt

# Make a 3d point, with uniquely named scalars name_x, name_y, name_z, from frame N
def make_point(N, name):
    coords = ('x', 'y', 'z')
    coords = map(lambda c: name + "_" + c + " ", coords)
    coords = reduce(lambda a, b: a + b, coords)
    x, y, z = symbols(coords)
    return x*N.i + y*N.j + z*N.k

def bezier_bases(n, param):
    bases = []
    for i in range(0, n+1):
        bases.append(bernstein_basis(n, i, param))

    return bases

def bernstein_basis(n, i, param):
    basis = binomial(n, i) * param**i * (1 - param)**(n-i)

    return basis


def make_bezier_expr(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i" % (
            len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    expr = reduce((lambda x, y: x + y), terms)
    return lambda t: expr

def bezier_curvature_2d():
    symbs = symbols('t x1 x2 x3 x4 y1 y2 y3 y4')
    t, x1, x2, x3, x4, y1, y2, y3, y4 = symbs

    a = 3 * (x2 - x1)
    b = 3 * (x3 - x2)
    c = 3 * (x4 - x3)
    u = 2 * (b - a)
    v = 2 * (c - b)

    d = 3 * (y2 - y1)
    e = 3 * (y3 - y2)
    f = 3 * (y4 - y3)
    w = 2 * (e - d)
    z = 2 * (f - e)
    
    bases_1st_deriv = bezier_bases(2, t)
    bases_2nd_deriv = bezier_bases(1, t)

    points_x_1st = (a, b, c)
    points_x_2nd = (u, v)
    points_y_1st = (d, e, f)
    points_y_2nd = (w, z)

    bx1 = make_bezier_expr(points_x_1st, bases_1st_deriv)
    bx2 = make_bezier_expr(points_x_2nd, bases_2nd_deriv)
    by1 = make_bezier_expr(points_y_1st, bases_1st_deriv)
    by2 = make_bezier_expr(points_y_2nd, bases_2nd_deriv)

    curvature = bx1(t) * by2(t) - by1(t) * bx2(t)
    result = expand(curvature)
    return (symbs, result)

# Assume a given cubic curve starts at [0,0] and ends at [x,0]
# leads to x1, y1, y2 = 0
def to_oriented_curve_2d(expr):
    x1, y1, y4 = symbols('x1, y1, y4')
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(y4, 0)

    return (expr_subbed)


def to_oriented_curve_3d(expr):
    x1, y1, z1, y4, z4 = symbols('x1, y1, z1, y4, z4')
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(z1, 0) \
        .subs(y4, 0) \
        .subs(z4, 0) \

    return (expr_subbed)

# Replace repeaded terms with variable names to emphasize their cachable nature
#
# todo: This is redundant, could probably use cse()
# https://docs.sympy.org/latest/modules/rewriting.html
#
# or: at least generalize by matching scalar * scalar
# https: // docs.sympy.org/latest/modules/utilities/iterables.html
def cache_variables(symbs, expr):
    t, x1, x2, x3, x4, y1, y2, y3, y4 = symbs

    sub_symbs = symbols('a b c d')
    a, b, c, d = sub_symbs

    substitutions = {
        a: x3 * y2,
        b: x4 * y2,
        c: x2 * y3,
        d: x4 * y3
    }

    substitutions_inv = { v: k for k, v in substitutions.items() }

    expr_subbed = expr.subs(substitutions_inv)

    new_symbs = symbs + sub_symbs
    return (new_symbs, expr_subbed, substitutions)

def solve_quadratic(expr, t):
    expr = simplify(expr)
    expr = collect(expr, t)  # collect in terms of a*t^0, b*t^1, c*t^2, ...

    poly = Poly(expr, t)
    coeffs = poly.coeffs()
    print("Got polynomial of degree: " + str(poly.degree()))

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


def print_pretty(common, exprs):
    print("\n---------------terms-----------------\n")

    for t in common:
        pprint(t)

    print("\n----------------roots-------------------\n")

    for expr in exprs:
        pprint(expr)

def print_code(common, exprs):
    print("\n----------------terms-------------------\n")

    for t in common:
        code = ccode(t)
        code = code[1:-1]
        comma_idx = code.find(",")
        code = code[0:comma_idx] + " =" + code[comma_idx+1:]
        code = code.replace("pow", "math.pow")
        code = code.replace("sqrt", "math.sqrt")
        code = "float " + code + ";"
        print(code)

    print("\n----------------roots-------------------\n")

    for i, expr in enumerate(exprs):
        print("float root_%d = " % i + ccode(expr) + ";")

# Todo: really want to use linear algebra abstractions here
# more importantly: use proper formulations of curvature, lol
# https://en.wikipedia.org/wiki/Curvature (Space Curves)
# alternative: https://en.wikipedia.org/wiki/Radius_of_curvature
# radius of curvature and curvature both defined by arclenght, bah
#
# Also: SIGNED curvature. Otherwise we don't get the actual
# sign flip that we're looking for XD


def bezier_inflections_3d():
    t = symbols('t')

    # Using the below, I get quadratics

    # symbs = symbols('t x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4')
    # t, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4 = symbs

    # xd1 = 3 * (x2 - x1)
    # xd2 = 3 * (x3 - x2)
    # xd3 = 3 * (x4 - x3)
    # xdd1 = 2 * (xd2 - xd1)
    # xdd2 = 2 * (xd3 - xd2)

    # yd1 = 3 * (y2 - y1)
    # yd2 = 3 * (y3 - y2)
    # yd3 = 3 * (y4 - y3)
    # ydd1 = 2 * (yd2 - yd1)
    # ydd2 = 2 * (yd3 - yd2)

    # zd1 = 3 * (z2 - z1)
    # zd2 = 3 * (z3 - z2)
    # zd3 = 3 * (z4 - z3)
    # zdd1 = 2 * (zd2 - zd1)
    # zdd2 = 2 * (zd3 - zd2)


    # using these though, I get cubics. Why? Something in the above that cancels automatically?
    # I guess...
    # Idea: after performing this expansion, sub the intermediate terms by x1..z4, and eliminate
    # the zero terms. Will not change order of curve, but still simplify the resulting equation.

    symbs_d = symbols('xd1 xd2 xd3 yd1 yd2 yd3 zd1 zd2 zd3')
    symbs_dd = symbols('xdd1 xdd2 ydd1 ydd2 zdd1 zdd2')
    xd1, xd2, xd3, yd1, yd2, yd3, zd1, zd2, zd3, = symbs_d
    xdd1, xdd2, ydd1, ydd2, zdd1, zdd2 = symbs_dd

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    points_x_d = (xd1, xd2, xd3)
    points_y_d = (yd1, yd2, yd3)
    points_z_d = (zd1, zd2, zd3)

    points_x_dd = (xdd1, xdd2)
    points_y_dd = (ydd1, ydd2)
    points_z_dd = (zdd1, zdd2)

    xd = make_bezier_expr(points_x_d, bases_d)
    yd = make_bezier_expr(points_y_d, bases_d)
    zd = make_bezier_expr(points_z_d, bases_d)
    xdd = make_bezier_expr(points_x_dd, bases_dd)
    ydd = make_bezier_expr(points_y_dd, bases_dd)
    zdd = make_bezier_expr(points_z_dd, bases_dd)

    # this gives us 3 polys, for potentially 6 roots

    # N = CoordSys3D('N')
    # b1 = bx1(t) * N.i + by1(t) * N.j + bz1(t) * N.j
    # b2 = bx2(t) * N.i + by2(t) * N.j + bz2(t) * N.j
    # curvature = cross(b2, b1)

    # store resulting coefficients for each spatial basis separately
    result = [
        yd(t) * zdd(t) - zd(t) * ydd(t),
        zd(t) * xdd(t) - xd(t) * zdd(t),
        xd(t) * ydd(t) - yd(t) * xdd(t)
    ]

    return (t, result)

def substitute_coeffs(expr):
    symbs = symbols('t x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4')
    symbs_d = symbols('xd1 xd2 xd3 yd1 yd2 yd3 zd1 zd2 zd3')
    symbs_dd = symbols('xdd1 xdd2 ydd1 ydd2 zdd1 zdd2')
    t, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4 = symbs
    xd1, xd2, xd3, yd1, yd2, yd3, zd1, zd2, zd3, = symbs_d
    xdd1, xdd2, ydd1, ydd2, zdd1, zdd2 = symbs_dd

    substitutions = {
        xd1: 3 * (x2 - x1),
        xd2: 3 * (x3 - x2),
        xd3: 3 * (x4 - x3),
        yd1: 3 * (y2 - y1),
        yd2: 3 * (y3 - y2),
        yd3: 3 * (y4 - y3),
        zd1: 3 * (z2 - z1),
        zd2: 3 * (z3 - z2),
        zd3: 3 * (z4 - z3),
        xdd1: 2 * (xd2 - xd1),
        xdd2: 2 * (xd3 - xd2),
        ydd1: 2 * (yd2 - yd1),
        ydd2: 2 * (yd3 - yd2),
        zdd1: 2 * (zd2 - zd1),
        zdd2: 2 * (zd3 - zd2)
    }

    return expr.subs(substitutions)

def curvature_3d():
    # Todo: formulate entirely using linear algebra
    # will simplify the code, shrink it.
    #
    # Also, optimize terms to yield efficient operations
    # on vector quantities, as that will work well with
    # SIMD and whatnot

    t, exprs = bezier_inflections_3d()

    for i in range(0, len(exprs)):
        # recursively replace xddi -> xdi -> xi
        for l in range(0, 2):
            exprs[i] = substitute_coeffs(exprs[i])

        exprs[i] = to_oriented_curve_3d(exprs[i])
        exprs[i] = simplify(exprs[i])

    # todo: IF WE GET CUBICS HERE, our formulation of the
    # solution is wrong. Print error.

    a, b = solve_quadratic(exprs[0], t)
    c, d = solve_quadratic(exprs[1], t)
    e, f = solve_quadratic(exprs[2], t)

    common, exprs = cse([a, b, c, d, e, f], numbered_symbols('a'))
    
    print_pretty(common, exprs)
    print_code(common, exprs)

def curvature_2d():
    symbs, expr = bezier_curvature_2d()

    t = symbs[0]

    symbs, expr = to_oriented_curve_2d(expr)
    symbs, expr, subst = cache_variables(symbs, expr)  # substitute with a,b,c,d
    # pprint(expr)

    expr = collect(expr, t)  # collect in terms of a*t^0, b*t^1, c*t^2, ...

    # expr = factor(expr) # factor out common 18
    # todo store factor, and remove it from expr temporarily
    # in a way that works with the polynomical coefficients below

    a, b = solve_quadratic(expr, t)
    common, exprs = cse([a, b], numbered_symbols('a'))
    print_code(common, exprs)

def main():
    curvature_3d()
    

if __name__ == "__main__":
    main()
