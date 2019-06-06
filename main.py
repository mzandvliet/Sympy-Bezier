from sympy import *
from sympy.vector import *
import matplotlib.pyplot as plt
import math
from functools import reduce

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

# Todo: really want to use linear algebra abstractions here
def bezier_curvature_3d():
    symbs = symbols('t x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4')
    t, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4 = symbs

    xd1 = 3 * (x2 - x1)
    xd2 = 3 * (x3 - x2)
    xd3 = 3 * (x4 - x3)
    xdd1 = 2 * (xd2 - xd1)
    xdd2 = 2 * (xd3 - xd2)

    yd1 = 3 * (y2 - y1)
    yd2 = 3 * (y3 - y2)
    yd3 = 3 * (y4 - y3)
    ydd1 = 2 * (yd2 - yd1)
    ydd2 = 2 * (yd3 - yd2)

    zd1 = 3 * (z2 - z1)
    zd2 = 3 * (z3 - z2)
    zd3 = 3 * (z4 - z3)
    zdd1 = 2 * (zd2 - zd1)
    zdd2 = 2 * (zd3 - zd2)

    bases_1st_deriv = bezier_bases(2, t)
    bases_2nd_deriv = bezier_bases(1, t)

    points_x_1st = (xd1, xd2, xd3)
    points_y_1st = (yd1, yd2, yd3)
    points_z_1st = (zd1, zd2, zd3)

    points_x_2nd = (xdd1, xdd2)
    points_y_2nd = (ydd1, ydd2)
    points_z_2nd = (zdd1, zdd2)

    bx1 = make_bezier_expr(points_x_1st, bases_1st_deriv)
    by1 = make_bezier_expr(points_y_1st, bases_1st_deriv)
    bz1 = make_bezier_expr(points_z_1st, bases_1st_deriv)
    bx2 = make_bezier_expr(points_x_2nd, bases_2nd_deriv)
    by2 = make_bezier_expr(points_y_2nd, bases_2nd_deriv)
    bz2 = make_bezier_expr(points_z_2nd, bases_2nd_deriv)

    N = CoordSys3D('N')

    # 3d cross product
    curvature = \
        (by1(t) * bz2(t) - bz1(t) * by2(t)) * N.i + \
        (bz1(t) * bx2(t) - bx1(t) * bz2(t)) * N.j + \
        (bx1(t) * by2(t) - by1(t) * bx2(t)) * N.k

    curvature = dot(curvature, curvature)

    result = expand(curvature)
    return (symbs, result)


def make_bezier_expr(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i"%(len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    expr = reduce((lambda x, y: x + y), terms)
    return lambda t: expr

# Assume a given cubic curve starts at [0,0] and ends at [x,0]
# leads to x1, y1, y2 = 0
def simplify_curvature_2d(symbs, expr):
    t, x1, x2, x3, x4, y1, y2, y3, y4 = symbs
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(y4, 0)

    return (symbs, expr_subbed)

# Assume a given cubic curve starts at [0,0,0] and ends at [x,0,0]
# leads to x1, y1, y2 = 0
def simplify_curvature_3d(symbs, expr):
    t, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4 = symbs
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(z1, 0) \
        .subs(y4, 0) \
        .subs(z4, 0)

    return (symbs, expr_subbed)


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


def invert_dict(dict):
    return {v: k for k, v in dict.items()}

def curvature_2d():
    symbs, expr = bezier_curvature_2d()

    t = symbs[0]

    symbs, expr = simplify_curvature_2d(symbs, expr)
    symbs, expr, subst = cache_variables(symbs, expr) # substitute with a,b,c,d
    # pprint(expr)

    expr = collect(expr, t)  # collect in terms of a*t^0, b*t^1, c*t^2, ...

    # expr = factor(expr) # factor out common 18
    # todo store factor, and remove it from expr temporarily
    # in a way that works with the polynomical coefficients below

    poly = Poly(expr, t)
    coeffs = poly.coeffs()

    v1, v2, v3 = symbols('v1 v2 v3')
    substitutions = {
        v1: coeffs[0],
        v2: coeffs[1],
        v3: coeffs[2],
    }

    expr_v = expr.subs(invert_dict(substitutions))

    # solve v1*t^2 + v2*t + v3
    roots = solve(expr_v, t)

    root_a = roots[0].subs(substitutions)
    root_b = roots[1].subs(substitutions)

    common, expr = cse([root_a, root_b], numbered_symbols('a'))
    print("---------------formula-----------------")
    pprint(expr)
    print("----------------terms-------------------")
    for t in common:
        pprint(t)

def curvature_3d():
    # Todo: formulate entirely using linear algebra
    # will simplify the code, shrink it.
    #
    # Also, optimize terms to yield efficient operations
    # on vector quantities, as that will work well with
    # SIMD and whatnot


    symbs, expr = bezier_curvature_3d()
    symbs, expr = simplify_curvature_3d(symbs, expr)
    t = symbs[0]
    # terms, expr_optimal = cse(expr)
    
    # pprint(expr_optimal)
    # print("----------------------------------------")
    # pprint(terms)

    expr = collect(expr, t)  # collect in terms of a*t^0, b*t^1, c*t^2, ...

    # expr = factor(expr) # factor out common 18
    # todo store factor, and remove it from expr temporarily
    # in a way that works with the polynomical coefficients below

    poly = Poly(expr, t)
    # coeffs = poly.coeffs()
    print("Got polynomial of degree: " + str(poly.degree()))

    # v1, v2, v3 = symbols('v1 v2 v3')
    # substitutions = {
    #     v1: coeffs[0],
    #     v2: coeffs[1],
    #     v3: coeffs[2],
    # }

    # expr_v = expr.subs(invert_dict(substitutions))

    # # solve v1*t^2 + v2*t + v3
    # roots = solve(expr_v, t)

    # root_a = roots[0].subs(substitutions)
    # root_b = roots[1].subs(substitutions)

    # common, expr = cse([root_a, root_b], numbered_symbols('a'))
    # print("---------------formula-----------------")
    # pprint(expr)
    
    # for t in common:
    #     pprint(t)

    # print("----------------code-------------------")

    # print("root_a = " + ccode(expr[0]) + ";")
    # print("root_b = " + ccode(expr[1]) + ";")

    # print("----------------parts-------------------")

    # for t in common:
    #     code = ccode(t)
    #     print(len(code))


def main():
    curvature_3d()
    

if __name__ == "__main__":
    main()
