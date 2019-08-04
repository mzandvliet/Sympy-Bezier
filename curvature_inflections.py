import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *


def prove_curve_derives():
    t = symbols('t')

    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')

    points = (p1, p2, p3)
    pos = make_bezier(points, bezier_bases(2, t))(t)

    tangents = differentiate_curve_points(points)
    tangents_a = expand(make_bezier(tangents, bezier_bases(1, t))(t))

    tangents_b = diff(pos, t)

    tangents_a = expand(tangents_a)
    tangents_b = expand(tangents_b)

    print("\nTangents A:\n")
    pprint(tangents_a)
    print("\nTangents B:\n")
    pprint(tangents_b)
    print("\nDifference:\n")
    pprint(tangents_b - tangents_a)

# Assume a given cubic curve starts at [0,0] and ends at [x,0]
# leads to x1, y1, y2 = 0

def to_oriented_curve_2d(expr):
    x1, y1, y4 = symbols('x1, y1, y4')
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(y4, 0)

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

    substitutions_inv = {v: k for k, v in substitutions.items()}

    expr_subbed = expr.subs(substitutions_inv)

    new_symbs = symbs + sub_symbs
    return (new_symbs, expr_subbed, substitutions)


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

    # recursively replace xddi -> xdi -> xi
    for l in range(0, 2):
        expr = expr.subs(substitutions)

    return expr


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

    bx1 = make_bezier(points_x_1st, bases_1st_deriv)
    bx2 = make_bezier(points_x_2nd, bases_2nd_deriv)
    by1 = make_bezier(points_y_1st, bases_1st_deriv)
    by2 = make_bezier(points_y_2nd, bases_2nd_deriv)

    curvature = bx1(t) * by2(t) - by1(t) * bx2(t)
    result = expand(curvature)
    return (symbs, result)


def bezier_curvature_partials_3d():
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

    xd = make_bezier(points_x_d, bases_d)
    yd = make_bezier(points_y_d, bases_d)
    zd = make_bezier(points_z_d, bases_d)
    xdd = make_bezier(points_x_dd, bases_dd)
    ydd = make_bezier(points_y_dd, bases_dd)
    zdd = make_bezier(points_z_dd, bases_dd)

    # this gives us 3 polys, for potentially 6 roots

    # N = CoordSys3D('N')
    # b1 = bx1(t) * N.i + by1(t) * N.j + bz1(t) * N.j
    # b2 = bx2(t) * N.i + by2(t) * N.j + bz2(t) * N.j
    # curvature = cross(b2, b1)

    # store resulting coefficients for each spatial basis separately
    result = [
        ydd(t) * zd(t) - zdd(t) * yd(t),
        zdd(t) * xd(t) - xdd(t) * zd(t),
        xdd(t) * yd(t) - ydd(t) * xd(t)
    ]

    return (t, result)


def bezier_ddt_partials():
    t = symbols('t')

    symbs_dd = symbols('xdd1 xdd2 ydd1 ydd2 zdd1 zdd2')
    xdd1, xdd2, ydd1, ydd2, zdd1, zdd2 = symbs_dd

    bases_dd = bezier_bases(1, t)

    points_x_dd = (xdd1, xdd2)
    points_y_dd = (ydd1, ydd2)
    points_z_dd = (zdd1, zdd2)

    xdd = make_bezier(points_x_dd, bases_dd)
    ydd = make_bezier(points_y_dd, bases_dd)
    zdd = make_bezier(points_z_dd, bases_dd)

    result = [
        xdd(t),
        ydd(t),
        zdd(t)
    ]

    return (t, result)


def inflections_cubic_3d():
    t = symbols('t')

    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')

    points = [p1, p2, p3, p4]
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    pd = make_bezier(points_d, bases_d)
    pdd = make_bezier(points_dd, bases_dd)

    curvature = pd(t).cross(pdd(t))
    curvature = expand(curvature)
    solutions = map(lambda partial: solveset(partial, t).args[0], curvature)

    common, exprs = cse(solutions, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)


def curvature_maxima_3d():
    t, exprs = bezier_curvature_partials_3d()

    for i in range(0, len(exprs)):
        exprs[i] = diff(exprs[i], t)
        exprs[i] = substitute_coeffs(exprs[i])
        # exprs[i] = to_oriented_curve_3d(exprs[i])
        exprs[i] = simplify(exprs[i])

    a = solveset(exprs[0], t)
    b = solveset(exprs[1], t)
    c = solveset(exprs[2], t)
    common, exprs = cse([a, b, c], numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)


def maxima_1st_cubic_2d():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')

    points = [p1, p2, p3, p4]
    points_d = differentiate_curve_points(points)

    bases_d = bezier_bases(2, t)

    pd = make_bezier(points_d, bases_d)(t)
    pd = expand(pd)

    solutions = map(lambda partial: solveset(partial, t).args[0], pd)

    common, exprs = cse(solutions, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)


def maxima_2nd_cubic_2d():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')

    points = [p1, p2, p3, p4]
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    bases_dd = bezier_bases(1, t)

    pdd = make_bezier(points_dd, bases_dd)(t)
    pdd = expand(pdd)

    solutions = map(lambda partial: solveset(partial, t).args[0], pdd)

    common, exprs = cse(solutions, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)


def inflections_cubic_2d():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')

    points = [p1, p2, p3, p4]
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    pd = make_bezier(points_d, bases_d)
    pdd = make_bezier(points_dd, bases_dd)

    curvature = cross_2d(pd(t), pdd(t))
    curvature = expand(curvature)

    a = solveset(Eq(curvature, 0), t)

    common, exprs = cse(a, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)


def inflections_deriv_cubic_2d():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')

    points = [p1, p2, p3, p4]
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    pd = make_bezier(points_d, bases_d)
    pdd = make_bezier(points_dd, bases_dd)

    curvature = cross_2d(pd(t), pdd(t))
    curvature = expand(curvature)

    curvature_deriv = diff(curvature, t)

    a = solveset(Eq(curvature_deriv, 0), t)

    common, exprs = cse(a, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)


def curvature_maxima_cubic_2d():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')

    p1 = Matrix([0, 0])
    p4[1] = 0

    points = [p1, p2, p3, p4]
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    pd = make_bezier(points_d, bases_d)
    pdd = make_bezier(points_dd, bases_dd)

    # Full curvature definition: https://pomax.github.io/bezierinfo/#curvature

    curvature = cross_2d(pd(t), pdd(t)) / \
        (pd(t)[0]**2 + pd(t)[1]**2)**Rational(3/2)
    curvature = simplify(curvature)
    curv_dif = diff(curvature, t)

    # this results in an absolutely soul-destroying mountain equation. Nevermind :P
