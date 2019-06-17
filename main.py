from sympy import *
from sympy.physics.vector import *
import matplotlib.pyplot as plt
import math
from functools import reduce

'''
    Todo:
    - Try a factoring approach to finding some solutions. Since Beziers are built up
    out of Berstein bases, and a lot of arithmetic is adding and multiplying those,
    I wouldn't be surprised if we can find at least some solutions through factoring,
    which would be far nicer than plugging everything into the quadratic formula.
'''


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

# Assume a given cubic curve starts at [0,0] and ends at [x,0]
# leads to x1, y1, y2 = 0
def to_oriented_curve_2d(expr):
    x1, y1, y4 = symbols('x1, y1, y4')
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(y4, 0)

    return (expr_subbed)


def to_oriented_cubic_curve_3d_xyz(expr):
    x1, y1, z1, y4, z4 = symbols('x1, y1, z1, y4, z4')
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(z1, 0) \
        .subs(y4, 0) \
        .subs(z4, 0) \

    return (expr_subbed)


def set_to_zero(symbs, expr):
    subs = { k: v for k, v in map(lambda s: [s, 0], symbs) }
    return expr.subs(subs)

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

def to_polynomial(expr, param):
    expr = simplify(expr)
    expr = collect(expr, param)  # collect in terms of a*t^0, b*t^1, c*t^2, ...
    poly = Poly(expr, param)
    return poly

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
        print(csharp(t))

    print("\n----------------roots-------------------\n")

    for i, expr in enumerate(exprs):
        print("float root_%d = " % i + replace_vector_vars(ccode(expr)) + ";")

'''
Poor man's CodeGen: take ccode output, and format it to valid
C# through very dodge string manipulation
'''

def csharp(code):
    code = ccode(code)
    code = code[1:-1]
    comma_idx = code.find(",")
    code = code[0:comma_idx] + " =" + code[comma_idx+1:]
    code = code.replace("pow", "math.pow")
    code = code.replace("sqrt", "math.sqrt")
    code = "float " + code + ";"

    code = replace_vector_vars(code)
    code = format_floats(code)

    return code


def replace_vector_vars(code):
    '''
    Scan through string finding occurances of 'p*_*'
    Replace each with 'curve[*].*'
    '''
    pos = 0
    while True:
        pos = code.find('p', pos)
        if pos == -1:
            break

        if code[pos+2] == '_':
            
            idx = int(code[pos+1])
            idx -= 1
            code = code[0:pos] + "curve[" + str(idx) + "]." + code[pos+3:]
        else:
            pos += 1
    return code


def format_floats(code):
    '''
    Scan through string finding occurances of 'n.n'
    Replace each with 'n.nf'
    '''
    pos = 0
    while True:
        pos = code.find('.', pos)
        if pos == -1:
            break

        if code[pos-1].isdigit() and code[pos+1].isdigit() and code[pos+2] != 'f':
            code = code[0:pos+2] + "f" + code[pos+2:]
        else:
            pos += 1
    return code

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

# Todo: really want to use linear algebra abstractions here
# more importantly: use proper formulations of curvature, lol
# https://en.wikipedia.org/wiki/Curvature (Space Curves)
# alternative: https://en.wikipedia.org/wiki/Radius_of_curvature
# radius of curvature and curvature both defined by arclenght, bah
#
# Also: SIGNED curvature. Otherwise we don't get the actual
# sign flip that we're looking for XD


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
 
    xdd = make_bezier_expr(points_x_dd, bases_dd)
    ydd = make_bezier_expr(points_y_dd, bases_dd)
    zdd = make_bezier_expr(points_z_dd, bases_dd)

    result = [
        xdd(t),
        ydd(t),
        zdd(t)
    ]

    return (t, result)

def bezier_height_dt_3d():
    t = symbols('t')

    symbs_d = symbols('xd1 xd2 xd3 yd1 yd2 yd3 zd1 zd2 zd3')
    xd1, xd2, xd3, yd1, yd2, yd3, zd1, zd2, zd3, = symbs_d

    bases_d = bezier_bases(2, t)

    points_y_d = (yd1, yd2, yd3)

    yd = make_bezier_expr(points_y_d, bases_d)

    return (t, yd(t))

def inflections_3d():
    # Todo: formulate entirely using linear algebra
    # will simplify the code, shrink it.
    #
    # Also, optimize terms to yield efficient operations
    # on vector quantities, as that will work well with
    # SIMD and whatnot

    t, exprs = bezier_curvature_partials_3d()

    for i in range(0, len(exprs)):
        exprs[i] = substitute_coeffs(exprs[i])
        exprs[i] = to_oriented_cubic_curve_3d_xyz(exprs[i])
        exprs[i] = simplify(exprs[i])
    
    # todo: IF WE GET CUBICS HERE, our formulation of the
    # solution is wrong. Print error.

    # a, b = solve_quadratic(exprs[0], t)
    # c, d = solve_quadratic(exprs[1], t)
    # e, f = solve_quadratic(exprs[2], t)
    # common, exprs = cse([a,b,c,d,e,f], numbered_symbols('a'))

    a = solveset(exprs[0], t)
    b = solveset(exprs[1], t)
    c = solveset(exprs[2], t)
    common, exprs = cse([a,b,c], numbered_symbols('a'))
    
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


def maxima_2nd_cubic_2d():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')

    points = [p1, p2, p3, p4]
    points_d = get_curve_point_deltas(points, 3)
    points_dd = get_curve_point_deltas(points_d, 2)

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    pdd = make_bezier_expr(points_dd, bases_dd)(t)
    pdd = expand(pdd)

    solutions = map(lambda partial: solveset(partial,t).args[0], pdd)

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
    points_d = get_curve_point_deltas(points, 3)
    points_dd = get_curve_point_deltas(points_d, 2)

    bases_d = bezier_bases(2, t)
    bases_dd = bezier_bases(1, t)

    pd = make_bezier_expr(points_d, bases_d)
    pdd = make_bezier_expr(points_dd, bases_dd)

    curvature = cross_2d(pd(t), pdd(t))
    curvature = expand(curvature)

    a = solveset(Eq(curvature,0), t)

    common, exprs = cse(a, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)

def get_curve_point_deltas(points, multiplier):
    deltas = []
    for i in range(0, len(points)-1):
        deltas.append(multiplier * (points[i+1] - points[i]))
    return deltas

def cross_2d(p1, p2):
    return p1[0] * p2[1] - p1[1] * p2[0]

def silhouette_cubic_2d():
    t = symbols('t')

    vx = 0
    vy = 0

    symbs = symbols('x1 x2 x3 x4 y1 y2 y3 y4')
    x1, x2, x3, x4, y1, y2, y3, y4 = symbs

    xd1 = 3 * (x2 - x1)
    xd2 = 3 * (x3 - x2)
    xd3 = 3 * (x4 - x3)
   
    yd1 = 3 * (y2 - y1)
    yd2 = 3 * (y3 - y2)
    yd3 = 3 * (y4 - y3)

    bases = bezier_bases(3, t)
    bases_d = bezier_bases(2, t)

    points_x = (x1, x2, x3, x4)
    points_y = (y1, y2, y3, y4)
    points_x_d = (xd1, xd2, xd3)
    points_y_d = (yd1, yd2, yd3)

    x = make_bezier_expr(points_x, bases)(t)
    y = make_bezier_expr(points_y, bases)(t)
    xd = make_bezier_expr(points_x_d, bases_d)(t)
    yd = make_bezier_expr(points_y_d, bases_d)(t)

    normal_x = -yd
    normal_y = xd

    viewdir_x = x - vx
    viewdir_y = y - vy

    solution = viewdir_x * normal_x + viewdir_y * normal_y
    solution = expand(solution)
    solution = to_oriented_cubic_curve_3d_xyz(solution)

    poly = to_polynomial(solution, t)
    print("Got polynomial of degree: " + str(poly.degree()))

    pprint(solution)

    # solution = solveset(solution, t)
    # common, exprs = cse(solution, numbered_symbols('a'))

    # print_code(common, exprs)

def silhouette_quadratic_2d():
    # The setup:

    # a curve (const), we need to express a point at t, and its normal
    # a view point (const)
    # direction from view point to curve point
    # dot product of normal and view direction
    # find t where that dot product = 0

    t = symbols('t')

    view_point = symbolic_vector_2d('viewPoint')
    # view_point = Matrix([0,0])

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')

    p1 = Matrix([0, 0])

    pd1 = 3 * (p2 - p1)
    pd2 = 3 * (p3 - p2)

    bases = bezier_bases(2, t)
    bases_d = bezier_bases(1, t)

    points = (p1, p2, p3)
    points_d = (pd1, pd2)

    p = make_bezier_expr(points, bases)(t)
    pd = make_bezier_expr(points_d, bases_d)(t)

    normal = Matrix([-pd[1], pd[0]])
    viewdir = p - view_point

    solution = viewdir.dot(normal)
    solution = expand(solution)

    poly = to_polynomial(solution, t)
    print("Got polynomial of degree: " + str(poly.degree()))

    solution = solveset(solution, t)
    common, exprs = cse(solution, numbered_symbols('a'))
    print_code(common, exprs)

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
            normal += tangents[0].cross(tangents[1]) * (bases_u[i] * bases_v[j])

    return normal


'''
Find silhouette of a square, quadratic bezier patch;
which is a net of 9 points.

--

Important: We want to find a silhouette CURVE, not
a poorly defined locus of silhouette points.

Therefore, grad(dot(viewdir(p(u,v)), n(u,v))) will not do,
as we would be left clueless as to which points to actually
construct our silhouette curve from.

--

We know that, supposing solution is a polynomial curve defined
by some p1, p2, ..., pn, that p1 and pn lie on patch's edge curves,
parameterized by some scalar params edge1_v, edge2_v.

It would therefore suffice to look exclusively along expected
patch edge for each of those, right?

Idea:

class square_patch_quadratic_3d
    9 control points
    eval_point(u,v)
    eval_tangent(u,v)
    eval_normal(u,v)

With u or v = 0, or 1, the equations would automatically
simplify, because of the p*t and p*(1-t) terms

Can then process the 4 edge curves using a subroutine:
eval_normal(u,0)
eval_normal(u,1)
eval_normal(0,v)
eval_normal(1,v)

Result: works, but we do get cubic solutions...

New idea: project view_pos, view-dir onto curve plane,
solve in 2d. Actual code would first need to check if
this is possible, with an early-rejection test.

--

Possible Optimizations
- Align view_point or p1 to [0,0,0]
- Align last curve point to [x,0,0]
- Align a prominent tangent to [x,0,0]
- Align whole quadratic curve to basis plane
- Solve 2d sub problems, instead of general 3d ones
etc.
'''
def silhouette_quadratic_patch_3d():
    u, v = symbols('u v')

    # bottom edge only
    v = 0

    view_point = symbolic_vector_3d('pview')
    # view_point = Matrix([0, 0, 0])
    view_point[2] = 0

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')],
    ]

    # patch_d = [
    #     [[symbolic_vector_3d('q1u'), symbolic_vector_3d('q1v')], [symbolic_vector_3d('q2u'), symbolic_vector_3d('q2v')]],
    #     [[symbolic_vector_3d('q3u'), symbolic_vector_3d('q3v')], [symbolic_vector_3d('q4u'), symbolic_vector_3d('q4v')]]
    # ]

    # todo: generate these derivative points, q, using a function
    patch_d = [
        [[3 * (patch[0][1] - patch[1][1]), 3 * (patch[0][1] - patch[0][2])], [3 * (patch[2][1] - patch[1][1]), 3 * (patch[1][2] - patch[1][1])]],
        [[3 * (patch[1][0] - patch[0][0]), 3 * (patch[0][1] - patch[0][0])], [3 * (patch[2][0] - patch[1][0]), 3 * (patch[1][1] - patch[1][0])]],
    ]

    pos = quadratic_patch_pos_3d(patch, u, v)
    normal = quadratic_patch_normal_3d(patch_d, u, v)

    viewdir = pos - view_point

    solution = viewdir.dot(normal)
    solution = expand(solution)

    # Todo: something like this, but for each separate edge curve
    # also not that we can lay any quadratic float on a basis plane, which
    # is even better
    solution = set_to_zero(symbols('p1_x p1_y p1_z p2_z p3_y p3_z'), solution)

    # pprint(solution)

    # poly = to_polynomial(solution, u)
    # print("Got polynomial of degree: " + str(poly.degree()))

    solution = solveset(solution, u)
    common, exprs = cse(solution, numbered_symbols('a'))

    # print_pretty(common, exprs)
    print_code(common, exprs)

def quadratic_2d_bezier():
    symbs = symbols('t, p1, p2, p3')
    t, p1, p2, p3 = symbs

    pd1 = 2 * (p2 - p1)
    pd2 = 2 * (p3 - p2)

    bases = bezier_bases(2, t)
    bases_d = bezier_bases(1, t)

    p = make_bezier_expr((p1, p2, p3), bases)(t)
    pd = make_bezier_expr((pd1, pd2), bases_d)(t)

    common, exprs = cse(p, numbered_symbols('a'))
    print("Point:")
    print_code(common, exprs)

    common, exprs = cse(pd, numbered_symbols('a'))
    print("Tangent:")
    print_code(common, exprs)

# def test_symbol_replacement():
#     p1 = symbolic_vector_3d("p1")
#     print("Before: ")
#     pprint(p1)
#     print("After: ")
#     p1 = set_to_zero(p1)
#     pprint(p1)

def main():
    # inflections_3d()
    # curvature_maxima_3d()    
    # height_maxima_3d()

    # silhouette_cubic_2d()
    # silhouette_quadratic_2d()
    # silhouette_quadratic_patch_3d()

    # quadratic_2d_bezier()

    maxima_2nd_cubic_2d()
    # inflections_cubic_2d()


if __name__ == "__main__":
    main()
