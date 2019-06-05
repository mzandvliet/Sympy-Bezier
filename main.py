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


def make_bezier_expr(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i"%(len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    expr = reduce((lambda x, y: x + y), terms)
    return lambda t: expr

# Assume a given cubic curve starts at [0,0] and ends at [x,0]
def simplify_curvature_2d(symbs, expr):
    t, x1, x2, x3, x4, y1, y2, y3, y4 = symbs
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(y4, 0)

    return (symbs, expr_subbed)


# Replace repeaded terms with variable names to emphasize their cachable nature
#
# todo: automatically determine repeated terms like these and
# generate a, b, c... to cache them in
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


def substitute(expr, substitutions):
    substitutions_inv = {v: k for k, v in substitutions.items()}
    expr_subbed = expr.subs(substitutions_inv)
    return expr_subbed


def main():
    symbs, expr = bezier_curvature_2d()

    t = symbs[0]

    symbs, expr = simplify_curvature_2d(symbs, expr)
    symbs, expr, subst = cache_variables(symbs, expr) # substitute with a,b,c,d
    expr = collect(expr, t)  # collect in terms of a*t^0, b*t^1, c*t^2, ...

    # expr = factor(expr) # factor out common 18
    # todo store factor, and remove it from expr temporarily
    # in a way that works with the polynomical coefficients below

    pprint(expr)
    print("----")

    poly = Poly(expr, t)
    coeffs = poly.coeffs()
    for c in enumerate(coeffs):
        pprint(c)

    v1, v2, v3 = symbols('v1 v2 v3')
    substitutions = {
        v1: coeffs[0],
        v2: coeffs[1],
        v3: coeffs[2],
    }
    
    expr_v = substitute(expr, substitutions)

    pprint(expr_v)

    root = solve(expr_v, t)
    pprint(root)

    # inflection = solve(Eq(expr, 0), symbs[0])
    # pprint(inflection)

    # symbs, expr, substitutions = cache_variables(symbs, expr)
    

    # pprint(expr)
    # print("\nNumber of terms: " + str(len(expr.args)))

    # show_expr_latex(inflection)
    

    
if __name__ == "__main__":
    main()
