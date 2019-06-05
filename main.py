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


def bezier_bases(n, param):
    bases = []
    for i in range(0, n+1):
        bases.append(bernstein_basis(n, i, param))

    return bases

def bernstein_basis(n, i, param):
    basis = binomial(n, i) * param**i * (1 - param)**(n-i)

    return basis

# straight port from Mathematica snippet by Pomax
# https://stackoverflow.com/questions/35901079/calculating-the-inflection-point-of-a-cubic-bezier-curve
def bezier_curvature():
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

    bx1 = lambda t : a * (1 - t)**2 + 2 * b * (1 - t) * t + c * t**2
    bx2 = lambda t : u * (1 - t) + v * t
    by1 = lambda t : d * (1 - t)**2 + 2 * e * (1 - t) * t + f * t**2
    by2 = lambda t : w * (1 - t) + z * t

    curvature = bx1(t) * by2(t) - by1(t) * bx2(t)

    result = expand(curvature)

    return (symbs, result)

def bezier_curvature_new():
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

    
    bases_2nd = bezier_bases(2, t)
    bases_3rd = bezier_bases(1, t)

    points_x_2nd = (a, b, c)
    points_x_3rd = (u, v)
    points_y_2nd = (d, e, f)
    points_y_3rd = (w, z)

    bx1 = produce_bezier(points_x_2nd, bases_2nd)
    bx2 = produce_bezier(points_x_3rd, bases_3rd)
    by1 = produce_bezier(points_y_2nd, bases_2nd)
    by2 = produce_bezier(points_y_3rd, bases_3rd)

    curvature = bx1(t) * by2(t) - by1(t) * bx2(t)

    result = expand(curvature)

    return (symbs, result)


def produce_bezier(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i"%(len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    bezier = reduce((lambda x, y: x + y), terms)
    func = lambda param: bezier
    return func

def simplify_curvature(symbs, expr):
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


def geometry():
    N = CoordSys3D('N')
    v = N.i * 3 + N.j * 2 + N.k * 1
    print(v)

def main():
    symbs, expr = bezier_curvature_new();
    symbs, expr = simplify_curvature(symbs, expr)

    symbs2, expr2 = bezier_curvature()
    symbs2, expr2 = simplify_curvature(symbs2, expr2)

    are_equal = simplify(expr - expr2) == 0
    pprint("Old and new equal? " + str(are_equal))

    # pprint(expr)
    # print("\nNumber of terms: " + str(len(expr.args)))

    # symbs, expr, subs = cache_variables(symbs, expr)
    # pprint(expr)

    # print(to_latex_docstring(expr))
    # # show_expr_latex(expr)

    # geometry()

    # bases = bezier_bases(3)
    # for b in bases:
    #     pprint(b)

    # a, b, c = symbols('a b c')
    # points = (a, b, c)
    # bases = bezier_bases(2)
    # bezier_terms = [p * b for p, b in zip(points, bases)]
    
    # for t in bezier_terms:
    #     pprint(t)

    # bezier = reduce((lambda x, y: x + y), bezier_terms)
    # pprint(bezier)
    

    
if __name__ == "__main__":
    main()
