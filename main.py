from sympy import *
from sympy.vector import *
import matplotlib.pyplot as plt
import math

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


def bezier_bases(n):
    bases = []
    for i in range(0, n+1):
        bases.append(bernstein_basis(n, i))

    return bases

def bernstein_basis(n, i):
    symbs = symbols('t')
    t = symbs

    basis = binomial(n, i) * t**i * (1 - t)**(n-i)

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
    # symbs, expr = bezier_curvature();
    # symbs, expr = simplify_curvature(symbs, expr)

    # pprint(expr)
    # print("\nNumber of terms: " + str(len(expr.args)))

    # symbs, expr, subs = cache_variables(symbs, expr)
    # pprint(expr)

    # print(to_latex_docstring(expr))
    # # show_expr_latex(expr)

    # geometry()

    bases = bezier_bases(3)
    for b in bases:
        pprint(b)

    
if __name__ == "__main__":
    main()
