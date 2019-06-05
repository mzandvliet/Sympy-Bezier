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

def bezier_curvature_linalg():
    symbs = symbols('t x y z')
    t = symbs

    N = CoordSys3D('N')
    
    # problem: defining 4 points like these makes them all the same, and they cancel to 0 on subtraction
    # p0 = x*N.i + y*N.j + z*N.k
    # p1 = x*N.i + y*N.j + z*N.k
    # p2 = x*N.i + y*N.j + z*N.k
    # p3 = x*N.i + y*N.j + z*N.k

    # but this is too explicit, p0..pn can be symbols without type for the initial expansion, no need to
    # make explicity that they have unique x, y, z
    p0 = make_point(N, 'p0')
    p1 = make_point(N, 'p1')
    p2 = make_point(N, 'p2')
    p3 = make_point(N, 'p3')
    
    a = 3 * (p1-p0)
    b = 3 * (p2-p1)
    c = 3 * (p3-p2)

    print("a: " + str(p0))

    u = 2 * (b-a)
    v = 2 * (c-b)

    bases_2nd = bezier_bases(2, t)
    bases_3rd = bezier_bases(1, t)

    points_2nd = (a, b, c)
    points_3rd = (u, v)

    b1 = produce_bezier(points_2nd, bases_2nd)
    b2 = produce_bezier(points_3rd, bases_3rd)

    curvature = cross(b2, b1)

    for t in b1:
        print(str(t))

    result = expand(curvature)

    return (symbs, result)

# Make a 3d point, with uniquely named scalars name_x, name_y, name_z, from frame N
def make_point(N, name):
    coords = ('x', 'y', 'z')
    coords = map(lambda c: name + "_" + c + " ", coords)
    coords = reduce(lambda a, b: a + b, coords)
    print(coords)
    x, y, z = symbols(coords)
    return x*N.i + y*N.j + z*N.k


def produce_bezier(points, bases):
    if len(points) != len(bases):
        raise Exception("Number of points %i should be equal to number of bases %i"%(len(points), len(bases)))

    terms = [p * b for p, b in zip(points, bases)]
    bezier = reduce((lambda x, y: x + y), terms)
    return bezier

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


def main():
    symbs, expr = bezier_curvature_linalg()
    # symbs, expr = simplify_curvature_2d(symbs, expr)

    pprint(expr)
    print("\nNumber of terms: " + str(len(expr.args)))

    # print(ccode(expr))

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
