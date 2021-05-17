import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *

'''
So bicubic filtering is equivalent to interpreting
as a cubic bezier quad.

Quadratic filtering would be the 2nd quad
Bilinear filtering the 1st order quad

Now then, switching to equilateral triangle tesselation,
shouldn't we get that barycentric interpolation also 
extends to higher orders, pulling in more points and
more derivatives?

But, unlike the quad setup, the degree should stay lower?
''' 

def symbol_matrix_4x4():
    mat = []
    for y in range(4):
        row = []
        for x in range(4):
            row.append(symbols('p%i_%i'%(x, y)))
        mat.append(row)
    return mat                

def create_cubic_from_vec4(mat, x):
    return create_cubic(mat[0], mat[1], mat[2], mat[3], x)

def create_cubic(p0, p1, p2, p3, x):
    a = Rational(-1, 2) * p0 + Rational(3, 2) * p1 + Rational(-3, 2) * p2 + Rational(1, 2) * p3
    b = p0 + Rational(-5, 2) * p1 + 2 * p2 + Rational(-1, 2) * p3
    c = Rational(-1, 2) * p0 + Rational(1, 2) * p2
    d = p1

    return a * x**3 + b * x**2 + c*x + d


def cubic_interp():
    print("cubic_interp:\n")

    p0, p1, p2, p3, x = symbols('p0 p1 p2 p3, x')
    cubic = create_cubic(p0, p1, p2, p3, x)
    cubic_dx = diff(cubic, x)
    print("Cubic")
    pprint(cubic)
    print("Cubic_dx")
    pprint(cubic_dx)

    common, exprs = cse(cubic_dx, numbered_symbols('a'))
    print_code_with_vars(common, exprs, (x))

def cubic_interp_neighbors():
    print("cubic_interp_neighbors:\n")

    p0, p1, p2, p3, p4, x = symbols('p0 p1 p2 p3, p4, x')
    cubic_l = create_cubic(p0, p1, p2, p3, x)
    cubic_r = create_cubic(p1, p2, p3, p4, x)

    common, exprs = cse((cubic_l, cubic_r), numbered_symbols('a'))
    print_code_with_vars(common, exprs, (x))

def bicubic_interp():
        # Todo: derive both the polynomial matrix form and the recursive form

    mat = symbol_matrix_4x4()
    u, v = symbols('u v')

    bicubic_u = []
    for i in range(4):
        bicubic_u.append(create_cubic_from_vec4(mat[i], u))

    bicubic = create_cubic_from_vec4(bicubic_u, v)

    pprint(bicubic)

    common, exprs = cse(bicubic, numbered_symbols('a'))
    print_code(common, exprs)

def bicubic_interpolation_tangents():
    mat = symbol_matrix_4x4()

    u, v = symbols('u v')

    bicubic_u = []
    for i in range(4):
        bicubic_u.append(create_cubic_from_vec4(mat[i], u))

    bicubic = create_cubic_from_vec4(bicubic_u, v)

    bicubic_du = diff(bicubic, u)
    bicubic_dv = diff(bicubic, v)

    common, exprs = cse((bicubic_du, bicubic_dv), numbered_symbols('a'))
    print_code(common, exprs)


def bicubic_interpolation_laplacian():
    mat = symbol_matrix_4x4()

    u, v = symbols('u v')

    bicubic_u = []
    for i in range(4):
        bicubic_u.append(create_cubic_from_vec4(mat[i], u))

    bicubic = create_cubic_from_vec4(bicubic_u, v)

    bicubic_du = diff(bicubic, u)
    bicubic_dv = diff(bicubic, v)
    bicubic_ddu = diff(bicubic_du, u)
    bicubic_ddv = diff(bicubic_dv, v)

    laplacian = bicubic_ddu + bicubic_ddv

    common, exprs = cse(laplacian, numbered_symbols('a'))
    print_code(common, exprs)


def bicubic_interpolation_curvature():
    mat = symbol_matrix_4x4()
    pprint(mat)

    u, v = symbols('u v')

    bicubic_u = []
    for i in range(4):
        bicubic_u.append(create_cubic_from_vec4(mat[i], u))

    bicubic = create_cubic_from_vec4(bicubic_u, v)

    du = diff(bicubic, u)
    dv = diff(bicubic, v)

    ddu = diff(bicubic, u, 2)
    ddv = diff(bicubic, v, 2)

    k = (du*ddv - dv*ddu) / ((du*du + dv*dv)**Rational(3,2))

    common, exprs = cse(k, numbered_symbols('a'))
    print_code(common, exprs)


def bicubic_interpolation_combined():
    mat = symbol_matrix_4x4()
    pprint(mat)

    u, v = symbols('u v')

    bicubic_u = []
    for i in range(4):
        bicubic_u.append(create_cubic_from_vec4(mat[i], u))

    c = create_cubic_from_vec4(bicubic_u, v)

    du = diff(c, u)
    dv = diff(c, v)

    ddu = diff(c, u, 2)
    ddv = diff(c, v, 2)

    k = (du*ddv - dv*ddu) / ((du*du + dv*dv)**Rational(3, 2))

    common, exprs = cse((c, du, dv, ddu, ddv, k), numbered_symbols('a'))
    print_code(common, exprs)

'''
The steps where we set first order and second order ends of the curve to zero
are like enforcing boundary conditions.

Thus, we can note that for finding unique solutions to higher-degree polynomials
of this kind, we'll have to supply more boudary conditions. Without them, solutions
become manifold.

What can we say about the spaces of solutions for non-unique variations? Can we still
find effective descriptions of the geometry of the valid space? Can we generate it?
'''

def interp_linear():
    a, b, x = symbols('a b x')
    q = a * x + b

    pprint(q)

    q_0 = q.subs('x', 0)
    q_1 = q.subs('x', 1)

    print("q_0:")
    pprint(q_0)
    print("q_1:")
    pprint(q_1)

    """
    This part is often not explained. Assuming
    we have enough boundary equations set up to
    ensure the solution will be unique, we can
    from the boundary solutions derive equations
    for what the points must be.
    """

    c0 = solve(q_0, b)
    c1 = solve(q_1, a)

    f = q.subs({a:c0, b:c1})
    print("interpolator:")
    pprint(f)

def interp_quadratic():
    a, b, c, x = symbols('a b c x')
    q = a * x**2 + b * x + c

    pprint(q)

    q_dx = diff(q, x)

    """
    So here, what do we do? We can set boundary
    equation, ensuring 1st order continuity.
    """

    q_0 = q.subs('x', 0)
    q_1 = q.subs('x', 1)
    q_dx_0 = q_dx.subs('x', 0)
    q_dx_1 = q_dx.subs('x', 1)

    print("q")
    pprint(q_0)
    pprint(q_1)
    print("q_dx")
    pprint(q_dx_0)
    pprint(q_dx_1)

    """
    But after that?
    """

    c0 = solve(q_1 + q_dx_1, a)
    c1 = solve(q_1 + q_dx_0, b)
    c2 = solve(q_0 + q_dx_1, c)

    f = q.subs({a: c0, b: c1, c:c2})
    print("interpolator:")
    pprint(f)

def create_quartic(p0, p1, p2, p3, p4, x):
    
    a,b,c,d,e = symbols('a b c d e')
    q = a * x**4 + b * x**3 + c * x**2 + d * x + e

    pprint(q)

    q_dx = diff(q, x)
    q_dxx = diff(q_dx, x)

    # pprint(q_dx)

    q_0 = q.subs('x', 0)
    q_1 = q.subs('x', 1)
    q_dx_0 = q_dx.subs('x', 0)
    q_dx_1 = q_dx.subs('x', 1)
    q_dxx_0 = q_dxx.subs('x', 0)
    q_dxx_1 = q_dxx.subs('x', 1)

    pprint(q_0)
    pprint(q_1)
    pprint(q_dx_0)
    pprint(q_dx_1)
    pprint(q_dxx_0)
    pprint(q_dxx_1)

    return 

def quartic_interp():
    p0, p1, p2, p3, p4, x = symbols('p0 p1 p2 p3, p4, x')
    curve = create_quartic(p0, p1, p2, p3, p4, x)
    

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)
    
    # interp_linear()
    # interp_quadratic()

    # cubic_interp();

    # cubic_interp_neighbors();

    # bicubic_interp()
    # bicubic_interpolation_tangents()
    # bicubic_interpolation_laplacian()
    # bicubic_interpolation_curvature()
    # bicubic_interpolation_combined()

    # quartic_interp()

    t = symbols('t')
    n = Rational(6)
    f = Rational(10)
    s = n / (pi * 2 * f)
    
    w = sin(pi * 2 * t) * exp(-(t**2) / (2 * s**2))
    pprint(w)

    w_1 = w.subs(t, Rational(1,2))
    pprint(w_1)

    w_1_int = integrate(w_1, (t, -1, 1))
    pprint(w_1_int)


if __name__ == "__main__":
    main()
