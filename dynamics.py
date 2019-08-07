import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *


def intersect_quadratic_with_line():
    t = symbols('t')
    p1 = symbolic_vector_2d('p1')
    v1 = symbolic_vector_2d('v1')
    g = symbolic_vector_2d('g')

    s = symbols('s')
    l1 = symbolic_vector_2d('l1')
    l2 = symbolic_vector_2d('l2')

    trajectory = p1 + v1 * t + Rational(1,2) * g * t**2
    line = l1 * (1-s) + l2 * s

    solution = solve(trajectory[1], t)
    pprint(solution)

    common, exprs = cse(solution, numbered_symbols('a'))
    print_code(common, exprs)

def quartic_bezier_wave_equation():
    # q = curve parameter
    # t = time
    # c = wave equation constant
    q, t, c = symbols('q t c')

    # Quartic curve and its second derivative
    p1 = symbolic_vector_2d('p_1')
    p2 = symbolic_vector_2d('p_2')
    p3 = symbolic_vector_2d('p_3')
    p4 = symbolic_vector_2d('p_4')
    p5 = symbolic_vector_2d('p_5')
    bases = bezier_bases(4, q)
    bases_dd = bezier_bases(2, q)

    points = (p1, p2, p3, p4, p5)
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    p = make_bezier(points, bases)(q)
    pdd = make_bezier(points_dd, bases_dd)(q)

    # Since laplacian is the second derivative of this curve
    # our wave equation turns out to be this:
    ddp_ddt = c**2 * pdd
    dp_dt = integrate(ddp_ddt, t)
    # pprint(dp_dt)

    # Now express that in terms of motion of the control points
    # using the chain rule:
    # dp3_dt = dp/dp3 * dp_dt

    dp_dp3 = diff(p, p3)
    dp3_dt = tensorproduct(dp_dp3, dp_dt)
    # pprint(dp_dp3)
    pprint(dp3_dt)

def quartic_bezier_wave_equation_1d():
    # q = curve parameter
    # t = time
    # c = wave equation constant
    q, t, c = symbols('q t c')

    # Quartic curve and its second derivative
    p1 = symbolic_vector_2d('p_1')
    p2 = symbolic_vector_2d('p_2')
    p3 = symbolic_vector_2d('p_3')
    p4 = symbolic_vector_2d('p_4')
    p5 = symbolic_vector_2d('p_5')
    bases = bezier_bases(4, q)
    bases_dd = bezier_bases(2, q)

    points = (p1, p2, p3, p4, p5)
    points_d = differentiate_curve_points(points)
    points_dd = differentiate_curve_points(points_d)

    p = make_bezier(points, bases)(q)
    pdd = make_bezier(points_dd, bases_dd)(q)

    # Since laplacian is the second derivative of this curve
    # our wave equation turns out to be this:
    ddp_ddt = c**2 * pdd[1]
    dp_dt = integrate(ddp_ddt, t)
    # pprint(dp_dt)

    # Now express that in terms of motion of the control points
    # using the chain rule:
    # dp3_dt = dp/dp3 * dp_dt

    dp_dp3 = diff(p[1], p3[1])
    dp3_dt = dp_dp3 * dp_dt
    # pprint(dp3_dt)

    common, exprs = cse(dp3_dt, numbered_symbols('a'))
    print_code(common, exprs)

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)
    # intersect_quadratic_with_line()
    quartic_bezier_wave_equation_1d()

if __name__ == "__main__":
    main()
