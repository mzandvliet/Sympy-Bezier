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
    '''
    Like the below 1-d version, but now with full
    motion in the plane
    '''
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
    '''
    Todo:
    - Prove that for each control point you can chose a specific Qi,
    such that the calculated velocity for Pi represents the change
    with respect to the entire wave section.

    See Readability notes for details, but one way to try it is through
    defining the gradient as the integral of gradients along all Q, and
    seeing whether something neat pops out.

    The hunch would be that terms left and right of such a chosen Q
    would cancel each other out.

    - If that proof pans out, derive a formula that expresses such a Q_i
    for any P_i

    Factoid: End points P_1 and P_5 will have Q_1=0 and Q_5=1, so we
    conclude that at least for those points,

    Another hunch is that if a representative Q_i exists, it would probably
    be the point Q on the curve that is **most affected** by change in P_i
    '''

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

    dp_dp1 = diff(p[1], p1[1])
    dp_dp2 = diff(p[1], p2[1])
    dp_dp3 = diff(p[1], p3[1])
    dp_dp4 = diff(p[1], p4[1])
    dp_dp5 = diff(p[1], p5[1])
    dp1_dt = dp_dp1 * dp_dt
    dp2_dt = dp_dp2 * dp_dt
    dp3_dt = dp_dp3 * dp_dt
    dp4_dt = dp_dp4 * dp_dt
    dp5_dt = dp_dp5 * dp_dt

    q1, q2, q3, q4, q5 = symbols('q1 q2 q3 q4 q5')
    dp1_dt = dp1_dt.subs(q, q1)
    dp2_dt = dp2_dt.subs(q, q2)
    dp3_dt = dp3_dt.subs(q, q3)
    dp4_dt = dp4_dt.subs(q, q4)
    dp5_dt = dp5_dt.subs(q, q5)

    ddp1_ddt = diff(dp1_dt, t)
    ddp2_ddt = diff(dp2_dt, t)
    ddp3_ddt = diff(dp3_dt, t)
    ddp4_ddt = diff(dp4_dt, t)
    ddp5_ddt = diff(dp5_dt, t)

    common, exprs = cse((ddp1_ddt, ddp2_ddt, ddp3_ddt, ddp4_ddt, ddp5_ddt), numbered_symbols('a'))
    print_code(common, exprs)

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)
    # intersect_quadratic_with_line()
    quartic_bezier_wave_equation_1d()

if __name__ == "__main__":
    main()
