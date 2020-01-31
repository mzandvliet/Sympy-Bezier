import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *

def join_meet(a,b):
    return a + b

def intersect_line_segments_projetive():
    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')

    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')



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

    --

    Boundary conditions:

    Since the endpoints of the string do not move, reason about where the
    energy goes, and whether we are currently even accounting for it.
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


def quadratic_piecewize_bezier_wave_equation_1d():
    '''
    Could we use quadratics? They're the simplest.
    But if you take their second derivative, it
    is a constant. Maybe that's good? The constant
    can still vary between individual sections.

    This might be easier with a B-spline formulation
    '''
    # q = curve parameter
    # t = time
    # c = wave equation constant
    q_a, q_b, t, c = symbols('q_a q_b t c')

    # Quartic curve and its second derivative
    p1 = symbolic_vector_2d('p_1')
    p2 = symbolic_vector_2d('p_2')
    p3 = symbolic_vector_2d('p_3')
    p4 = symbolic_vector_2d('p_4')
    p5 = symbolic_vector_2d('p_5')
    bases = bezier_bases(2, q_a)
    bases_dd = bezier_bases(0, q_b)

    points_a = (p1, p2, p3)
    points_b = (p3, p4, p5)
    points_a_d = differentiate_curve_points(points_a)
    points_a_dd = differentiate_curve_points(points_a_d)
    points_b_d = differentiate_curve_points(points_b)
    points_b_dd = differentiate_curve_points(points_b_d)

    p_a = make_bezier(points_a, bases)(q_a)
    p_a_dd = make_bezier(points_a_dd, bases_dd)(q_a)
    p_b = make_bezier(points_b, bases)(q_b)
    p_b_dd = make_bezier(points_b_dd, bases_dd)(q_b)

    # Since laplacian is the second derivative of this curve
    # our wave equation turns out to be this:
    ddp_a_ddt = c**2 * p_a_dd[1]
    dp_a_dt = integrate(ddp_a_ddt, t)
    ddp_b_ddt = c**2 * p_b_dd[1]
    dp_b_dt = integrate(ddp_b_ddt, t)

    # Now express that in terms of motion of the control points
    # using the chain rule:
    # dp3_dt = dp/dp3 * dp_dt

    '''
    Hmm, at this point, need to carefully consider how
    to join the two pieces together. What are the
    equations of motion for the middle bit?
    Does control point B still experience weak influence
    from the neighboring piece? Yes. How do we represent it?
    '''

    dp_dp1 = diff(p_a[1], p1[1])
    dp_dp2 = diff(p_a[1], p2[1])
    dp_dp3 = diff(p_a[1], p3[1]) + diff(p_b[1], p3[1])
    dp_dp4 = diff(p_b[1], p4[1])
    dp_dp5 = diff(p_b[1], p5[1])
    dp1_dt = dp_dp1 * dp_a_dt
    dp2_dt = dp_dp2 * dp_a_dt
    dp3_dt = dp_dp3 * dp_a_dt
    dp4_dt = dp_dp4 * dp_b_dt
    dp5_dt = dp_dp5 * dp_b_dt

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

    common, exprs = cse((ddp1_ddt, ddp2_ddt, ddp3_ddt,
                         ddp4_ddt, ddp5_ddt), numbered_symbols('a'))
    print_code(common, exprs)

'''
Implement de Carpentier's approximate linear drag ballistics, and
generate Unity code for it. Make a cool action game with it.
'''
def magnitude(vec):
    return sqrt(vec.dot(vec))

def exp_rat(x):
    '''
    First order rational approximation of e^x around x=0
    '''
    return (2+x)/(2-x)

def linear_drag():
    mag = Function('mag')

    t = symbols('t')
    g = symbols('g') #symbolic_vector_2d('g')
    v_medium = symbols('v-medium')  # symbolic_vector_2d('v-medium')

    p_0 = symbols('p0') #symbolic_vector_2d("p0")
    v_0 = symbols('v0') #symbolic_vector_2d("v0")
    v_terminal = symbols('v-terminal')  # symbolic_vector_2d('v-terminal')
    
    v_inf = symbols('v-inf') #v_terminal + v_medium
    k = symbols('k') #Rational(1, 2) * (g / mag(v_terminal))
    
    p_linear = (1/(2*k)) * (v_0 - v_inf)*(1 - exp(-2*k*t)) + v_inf * t + p_0

    # remove p_0 so it doesn't get caught up in the following step
    # todo: why does sympy work p0 into the resulting ratio?
    p_linear = p_linear - p_0
    # Replace exp by the rational approximation
    p_linear = p_linear.replace(exp, exp_rat).replace(mag, magnitude)
    p_linear = simplify(p_linear)
    # add p_0 back
    p_linear = p_linear + p_0
    # we have now derived the basic formula from the paper:
    # p_linear = (((v_0 + k*t*v_inf) * t) / (1 + k*t)) + p0

    ppprint(p_linear)


def ppprint(expr):
    '''Pretty printing function that auomatically passes formatting config, because I cannot
    get init_printing() to work'''
    pprint(expr, use_unicode=True, num_columns=240)

def main():
    # intersect_quadratic_with_line()
    # quartic_bezier_wave_equation_1d()
    linear_drag()

if __name__ == "__main__":
    main()
