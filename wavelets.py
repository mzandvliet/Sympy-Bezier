import math
from sympy import *
from ramjet.math import *
from ramjet.util import *

def morlet_orthogonality():
    # pprint(integrate(sin(t) * sin(t+pi/2), (t, -1, 1)).evalf())

    x, t = symbols('x t')
    n = Rational(6)
    f = Rational(10)
    s = n / (pi * 2 * f)

    w = cos(pi * 2 * f * (x+t)) * exp(-(t**2) / (2 * s**2))
    # pprint(w)

    # w_1 = w.subs([(t, Rational(1, 2)), (x, Rational(0))])
    # pprint(w_1)

    # w_1_int = integrate(w.subs(x, Rational(0)), (t, -1, 1))
    # pprint(w_1_int.evalf())

    w_a = w.subs([(x, Rational(0))])
    w_b = w.subs([(x, Rational(1, 2))])
    pprint(integrate(w_a * w_b, (t, -1, 1)).evalf())

def sum_test():
    x = symbols('x')
    test = lambda t, y: t**y
    my_sum = Sum(test(x,2), (x, -16, +16))
    pprint(my_sum)
    pprint(my_sum.doit())

def morlet(f,x,n,t):
    s = n / (pi * 2 * f)
    return cos(pi * 2 * f * (x+t)) * exp(-(t**2) / (2 * s**2))

def discrete_convolution_derivatives():
    # Goal: find gradient of convolution s.t. we can maximize dot product bewteen signal and wavelet

    x, t, wt = symbols('x t wt')
    n = 6

    f_a = 10
    f_b = 9.9

    w_a = morlet(f_a, 0, n, t)
    w_b = morlet(f_b, 0, n, t)

    convolution = Sum(w_a * w_b, (t, -16, +16))
    pprint(convolution)
    pprint(convolution.evalf())


def freq_modulation_2():
    q = symbols('q')

    p0 = Matrix([0, 1])
    p1 = Matrix([1, 2])
    p2 = Matrix([2, 2])

    # Todo: use a bezier curve here, so we can fit many points in time-frequency space?
    # CAREFUL: don't mix up linear time and the spline parameter 't'

    # options:
    # integrate and maximize numerically, which is the easiest

    points = (p0, p1, p2)
    bases = bezier_bases(2, q)
    p = make_bezier(points, bases)(q)

    # osc = sin(f * t * pi * 2)
    # period = 1 / f

    f = p[1]
    t = p[0]

    phase = f * t

    pprint(phase)

    phaseStep = solve(phase-1, q)

    pprint(phaseStep)

def freq_modulation():
    t = symbols('t')

    t0 = Rational(0)
    f0 = Rational(1)

    t1 = Rational(1)
    f1 = Rational(2)

    # Todo: use a bezier curve here, so we can fit many points in time-frequency space?
    # CAREFUL: don't mix up linear time and the spline parameter 't'

    # options:
    # simplify by inverting 1/f and expressing points in time-wavelength space?
    # integrate and maximize numerically, which is the easiest

    f = f0 * (1-t) + f1 * t

    # osc = sin(f * t * pi * 2)
    # period = 1 / f

    phase = f * t

    pprint(phase)

    solution = solve(phase-1, t)

    pprint(solution)

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)

    # sum_test()

    # morlet_orthogonality()
    discrete_convolution_derivatives()

if __name__ == "__main__":
    main()
