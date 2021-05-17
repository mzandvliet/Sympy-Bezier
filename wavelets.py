import math
from sympy import *

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

    # morlet_orthogonality()
    freq_modulation()

if __name__ == "__main__":
    main()
