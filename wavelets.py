from sympy import *
from ramjet.math import *
from ramjet.util import *


def LinearWavelengthStep(f0, a):
    df = f0 * (-1 + pow(2, sin(a)))
    return Matrix([cos(a) / (f0 + df), df])


def TimeFreqStep(t, f):
    df = 1 * (-1 + pow(2, f))
    return Matrix([t, df])

def IntegrateWavelengthStep():
    t0, f0, a = symbols('t0, f0, a')
    step = LinearWavelengthStep(f0, a)

    print("Linear discrete wave step: ")
    pprint(step)
    
    q = symbols('q')
    # p1 = symbolic_vector_2d('p1')
    # p2 = symbolic_vector_2d('p2')
    # p1 = Matrix([t0, f0])
    # p2 = p1 + step
    p1 = Matrix([t0, f0])
    p2 = Matrix([t0+1, f0+1])
    points = (p1, p2)
    bases = bezier_bases(1, q)
    p = make_bezier(points, bases)(q)

    dpdt = diff(p, q)
    tfdpdt = TimeFreqStep(dpdt[0], dpdt[1])

    integratedStep = integrate(tfdpdt, (q, 0, 1))

    print("Linear wave step integrated along 1st order bezier curve (should be equivalent): ")
    pprint(integratedStep)

def cubic_point_fit_gradient_time_frequency():
    t = symbols('t')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')
    p4 = symbolic_vector_2d('p4')
    q = symbolic_vector_2d('q')
    points = (p1, p2, p3, p4)
    bases = bezier_bases(3, t)
    p = make_bezier(points, bases)(t)

    delta = q - p
    error = delta.dot(delta)**2
    grad_p2 = diff(error, p2)
    grad_p3 = diff(error, p3)

    common, exprs=cse((grad_p2, grad_p3), numbered_symbols('a'))
    print(ccode(common))
    print(ccode(exprs))

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

def discrete_convolution_derivatives():
    # Goal: find gradient of convolution s.t. we can maximize dot product bewteen signal and wavelet

    x, t, wt=symbols('x t wt')
    n=6

    f_a=10
    f_b=9.9

    w_a=morlet(f_a, 0, n, t)
    w_b=morlet(f_b, 0, n, t)

    convolution=Sum(w_a * w_b, (t, -16, +16))
    pprint(convolution)
    pprint(convolution.evalf())

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)

    # sum_test()

    # morlet_orthogonality()
    # discrete_convolution_derivatives()
    # cubic_point_fit_gradient_time_frequency()
    IntegrateWavelengthStep()

if __name__ == "__main__":
    main()
