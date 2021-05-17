import math
from sympy import *

def wave_equation():
    x, t, c = symbols('x t c')

    u = Function('u')
    
    wave_eq = Eq(diff(u(t,x), t, t), c**2 * diff(u(t,x), x, 2))
    pprint(wave_eq)

def wave_equation_2():
    x, t, c = symbols('x t c')

    u = Function('u')(x, t) # object representing function u evaluated at x and t

    wave_eq = Eq(diff(u, t, 2), c**2 * diff(u, x, 2))
    pprint(wave_eq)

def matrices():
    m = Matrix([
        [1, 2], # row at top
        [3, 4]  # row at bottom
    ])

    pprint(m)

    m = Matrix([1,2,3]) # a column vector

    pprint(m)

    x, y, z = symbols('x y z')

    mat_a = Matrix([
        [1, 0, 1],
        [-1, 2, 3],
        [1, 2, 3]
    ])

    mat_b = Matrix([x, y, z])

    mat_c = mat_a * mat_b
    pprint(mat_c)

    pprint(mat_c.jacobian([x,y,z]))


def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)
    
    # wave_equation_2()
    matrices()


if __name__ == "__main__":
    main()
