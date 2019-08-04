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

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)
    intersect_quadratic_with_line()

if __name__ == "__main__":
    main()
