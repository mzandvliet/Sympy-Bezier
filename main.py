from sympy import *
import matplotlib.pyplot as plt
import math

# convert sympy expression to inline latex ready for showing in MatPlotLib
def show_expr_latex(expr):
    show_latex(latex(expr, mode='inline'))

    # https://docs.sympy.org/latest/modules/printing.html#sympy.printing.latex.latex

# Render some latex using MatPlotLib
def show_latex(latex_str):
    plt.rc('text', usetex=True)
    plt.title(latex_str)
    plt.show()

    # https://matplotlib.org/users/usetex.html
    # Note: requires working install of LateX on PATH


def literal_sqrt():
    result = math.sqrt(9)
    print("Result: " + str(result))

def symbolic_sqrt():
    result = sqrt(8)
    print("Result: " + str(result)) # auto-simplifies to 2*sqrt(2)

def simple_expr():
    x, y = symbols('x y')
    expr = x + 2*y
    pprint(expr)

    expr2 = expr + 1
    pprint(expr2)

    expr3 = x * expr  # doesn't expand this automatically
    pprint(expr3) # unicode pretty print

    show_expr_latex(expr3) # Show latex



def main():
    # literal_sqrt()
    # symbolic_sqrt()
    simple_expr()
    
if __name__ == "__main__":
    main()
