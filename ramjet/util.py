from sympy import *
from sympy.physics.vector import *
import matplotlib.pyplot as plt


def invert_dict(dict):
    return {v: k for k, v in dict.items()}

# convert sympy expression to inline latex ready for showing in MatPlotLib


def show_expr_latex(expr):
    show_latex_str(latex(expr, mode='inline'))

    # https://docs.sympy.org/latest/modules/printing.html#sympy.printing.latex.latex

# Render a latex string using MatPlotLib


def show_latex_str(latex_str):
    plt.rc('text', usetex=True)
    plt.title(latex_str)
    plt.show()

    # https://matplotlib.org/users/usetex.html
    # Note: requires working install of LateX on PATH


def to_latex_docstring(expr):
    doc_srt = r"""
\documentclass{article}
\begin{document}
    %s
\end{document}""" % (
        latex(expr, mode='inline')
    )

    return doc_srt

def print_pretty(common, exprs):
    print("\n---------------terms-----------------\n")

    for t in common:
        pprint(t)

    print("\n----------------roots-------------------\n")

    for expr in exprs:
        pprint(expr)


def print_code(common, exprs):
    print("\n----------------terms-------------------\n")

    for t in common:
        print(csharp(t))

    print("\n----------------roots-------------------\n")

    for i, expr in enumerate(exprs):
        print("float root_%d = " % i + replace_vector_vars(ccode(expr)) + ";")


'''
Poor man's CodeGen: take ccode output, and format it to valid
C# through very dodge string manipulation

Todo: make a proper sympy codegen module for C#, or pipe results
into Roslyn AST using F# or something...
'''


def csharp(code):
    code = ccode(code)
    code = code[1:-1]

    comma_idx = code.find(",")
    code = code[0:comma_idx] + " =" + code[comma_idx+1:]
    code = "float " + code + ";"

    code = format_math_funcs(code)
    code = replace_vector_vars(code)
    code = format_floats(code)

    return code

def format_math_funcs(code):
    code = code.replace("pow", "math.pow")
    code = code.replace("sqrt", "math.sqrt")
    return code

def replace_vector_vars(code):
    '''
    Scan through string finding occurances of 'identifier{n}_{x,y,z}'
    Replace each with 'curve[*].*'
    '''

    xyz = ['x', 'y', 'z']

    pos = 0
    while True:
        pos = code.find('_', pos)
        if pos == -1:
            break

        if code[pos+1] in xyz:
            idx = int(code[pos-1])
            idx -= 1

            name_start = find_identifier_backwards_from(code, pos-2)
            name = code[name_start:pos-1]

            code = code[0:name_start] + name + "[" + str(idx) + "]." + code[pos+1:]
        else:
            pos += 1
    return code

def find_identifier_backwards_from(code, end):
    for i in range(end, -1, -1):
        if not code[i].isalpha():
            return i+1
    return 0

def format_floats(code):
    '''
    Scan through string finding occurances of 'n.n'
    Replace each with 'n.nf'
    '''
    pos = 0
    while True:
        pos = code.find('.', pos)
        if pos == -1:
            break

        if code[pos-1].isdigit() and code[pos+1].isdigit() and code[pos+2] != 'f':
            code = code[0:pos+2] + "f" + code[pos+2:]
        else:
            pos += 1
    return code
