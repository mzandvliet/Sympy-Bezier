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

    print("\n--------------solutions----------------\n")

    for expr in exprs:
        pprint(expr)


def print_code(common, exprs):
    print("\n/*----------------terms-------------------*/\n")

    for t in common:
        print(csharp(t))

    print("\n/*--------------solutions------------------*/\n")

    for i, expr in enumerate(exprs):
        expr = ccode(expr);
        expr = format_math_funcs(expr)
        expr = replace_matrix_vars(expr)
        expr = format_floats(expr)
        print("float output_%d = " % i + expr + ";")


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
    # code = replace_vector_vars(code)
    code = replace_matrix_vars(code)
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

    # tri_quad_indices = {
    #     "002" : 0,
    #     "011" : 1,
    #     "020" : 2,
    #     "101" : 3,
    #     "110" : 4,
    #     "200" : 5,
    # }

    tri_quart_indices = {
        "004": 0,
        "013": 1,
        "022": 2,
        "031": 3,
        "040": 4,
        "103": 5,
        "112": 6,
        "121": 7,
        "130": 8,
        "202": 9,
        "211": 10,
        "220": 11,
        "301": 12,
        "310": 13,
        "400": 14,
    }

    bases = ['x', 'y', 'z', 'w']

    pos = 0
    while True:
        pos = code.find('_', pos)
        if pos <= 1:
            break

        # in bases and code[pos-1].isdigit()
        if code[pos-1] == 'p' and code[pos+1].isdigit():
            # Quad point indexing, e.g. p_0.x -> p[0].x

            numDigits = 1
            while (code[pos+numDigits+1].isdigit()):
                numDigits+=1

            index = int(code[pos+1:pos+1+numDigits]) - 1

            code = code[0:pos-1] + "p[" + str(index) + "]." + code[pos+1+numDigits+1:]
        
        # elif code[pos-1] == 't':
        #     # Tri point indexing, e.g. t_101.x -> p[3].x
        #     # Todo: find out which degree of patch we're dealing with

        #     idxString = code[pos+2:pos+5]
        #     idx = tri_quart_indices[idxString]

        #     code = code[pos] + "p[" + idx + "]." + code[pos+5:]

        else:
            pos += 1
    return code

def find_identifier_backwards_from(code, end):
    for i in range(end, -1, -1):
        if code[i].isalpha():
            return i
    return 0


def replace_matrix_vars(code):
    pos = 0
    while True:
        pos = code.find('p', pos) # returns -1 if none found

        if pos == -1:
            break

        if code[pos+2] == '_' and code[pos+1].isdigit() and code[pos+3].isdigit():
            x = int(code[pos+1])
            y = int(code[pos+3])

            code = code[0:pos] + "p[%i][%i]" % (x, y) + code[pos+4:]
            
        pos += 1

    return code

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
