import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *

'''
Todo:

Vlachos PN Triangles paper

I don't understand a specific thing:
We're generating, I think, in the exact same way that PN triangle paper suggests,
yet for [i,j,k] = [3,0,0] the get term:

p300*w^3

whereas I get:
p300*u^3

Because u^i...

Why? Might be a printing mistake or something.

With ijk cycled backwards one tick, I get:
(symbols[0]**j * symbols[1]**k * symbols[2]**i)
The results of which match the paper. Meh.
'''

# de Casteljau recursive evaluation

def barycentric_triangle(a, b, c, params):
    p = a * params[2] + b * params[1] + c * params[0]
    return p

def barycentric_triangle_2(points, params):
    t0_0 = barycentric_triangle(points[0], points[1], points[3], params)
    t0_1 = barycentric_triangle(points[1], points[2], points[4], params)
    t0_2 = barycentric_triangle(points[3], points[4], points[5], params)

    t1_0 = barycentric_triangle(t0_0, t0_1, t0_2, params)

    return t1_0


def triangular_patch_symbolic(points, params, degree):
    patch = 0

    s = 0
    for i in range(0, degree+1):
        for j in range(0, degree+1):
            for k in range(0, degree+1):
                if (i+j+k != degree):
                    continue

                tri = trinomial(degree, i, j, k)
                p_ijk = points[s] * tri * (params[0]**i * params[1]**j * params[2]**k)
                patch += p_ijk
                s += 1

    return patch

def triangular_patch(params, degree, basis):
    patch = Matrix([0]*len(basis))

    for i in range(0, degree+1):
        for j in range(0, degree+1):
            for k in range(0, degree+1):
                if (i+j+k != degree):
                    continue

                tri = trinomial(degree, i, j, k)
                p_ijk = symbolic_vector('p%i%i%i' % (
                    i, j, k), basis) * tri * (params[0]**i * params[1]**j * params[2]**k)
                print('%i%i%i' % (i, j, k))
                patch += p_ijk

    return patch

# polynomial evaluation

def triangular_patch_with_points(points, symbols, degree):
    dimensions = points[0].shape[0]
    patch = Matrix([0]*dimensions)

    pointIdx = 0

    for i in range(0, degree+1):
        for j in range(0, degree+1):
            for k in range(0, degree+1):
                if (i+j+k != degree):
                    continue

                tri = trinomial(degree, i, j, k)
                p_ijk = points[pointIdx] * tri * \
                    (symbols[0]**i * symbols[1]**j * symbols[2]**k)
                pointIdx += 1
                patch += p_ijk

    return patch


def triangular_patch_3d_du(symbols, degree):
    patch = Matrix([0, 0, 0])

    for i in range(0, degree):
        for j in range(0, degree):
            for k in range(0, degree):
                if (i+j+k != degree-1):
                    continue

                tri = trinomial(degree, i, j, k)
                p_ijk = symbolic_vector_3d('p%i%i%i' % (
                    i+1, j, k)) * tri * (symbols[0]**i * symbols[1]**j * symbols[2]**k)
                patch += p_ijk

    return patch


def triangular_patch_3d_dv(symbols, degree):
    patch = Matrix([0, 0, 0])

    for i in range(0, degree):
        for j in range(0, degree):
            for k in range(0, degree):
                if (i+j+k != degree-1):
                    continue

                tri = trinomial(degree, i, j, k)
                p_ijk = symbolic_vector_3d('p%i%i%i' % (
                    i, j+1, k)) * tri * (symbols[0]**i * symbols[1]**j * symbols[2]**k)
                patch += p_ijk

    return patch


def quadratic_triangular_patch_3d():
    u, v, w = symbols('u v w')
    # w = 1 - u - v

    patch = triangular_patch((u, v, w), 2, BASIS_3D)

    common, exprs = cse(patch, numbered_symbols('a'))
    print_code(common, exprs)


def quadratic_triangular_patch_3d_prove_derivatives():
    u, v, w = symbols('u v w')
    # w = 1 - u - v

    patch = triangular_patch((u, v, w), 2, BASIS_3D)
    pprint(diff(patch, u) - triangular_patch_3d_du((u, v, w), 2))


def quadratic_triangular_patch_3d_silhouette():
    u, v, w = symbols('u v w')
    # v = 1-u
    # w = 0
    '''
    Trying to get closed form solution for silhouette along
    the u-v edge. But pluggin the above into differentiation
    goes wrong of course.

    with prepared du & dv patches for tangent calculations,
    this would work.
    '''

    patch = triangular_patch((u, v, w), 2, BASIS_3D)
    patch_du = triangular_patch_3d_du((u, v, w), 2)
    patch_dv = triangular_patch_3d_dv((u, v, w), 2)

    patch_normal = patch_du.cross(patch_dv)
    pprint(patch_du)

    # viewpoint = Matrix([0,0,0])
    # viewdir = patch - viewpoint

    # silhouette = viewdir.dot(patch_normal)**2
    # silhouette = diff(silhouette, u)
    # silhouette = solveset(silhouette, u, domain=S.Reals)
    # print(silhouette)

    # solutions = list(map(lambda s: s.subs(v, 1-u).subs(w, 0), silhouette))

    # silhouette = silhouette.subs(v, 1-u).subs(w, 0)
    # silhouette = simplify(silhouette)
    # print(solutions[0])

    # common, exprs = cse(silhouette, numbered_symbols('a'))
    # print_code(common, exprs)


def quadratic_triangular_patch_3d_silhouette_gradient():
    u, v, w = symbols('u v w')
    # w = 1 - u - v

    patch = triangular_patch((u, v, w), 2, BASIS_3D)
    patch_du = diff(patch, u)
    patch_dv = diff(patch, v)

    patch_normal = patch_du.cross(patch_dv)

    viewpoint = Matrix([0, 0, 0])
    viewdir = patch - viewpoint

    silhouette = viewdir.dot(patch_normal)**2

    grad_u = diff(silhouette, u)
    grad_v = diff(silhouette, v)

    common, exprs = cse((grad_u, grad_v), numbered_symbols('a'))
    print_code(common, exprs)


def quadratic_rational_triangular_patch_3d():
    u, v, w = symbols('u v w')

    patch = triangular_patch((u, v, w), 2, BASIS_4D)
    pprint(diff(patch, u), use_unicode=True, num_columns=140)


def quartic_rational_triangular_patch_3d():
    '''
    Todo: Find biquadratic forms
    '''
    u, v, w = symbols('u v w')

    patch = triangular_patch((u, v, w), 4, BASIS_4D)

    common, exprs = cse(patch, numbered_symbols('a'))
    print_code(common, exprs)


def quadratic_rational_triangular_patch_3d_embedded_line():
    '''
    Full euclidean 3d evaluation of embedded linear
    segment. Had presumed this to be a geodesic
    on the sphere in case of rational octant, but
    the general result is a quartic.
    '''
    u, v, w, t = symbols('u v w t')

    patch = triangular_patch((u, v, w), 2, BASIS_4D)

    uvw1 = symbolic_vector_3d('uvw1')
    uvw2 = symbolic_vector_3d('uvw2')
    # uvw1[2] = 1 - uvw1[0] - uvw1[1]
    # uvw2[2] = 1 - uvw2[0] - uvw2[1]
    bases_uv = bezier_bases(1, t)
    uvw = make_bezier((uvw1, uvw2), bases_uv)(t)

    p = patch.subs(u, uvw[0]).subs(v, uvw[1]).subs(w, uvw[2])
    poly = to_polynomial(p[0], t)
    print("Got polynomial of degree: %i" % (poly.degree()))
    # pprint(expand(p[0]))

    common, exprs = cse(p, numbered_symbols('a'))
    print_code(common, exprs)


def quadratic_rational_sphere_octant_3d_embedded_line():
    u, v, w, t = symbols('u v w t')

    #halfsqrt2 = Rational(1,2) * 2**Rational(1,2)
    halfsqrt2 = Rational(1, 1)
    patch = triangular_patch_with_points([
        Matrix([0, 0, 1, 1]),  # 002
        Matrix([0, 1, 1, 1]) * halfsqrt2,  # 011
        Matrix([0, 1, 0, 1]),  # 020
        Matrix([1, 0, 1, 1]) * halfsqrt2,  # 101
        Matrix([1, 1, 0, 1]) * halfsqrt2,  # 110
        Matrix([1, 0, 0, 1]),  # 200
    ], (u, v, w), 2)

    uvw1 = symbolic_vector_3d('uvw1')
    uvw2 = symbolic_vector_3d('uvw2')
    # uvw1[2] = 1 - uvw1[0] - uvw1[1]
    # uvw2[2] = 1 - uvw2[0] - uvw2[1]
    bases_uv = bezier_bases(1, t)
    uvw = make_bezier((uvw1, uvw2), bases_uv)(t)

    p = patch.subs(u, uvw[0]).subs(v, uvw[1]).subs(w, uvw[2])
    poly = to_polynomial(p[0], t)
    print("Got polynomial of degree: %i" % (poly.degree()))

    pprint(p[0])

    # common, exprs = cse(p, numbered_symbols('a'))
    # print_code(common, exprs)


def quadratic_rational_sphere_octant_3d_norm_proof():
    '''
    Here we check the magnitude of the vector we
    get in the center. All vectors we get should
    have unit norm, of course. But!

    0.37037037037037⋅√2 + 0.666666666666667 == ~1.19

    So evidently, we do not actually have a sphere here.
    '''
    u, v, w = symbols('u v w')

    u = Rational(1, 3)
    v = Rational(1, 3)
    w = 1 - u - v

    halfsqrt2 = Rational(1, 2) * 2**Rational(1, 2)
    patch = triangular_patch_with_points([
        Matrix([0, 0, 1, 1]),  # 002
        Matrix([0, 1, 1, 1]) * halfsqrt2,  # 011
        Matrix([0, 1, 0, 1]),  # 020
        Matrix([1, 0, 1, 1]) * halfsqrt2,  # 101
        Matrix([1, 1, 0, 1]) * halfsqrt2,  # 110
        Matrix([1, 0, 0, 1]),  # 200
    ], (u, v, w), 2)

    p3d = Matrix([patch[0] / patch[3], patch[1] /
                  patch[3], patch[2] / patch[3]])
    norm = sqrt(p3d.dot(p3d))

    pprint(norm.evalf())


def quadratic_rational_sphere_octant_3d_norm_gradient():
    u, v, w = symbols('u v w')

    diagonalWeight = symbols('diag_w')
    patch = triangular_patch_with_points([
        Matrix([0, 0, 1, 1]),  # 002
        Matrix([0, 1, 1, 1]) * diagonalWeight,  # 011
        Matrix([0, 1, 0, 1]),  # 020
        Matrix([1, 0, 1, 1]) * diagonalWeight,  # 101
        Matrix([1, 1, 0, 1]) * diagonalWeight,  # 110
        Matrix([1, 0, 0, 1]),  # 200
    ], (u, v, w), 2)

    p3d = Matrix([patch[0] / patch[3], patch[1] /
                  patch[3], patch[2] / patch[3]])
    norm = sqrt(p3d.dot(p3d))
    normError = (Rational(1, 1) - norm)**2

    dErrordWeight = diff(normError, diagonalWeight)

    common, exprs = cse(dErrordWeight, numbered_symbols('a'))
    print_code(common, exprs)


def quadratic_rational_triangular_patch_3d_geodesic_gradient():
    '''
    Starting from a random curve, give dControls/GeodesicError,
    such that we can numerically find geodesics
    '''
    u, v, w, t = symbols('u v w t')

    patch = triangular_patch((u, v, w), 2, BASIS_4D)

    uvw1 = symbolic_vector_3d('uvw1')
    uvw2 = symbolic_vector_3d('uvw2')
    uvw3 = symbolic_vector_3d('uvw3')

    bases_uv = bezier_bases(2, t)
    uvw = make_bezier((uvw1, uvw2, uvw3), bases_uv)(t)

    p = patch.subs(u, uvw[0]).subs(v, uvw[1]).subs(w, uvw[2])

    # Piecewise Euclidean quadrance, stepping along curve
    # quadrance_steps = 8
    # quadrance = Rational(0,1)
    # t_step = Rational(1, quadrance_steps)
    # for i in range(0, quadrance_steps-1):
    #     t_i = t_step * i
    #     p_delta = p.subs(t, t_i+t_step) - p.subs(t, t_i)
    #     quadrance += p_delta.dot(p_delta)

    # Quadrance based on control point chord distance
    uvw_delta1 = uvw2 - uvw1
    uvw_delta2 = uvw3 - uvw2
    quadrance = uvw_delta1.dot(uvw_delta1) + uvw_delta2.dot(uvw_delta2)

    quadranceGradU = diff(quadrance, uvw2[0])
    quadranceGradV = diff(quadrance, uvw2[1])
    quadranceGradW = diff(quadrance, uvw2[2])

    common, exprs = cse((quadranceGradU, quadranceGradV,
                         quadranceGradW), numbered_symbols('a'))
    print_code(common, exprs)


def cubic_triangular_patch_3d():
    u, v, w = symbols('u v w')

    patch = triangular_patch((u, v, w), 3, BASIS_3D)
    # pprint(patch, use_unicode=True, num_columns=140)
    pprint(diff(patch, u), use_unicode=True, num_columns=140)


def cubic_triangular_patch_3d_silhouette_gradient():
    u, v, w = symbols('u v w')
    # w = 1 - u - v

    patch = triangular_patch((u, v, w), 3, BASIS_3D)
    patch_du = diff(patch, u)
    patch_dv = diff(patch, v)

    patch_normal = patch_du.cross(patch_dv)

    viewpoint = Matrix([0, 0, 0])
    viewdir = patch - viewpoint

    silhouette = viewdir.dot(patch_normal)**2

    grad_u = diff(silhouette, u)
    grad_v = diff(silhouette, v)

    common, exprs = cse((grad_u, grad_v), numbered_symbols('a'))
    print_code(common, exprs)


def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)
    
    points = symbols("a b c d e f")  # left-to-right, bottom-to-top
    params = symbols("u v w")

    print("Polynomial construction")
    poly = triangular_patch_symbolic(points, params, 2)
    pprint(poly)

    print("Casteljau construction")
    cast = barycentric_triangle_2(points, params)
    pprint(cast)

    diff = expand(poly) - expand(cast)
    print("Difference: " + str(diff))


if __name__ == "__main__":
    main()
