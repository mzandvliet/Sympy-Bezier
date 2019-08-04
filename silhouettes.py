import math
from sympy import *
from sympy.physics.vector import *
from functools import reduce
from ramjet.math import *
from ramjet.util import *

def set_to_zero(symbs, expr):
    subs = { k: v for k, v in map(lambda s: [s, 0], symbs) }
    return expr.subs(subs)

def to_oriented_cubic_curve_3d_xyz(expr):
    x1, y1, z1, y4, z4 = symbols('x1, y1, z1, y4, z4')
    expr_subbed = expr \
        .subs(x1, 0) \
        .subs(y1, 0) \
        .subs(z1, 0) \
        .subs(y4, 0) \
        .subs(z4, 0) \

    return (expr_subbed)

def silhouette_cubic_2d():
    t = symbols('t')

    vx = 0
    vy = 0

    symbs = symbols('x1 x2 x3 x4 y1 y2 y3 y4')
    x1, x2, x3, x4, y1, y2, y3, y4 = symbs

    xd1 = 3 * (x2 - x1)
    xd2 = 3 * (x3 - x2)
    xd3 = 3 * (x4 - x3)
   
    yd1 = 3 * (y2 - y1)
    yd2 = 3 * (y3 - y2)
    yd3 = 3 * (y4 - y3)

    bases = bezier_bases(3, t)
    bases_d = bezier_bases(2, t)

    points_x = (x1, x2, x3, x4)
    points_y = (y1, y2, y3, y4)
    points_x_d = (xd1, xd2, xd3)
    points_y_d = (yd1, yd2, yd3)

    x = make_bezier(points_x, bases)(t)
    y = make_bezier(points_y, bases)(t)
    xd = make_bezier(points_x_d, bases_d)(t)
    yd = make_bezier(points_y_d, bases_d)(t)

    normal_x = -yd
    normal_y = xd

    viewdir_x = x - vx
    viewdir_y = y - vy

    solution = viewdir_x * normal_x + viewdir_y * normal_y
    solution = expand(solution)
    solution = to_oriented_cubic_curve_3d_xyz(solution)

    poly = to_polynomial(solution, t)
    print("Got polynomial of degree: " + str(poly.degree()))

    pprint(solution)

    # solution = solveset(solution, t)
    # common, exprs = cse(solution, numbered_symbols('a'))

    # print_code(common, exprs)

def silhouette_quadratic_2d():
    # The setup:

    # a curve (const), we need to express a point at t, and its normal
    # a view point (const)
    # direction from view point to curve point
    # dot product of normal and view direction
    # find t where that dot product = 0

    t = symbols('t')

    view_point = symbolic_vector_2d('viewPoint')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')

    p1 = Matrix([0, 0])

    pd1 = 2 * (p2 - p1)
    pd2 = 2 * (p3 - p2)

    bases = bezier_bases(2, t)
    bases_d = bezier_bases(1, t)

    points = (p1, p2, p3)
    points_d = (pd1, pd2)

    p = make_bezier(points, bases)(t)
    pd = make_bezier(points_d, bases_d)(t)

    normal = Matrix([-pd[1], pd[0]])
    viewdir = p - view_point

    solution = viewdir.dot(normal)
    solution = expand(solution)

    poly = to_polynomial(solution, t)
    print("Got polynomial of degree: " + str(poly.degree()))

    solution = solveset(solution, t)
    common, exprs = cse(solution, numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_2d_gradient():
    t = symbols('t')

    view_point = symbolic_vector_2d('viewPoint')

    p1 = symbolic_vector_2d('p1')
    p2 = symbolic_vector_2d('p2')
    p3 = symbolic_vector_2d('p3')

    # pd1 = 2 * (p2 - p1)
    # pd2 = 2 * (p3 - p2)
    n1 = symbolic_vector_2d('n1')
    n2 = symbolic_vector_2d('n2')

    bases_p = bezier_bases(2, t)
    bases_n = bezier_bases(1, t)

    points = (p1, p2, p3)
    normals = (n1, n2)

    p = make_bezier(points, bases_p)(t)
    n = make_bezier(normals, bases_n)(t)

    viewdir = p - view_point

    solution = viewdir.dot(n)**2
    
    # poly = to_polynomial(solution, t)
    # print("Got polynomial of degree: " + str(poly.degree()))

    solution = diff(solution, t)
    solution = simplify(solution)

    common, exprs = cse(solution, numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_projected_2d():
    '''

    Like above, but for an aligned slice of a 3d patch.
    Edges, and middles.

    Can still transform to 2d plane and solve there,
    but now need full normal information around the
    slice curve.

    --

    Edit:

    Here, I wrongfully assumed a single 2x2 patch could 
    cache normals of a 3x3 quadratic patch.

    '''

    u, v = symbols('u v')
    v = 0

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')],
    ]
    patch[0][0] = Matrix([0,0,0])
    patch[0][1][2] = 0
    patch[0][2][1] = 0
    patch[0][2][2] = 0

    bases = bezier_bases(2, u)
    points = (patch[0][0], patch[0][1], patch[0][2])
    p = make_bezier(points, bases)(u)

    patch_n = [
        [symbolic_vector_3d('normal1'), symbolic_vector_3d('normal2')],
        [symbolic_vector_3d('normal2'), symbolic_vector_3d('normal3')],
    ]
    normal = make_bezier_patch_with_points(patch_n, u, v)
    normal[2] = 0

    viewpos = symbolic_vector_3d('viewpos')
    viewdir = p - viewpos
    viewdir[2] = 0

    solution = viewdir.dot(normal)

    pprint(solution)

    '''
    Let's see. We now have dot(bezier, bezier)
    bezier * bezier + bezier * bezier
    '''

    # pprint(solution)

    # poly = to_polynomial(solution, u)
    # print("Got polynomial of degree: " + str(poly.degree()))

    # solution = solveset(solution, u)

    # common, exprs = cse(solution, numbered_symbols('a'))
    # print_code(common, exprs)

def silhouette_quadratic_3d_gradient():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]
   
    pos = make_bezier_patch_with_points(patch, u, v)

    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]
 
    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)**2

    partial_u = diff(solution, u)
    partial_v = diff(solution, v)

    common, exprs = cse((partial_u, partial_v), numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_3d_gradient_wrt_embedded_cubic():
    u, v, t = symbols('u v t')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]
   
    pos = make_bezier_patch_with_points(patch, u, v)

    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]
 
    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)**2

    uv1 = symbolic_vector_2d('uv1')
    uv2 = symbolic_vector_2d('uv2')
    uv3 = symbolic_vector_2d('uv3')
    uv4 = symbolic_vector_2d('uv4')
    uvs = (uv1, uv2, uv3, uv4)

    bases_uv = bezier_bases(3, t)
    uv = make_bezier(uvs, bases_uv)(t)

    solution = solution.subs(u, uv[0]).subs(v, uv[1])

    partials = []
    for p in uvs:
        partials.append(diff(solution, p))
        # partials.append(diff(solution, p[0]))
        # partials.append(diff(solution, p[1]))

    common, exprs = cse(partials, numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_3d_gradient_wrt_embedded_cubic_tangents():
    u, v, t, tuv2, tuv3 = symbols('u v t tuv2 tuv3')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]
   
    pos = make_bezier_patch_with_points(patch, u, v)

    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]
 
    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)**2

    uv1 = symbolic_vector_2d('uv1')
    uv4 = symbolic_vector_2d('uv4')
    uv2 = uv1 + symbolic_vector_2d('duv2') * tuv2
    uv3 = uv4 + symbolic_vector_2d('duv3') * tuv3
    
    uvs = (uv1, uv2, uv3, uv4)
    bases_uv = bezier_bases(3, t)
    uv = make_bezier(uvs, bases_uv)(t)

    solution = solution.subs(u, uv[0]).subs(v, uv[1])

    partials = [
        diff(solution, tuv2),
        diff(solution, tuv3)
    ]

    common, exprs = cse(partials, numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_3d_gradient_2nd():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]
   
    pos = make_bezier_patch_with_points(patch, u, v)

    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]
 
    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewdir = pos - viewpos

    '''
    Idea: derive 2nd order gradient in order to converge way faster

    Todo: this involves the Hessian, which I need to study up on.
    '''

    solution = viewdir.dot(normal)**2

    partial_u = diff(solution, u)
    partial_v = diff(solution, v)
    partial_uu = diff(partial_u, u)
    partial_vv = diff(partial_v, v)

    common, exprs = cse((partial_u, partial_uu, partial_v, partial_vv), numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_3d_edge():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]

    pos = make_bezier_patch_with_points(patch, u, v)
 
    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]

    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewpos = Matrix([0,0,0])
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)
    solution = Eq(solution, 0)
    solution = solution.subs(v, 0)
    solution = expand(solution)

    poly = to_polynomial(solution, u)
    print("Got polynomial of degree: %i"%poly.degree())

    solution_horned = horner(solution, wrt=u)
    print("Horner solution:")
    print(solution_horned)

    # 5th degree poly in u, so we need to find more clever workarounds

    # solution = solve(solution, u)

    # common, exprs = cse(solution, numbered_symbols('a'))
    # for (i,t) in enumerate(common):
    #     print("%i, %s"%(i, t))

    # for (i, t) in enumerate(exprs):
    #     print("%i, %s" % (i, t))
    
    # common, exprs = cse(solution, numbered_symbols('a'))
    # print_code(common, exprs)

def silhouette_quadratic_3d_edge_u():
    u = symbols('u')

    edge = (symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3'))

    bases_1st = bezier_bases(1, u)
    bases_2nd = bezier_bases(2, u)
    pos = make_bezier(edge, bases_2nd)(u)
 
    edge_du = (symbolic_vector_3d('du1'), symbolic_vector_3d('du2'))
    edge_dv = (symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3'))

    tangents_u = make_bezier(edge_du, bases_1st)(u)
    tangents_v = make_bezier(edge_dv, bases_2nd)(u)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewpos = Matrix([0,0,0])
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)
    solution = expand(solution)
    solution = simplify(solution)
    # poly = to_polynomial(solution, )
    pprint(solution)
    '''
    Turns out that writing this specialized version yields the same
    result as evaluating the full thing with v=0 substitution. 
    
    Sympy is cool like that. :)
    '''

    # poly = to_polynomial(solution, u)
    # print("Got polynomial of degree: %i"%poly.degree())
  

def silhouette_quadratic_3d_doubledot():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]

    pos = make_bezier_patch_with_points(patch, u, v)
 
    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]
    # patch_du = differentiate_patch_points_u(patch)
    # patch_dv = differentiate_patch_points_v(patch)

    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewpos = Matrix([0,0,0])
    viewdir = pos - viewpos

    '''
    Idea: We're looking for condition where dot(viewdir,normal)=0. How else can we get
    at this information?

    We have 2 orthogonal tangents ready to go. If we dot them both, and add them, and
    they sum to zero, we know the viewdir is in the tangent plane. Which is what
    we wanted to know, except now there's no cross product.

    Edit: Ah, of course the identity that we're looking to use is sin^2+cos^2 = 1
    This changes things. 6th degree polynomial. Nevermind.

    There's some more identities to play with, see if anything simpler pops out.

    '''
    dota = viewdir.dot(tangents_u)
    dotb = viewdir.dot(tangents_v)

    dot = dota**2 - dotb**2
    solution = Eq(dot, 0)
    solution = solution.subs(v, 0)
    solution = expand(solution)
    print(solution)

def silhouette_quadratic_3d_quadratic_normals():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]

    pos = make_bezier_patch_with_points(patch, u, v)
 
    patch_n = [
        [symbolic_vector_3d('n1'), symbolic_vector_3d('n2')],
        [symbolic_vector_3d('n3'), symbolic_vector_3d('n4')],
        [symbolic_vector_3d('n5'), symbolic_vector_3d('n6')]
    ]

    normal = make_bezier_patch_with_points(patch_n, u, v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewpos = Matrix([0,0,0])
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)
    solution = Eq(solution, 0)
    solution = solution.subs(v, 0)
    solution = expand(solution)
    print(solution)

    poly = to_polynomial(solution, u)
    print("Got polynomial of degree: %i"%poly.degree())
    # 4th degree poly in u, so we need to find more clever workarounds

    # solution = solve(solution, u)

    # common, exprs = cse(solution, numbered_symbols('a'))
    # for (i,t) in enumerate(common):
    #     print("%i, %s"%(i, t))

    # for (i, t) in enumerate(exprs):
    #     print("%i, %s" % (i, t))
    
    # common, exprs = cse(solution, numbered_symbols('a'))
    # print_code(common, exprs)

def silhouette_quadratic_3d_edge_gradient():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]

    pos = make_bezier_patch_with_points(patch, u, v)
 
    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]
    # patch_du = differentiate_patch_points_u(patch)
    # patch_dv = differentiate_patch_points_v(patch)

    tangents_u = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v = make_bezier_patch_with_points(patch_dv, u, v)

    normal = tangents_u.cross(tangents_v)

    viewpos = symbolic_vector_3d('viewPoint')
    viewdir = pos - viewpos

    solution = viewdir.dot(normal)**2
    solution = solution.subs(v, 0)
    partial_u = diff(solution, u)

    common, exprs = cse(partial_u, numbered_symbols('a'))
    print_code(common, exprs)

def silhouette_quadratic_3d_homogeneous_edge():
    '''
    Here, we assume we already have tangents in projective
    space available to us.
    '''
    u, v = symbols('u v')
 
    patch_du = [
        [symbolic_vector_3d('du1'), symbolic_vector_3d('du2')],
        [symbolic_vector_3d('du3'), symbolic_vector_3d('du4')],
        [symbolic_vector_3d('du5'), symbolic_vector_3d('du6')]
    ]
    patch_dv = [
        [symbolic_vector_3d('dv1'), symbolic_vector_3d('dv2'), symbolic_vector_3d('dv3')],
        [symbolic_vector_3d('dv4'), symbolic_vector_3d('dv5'), symbolic_vector_3d('dv6')]
    ]

    tangents_u_hom = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v_hom = make_bezier_patch_with_points(patch_dv, u, v)

    tangents_u = Matrix((tangents_u_hom / tangents_u_hom[3])[0:3])
    tangents_v = Matrix((tangents_v_hom / tangents_v_hom[3])[0:3])

    # pprint(tangents_u[0])

    normal = tangents_u.cross(tangents_v)

    pprint(normal[2])

    # solution = Eq(normal[2], 0)
    # solution = solution.subs(v, 0)
    # solution = solve(solution, u)
    # print(solution)

    # common, exprs = cse(solution, numbered_symbols('a'))
    # print_code(common, exprs)

def silhouette_quadratic_3d_homogeneous_edge_explicit():
    '''
    Here, we assume we already have tangents in projective
    space available to us.
    '''
    u, v = symbols('u v')
 
    patch = [
        [symbolic_vector_4d('p1'), symbolic_vector_4d('p2'), symbolic_vector_4d('p3')],
        [symbolic_vector_4d('p4'), symbolic_vector_4d('p5'), symbolic_vector_4d('p6')],
        [symbolic_vector_4d('p7'), symbolic_vector_4d('p8'), symbolic_vector_4d('p9')]
    ]

    patch_du = differentiate_patch_points_u(patch)
    patch_dv = differentiate_patch_points_v(patch)

    tangents_u_hom = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v_hom = make_bezier_patch_with_points(patch_dv, u, v)

    tangents_u = Matrix((tangents_u_hom / tangents_u_hom[3])[0:3])
    tangents_v = Matrix((tangents_v_hom / tangents_v_hom[3])[0:3])

    normal = tangents_u.cross(tangents_v)

    solution = Eq(normal[2], 0)
    solution = solution.subs(v, 0)
    poly = to_polynomial(solution, u)
    print("Degree: " + str(poly.degree()))
    # Yields a degree 3 polynomial

    # solution = solve(solution, u)
    # print(solution)

    # common, exprs = cse(solution, numbered_symbols('a'))
    # print_code(common, exprs)

def silhouette_quadratic_3d_homogeneous_edge_explicit_gradient():
    '''
    Here, we assume we already have tangents in projective
    space available to us.
    '''
    u, v = symbols('u v')
 
    patch = [
        [symbolic_vector_4d('p1'), symbolic_vector_4d('p2'), symbolic_vector_4d('p3')],
        [symbolic_vector_4d('p4'), symbolic_vector_4d('p5'), symbolic_vector_4d('p6')],
        [symbolic_vector_4d('p7'), symbolic_vector_4d('p8'), symbolic_vector_4d('p9')]
    ]

    patch[0][0] = Matrix([0,0,0,1])

    patch_du = differentiate_patch_points_u(patch)
    patch_dv = differentiate_patch_points_v(patch)

    tangents_u_hom = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v_hom = make_bezier_patch_with_points(patch_dv, u, v)

    tangents_u = Matrix((tangents_u_hom / tangents_u_hom[3])[0:3])
    tangents_v = Matrix((tangents_v_hom / tangents_v_hom[3])[0:3])

    normal = tangents_u.cross(tangents_v)

    solution = normal[2]**2
    solution = solution.subs(v, 0)
    partial_u = diff(solution, u)

    common, exprs = cse(partial_u, numbered_symbols('a'))
    print_code(common, exprs)

def quadratic_patch_3d_normals():
    u, v = symbols('u v')

    patch = [
        [symbolic_vector_3d('p1'), symbolic_vector_3d('p2'), symbolic_vector_3d('p3')],
        [symbolic_vector_3d('p4'), symbolic_vector_3d('p5'), symbolic_vector_3d('p6')],
        [symbolic_vector_3d('p7'), symbolic_vector_3d('p8'), symbolic_vector_3d('p9')]
    ]

    pos = make_bezier_patch_with_points(patch, u, v)

    tangent_u = diff(pos, u)
    tangent_v = diff(pos, v)
    normal = tangent_u.cross(tangent_v)

    # pprint(normal[0])

    common, exprs = cse(normal, numbered_symbols('a'))
    print_code(common, exprs)

def quadratic_2d_bezier():
    symbs = symbols('t, p1, p2, p3')
    t, p1, p2, p3 = symbs

    pd1 = 2 * (p2 - p1)
    pd2 = 2 * (p3 - p2)

    bases = bezier_bases(2, t)
    bases_d = bezier_bases(1, t)

    p = make_bezier((p1, p2, p3), bases)(t)
    pd = make_bezier((pd1, pd2), bases_d)(t)

    common, exprs = cse(p, numbered_symbols('a'))
    print("Point:")
    print_code(common, exprs)

    common, exprs = cse(pd, numbered_symbols('a'))
    print("Tangent:")
    print_code(common, exprs)

def quartic_bezier_3d():
    '''
    Todo: Find biquadratic forms
    '''
    t = symbols('t')

    p1, p2, p3, p4, p5 = symbols('p1 p2 p3 p4 p5')
    # p2 = symbolic_vector_3d('p2')
    # p3 = symbolic_vector_3d('p3')
    # p4 = symbolic_vector_3d('p4')
    # p5 = symbolic_vector_3d('p5')

    points = (p1, p2, p3, p4, p5)
    bases = bezier_bases(4, t)
    p = make_bezier(points, bases)(t)

    common, exprs = cse(p, numbered_symbols('a'))
    print_code(common, exprs)

    # points_d = get_curve_point_deltas(points, 4)
    # bases_d = bezier_bases(3, t)
    # pd = make_bezier_expr(points_d, bases_d)

def diagonal_of_linear_patch():
    '''
    Essentially: bilinear interpolation across a quad
    We can see that we get solution that is quadratic in t
    
    todo:

    solve for a line through 2 points on surface

    Reformulate points as offsets from p1

    I still have a notion that rewriting should yield diagonal
    line as linear combination as the bottom left and
    top right control point. For this flat quad, anyway.
    '''
    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')

    patch = [
        [p1, p2],
        [p3, p4]
    ]

    u, v, t = symbols('u v t')

    patch = make_bezier_patch_with_points(patch, u, v)

    patch.subs(p2, p2 - p1).subs(p3, p3 - p1)
    patch = patch.subs(v, t).subs(u, t)

    pprint(patch)

def diagonal_of_quadratic_patch():
    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')
    p5 = symbolic_vector_3d('p5')
    p6 = symbolic_vector_3d('p6')
    p7 = symbolic_vector_3d('p7')
    p8 = symbolic_vector_3d('p8')
    p9 = symbolic_vector_3d('p9')

    patch = [
        [p1, p2, p3],
        [p4, p5, p6],
        [p7, p8, p9]
    ]

    u, v, t = symbols('u v t')

    patch = make_bezier_patch_with_points(patch, u, v)
    p = patch.subs(v, t).subs(u, t)

    # pprint(p)

    common, exprs = cse(p, numbered_symbols('a'))

    print_code(common, exprs)

def line_inside_quadratic_patch():
    '''
    Yield a quartic curve that spans between two
    points on the surface defined in uv space

    In uv space itself, the line is a linear combination,
    i.e. a first order bezier curve.
    '''

    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')
    p5 = symbolic_vector_3d('p5')
    p6 = symbolic_vector_3d('p6')
    p7 = symbolic_vector_3d('p7')
    p8 = symbolic_vector_3d('p8')
    p9 = symbolic_vector_3d('p9')

    # Todo: would express with abstract symbols, but patch need linalg.
    # Resulting polynomials are identical over [x,y,z] though
    # p1, p2, p3, p4, p5, p6, p7, p8, p9 = symbols('p1, p2, p3, p4, p5, p6, p7, p8, p9')

    patch = [
        [p1, p2, p3],
        [p4, p5, p6],
        [p7, p8, p9]
    ]

    u, v, t = symbols('u v t')

    patch = make_bezier_patch_with_points(patch, u, v)

    # now re-express u,v as linear functions of single param t
    uv1 = symbolic_vector_2d('uv1')
    uv2 = symbolic_vector_2d('uv2')
    bases_uv = bezier_bases(1, t)
    uv = make_bezier((uv1, uv2), bases_uv)(t)

    p = patch.subs(u, uv[0]).subs(v, uv[1])

    common, exprs = cse(p, numbered_symbols('a'))

    # for expr in exprs:
    #     pprint(expr)

    print_code(common, exprs)

def quadratic_curve_on_quadratic_patch():
    '''
    As of 23-06-19, this looks like the fullest description for a
    full silhouette curve on a quadratic patch.

    It's an 8th degree polynomial... Hoo boy.

    Of course, there's clever ways to evaluate it. It's just that
    these babies might get unwieldy?
    '''

    u, v, t = symbols('u v t')

    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')
    p5 = symbolic_vector_3d('p5')
    p6 = symbolic_vector_3d('p6')
    p7 = symbolic_vector_3d('p7')
    p8 = symbolic_vector_3d('p8')
    p9 = symbolic_vector_3d('p9')

    patch = [
        [p1, p2, p3],
        [p4, p5, p6],
        [p7, p8, p9]
    ]
    patch = make_bezier_patch_with_points(patch, u, v)
    
    uv1 = symbolic_vector_2d('uv1')
    uv2 = symbolic_vector_2d('uv2')
    uv3 = symbolic_vector_2d('uv3')

    bases_uv = bezier_bases(2, t)
    uv = make_bezier((uv1, uv2, uv3), bases_uv)(t)

    p = patch.subs(u, uv[0]).subs(v, uv[1])

    # poly = to_polynomial(p[0], t)
    # print("Got polynomial of degree: " + str(poly.degree()))

    common, exprs = cse(p, numbered_symbols('a'))
    print_code(common, exprs)

def cubic_curve_on_quadratic_patch():
    u, v, t = symbols('u v t')

    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')
    p5 = symbolic_vector_3d('p5')
    p6 = symbolic_vector_3d('p6')
    p7 = symbolic_vector_3d('p7')
    p8 = symbolic_vector_3d('p8')
    p9 = symbolic_vector_3d('p9')

    patch = [
        [p1, p2, p3],
        [p4, p5, p6],
        [p7, p8, p9]
    ]
    patch = make_bezier_patch_with_points(patch, u, v)
    
    uv1 = symbolic_vector_2d('uv1')
    uv2 = symbolic_vector_2d('uv2')
    uv3 = symbolic_vector_2d('uv3')
    uv4 = symbolic_vector_2d('uv4')

    bases_uv = bezier_bases(3, t)
    uv = make_bezier((uv1, uv2, uv3, uv4), bases_uv)(t)

    p = patch.subs(u, uv[0]).subs(v, uv[1])

    poly = to_polynomial(p[0], t)
    print("Got polynomial of degree: " + str(poly.degree()))

    common, exprs = cse(p, numbered_symbols('a'))
    print_code(common, exprs)

    # The above difference yields [0,0,0], so this checks out

def make_point_list(count, basis):
    points = []
    for i in range(0, count):
        points.append(symbolic_vector('p%i'%i, basis))
    return points

def make_patch(u, v, degree_u, degree_v, basis):
    dimensions = len(basis)
 
    bases_u = bezier_bases(degree_u, u)
    bases_v = bezier_bases(degree_v, v)

    pos = Matrix([0]*dimensions)

    i = 1
    for vIdx in range(0, degree_v+1):
        for uIdx in range(0, degree_u+1):
            pos += symbolic_vector('p_%i' % i, basis) * bases_u[uIdx] * bases_v[vIdx]
            i+=1

    return pos

def line_inside_quartic_patch():
    u, v, t = symbols('u v t')

    patch = make_patch(u, v, 4, 4, BASIS_3D)

    # now re-express u,v as linear functions of single param t
    uv1 = symbolic_vector_2d('uv1')
    uv2 = symbolic_vector_2d('uv2')
    bases_uv = bezier_bases(1, t)
    uv = make_bezier((uv1, uv2), bases_uv)(t)

    p = patch.subs(u, uv[0]).subs(v, uv[1])

    common, exprs = cse(p, numbered_symbols('a'))

    print_code(common, exprs)

def prove_patch_derives():
    u, v = symbols('u v')

    p1 = symbolic_vector_3d('p1')
    p2 = symbolic_vector_3d('p2')
    p3 = symbolic_vector_3d('p3')
    p4 = symbolic_vector_3d('p4')
    p5 = symbolic_vector_3d('p5')
    p6 = symbolic_vector_3d('p6')
    p7 = symbolic_vector_3d('p7')
    p8 = symbolic_vector_3d('p8')
    p9 = symbolic_vector_3d('p9')

    patch = [
        [p1, p2, p3],
        [p4, p5, p6],
        [p7, p8, p9]
    ]
    pos = make_bezier_patch_with_points(patch, u, v)

    patch_du = differentiate_patch_points_u(patch)
    patch_dv = differentiate_patch_points_v(patch)
    tangents_u_a = make_bezier_patch_with_points(patch_du, u, v)
    tangents_v_a = make_bezier_patch_with_points(patch_dv, u, v)

    tangents_u_b = diff(pos, u)
    tangents_v_b = diff(pos, v)

    print("\nTangents A:\n")
    pprint(expand(tangents_u_a[0]))
    print("\nTangents B:\n")
    pprint(expand(tangents_u_b[0]))
    print("\nDifference:\n")
    difference = expand(tangents_u_a[0]) - expand(tangents_u_b[0])
    pprint(difference)

    normals_a = tangents_u_a.cross(tangents_v_a)
    normals_b = tangents_u_b.cross(tangents_v_b)

    print("\nNormals A:\n")
    pprint(normals_a[0])
    print("\nNormals B:\n")
    pprint(normals_b[0])
    print("\nDifference:\n")
    difference = expand(normals_a[0]) - expand(normals_b[0])
    pprint(difference)

def main():
    init_printing(pretty_print=True, use_unicode=True, num_columns=180)

    line_inside_quartic_patch()

if __name__ == "__main__":
    main()
