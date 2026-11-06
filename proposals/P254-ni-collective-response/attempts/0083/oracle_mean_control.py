"""Exact new curvature/transpose block; original A29/A50 proofs are carried."""
import sympy as sp
from mean_control_blocks import gamma_tangent_increment, gamma_transpose_increment


def check():
    q0, q1, a, b, c, d = sp.symbols('q0 q1 a b c d')
    potential = (a*q0**3 + 3*b*q0**2*q1 + 3*c*q0*q1**2 + d*q1**3)/6
    curvature = sp.hessian(potential, (q0, q1))
    j = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1],
                   [-1, 0, 0, 0], [0, -1, 0, 0]])
    h = sp.zeros(4)
    h[:2, :2] = curvature
    x = sp.Matrix(sp.symbols('x0:4'))
    z = sp.Matrix(sp.symbols('z0:4'))
    expected = j*h*x
    forward = gamma_tangent_increment(j[:, :2], curvature, x[:2, :])
    dual = gamma_transpose_increment(j[:, :2], curvature, z)
    forward_defect = sp.simplify(forward-expected)
    transpose_defect = sp.expand((z.T*forward)[0]-(dual.T*x[:2, :])[0])
    if forward_defect != sp.zeros(4, 1) or transpose_defect != 0:
        raise AssertionError('new curvature map is not the scalar Hamiltonian tangent/true transpose')
    values = {q0: 2, q1: -1, a: 1, b: 2, c: 3, d: 4,
              **dict(zip(x, (1, 2, 3, 4))), **dict(zip(z, (5, 6, 7, 8)))}
    # A plausible wrong momentum/curvature sign must fail the actual identity.
    sign_mutation = sp.expand((z.T*(-forward))[0]-(dual.T*x[:2, :])[0]).subs(values)
    if sign_mutation == 0:
        raise AssertionError('sign mutation was not exposed')
    return dict(schema='p254-a83-new-block-symbolic-oracle-v1',
        scalar=str(potential), forward_defect=str(forward_defect),
        bilinear_transpose_defect=str(transpose_defect),
        sign_mutation_defect=str(sign_mutation),
        scope='Exact new additive Gamma curvature and nonconjugating transpose only; original full A29/A50 coupled tangent/CAR proofs carried, not replayed. No temporal or material accuracy certificate.')
