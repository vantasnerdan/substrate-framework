#!/usr/bin/env python3
"""Independent SymPy finite NC scalar oracle; no material files or DFT.

Run from any directory with the project dependencies available. The synthetic
finite basis has two off-diagonally coupled s projectors plus three p channels,
a nonzero atomic phase, a non-axis polarization and full PW complement actions.
This is an equation check, not a material current result.
"""
from pathlib import Path
import json
import sys
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3] / 'src'))

import numpy as np
import sympy as sp
from scipy.linalg import eigvalsh
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian, PWFactors, SpeciesProjectors
from substrate_framework.qe_material_current import MaterialOpticalCurrent
from substrate_framework.qe_nc_projectors import centered_projectors, place_projectors


def main():
    t, s = sp.symbols('t s', real=True)
    k = sp.Matrix([sp.Rational(1, 3), sp.Rational(1, 4), sp.Rational(1, 5)])
    e = sp.Matrix([sp.Rational(2, 3), -sp.Rational(1, 2), sp.Rational(1, 4)])
    f = sp.Matrix([sp.Rational(1, 5), sp.Rational(1, 3), -sp.Rational(1, 2)])
    tau = sp.Matrix([sp.Rational(1, 7), -sp.Rational(1, 8), sp.Rational(1, 9)])
    miller = sp.Matrix([[0, 1, -1], [0, 1, 0], [0, 0, 1]])
    d = sp.diag(sp.Rational(3, 5), sp.Rational(4, 5),
                sp.Rational(2, 7), sp.Rational(2, 7), sp.Rational(2, 7))
    d[0, 1] = d[1, 0] = sp.Rational(1, 4)
    momenta = [k+miller[:, i] for i in range(3)]
    # Independent Cartesian solid harmonics; no production interpolation,
    # harmonic routine, derivative stencil or product-rule helper is used.
    rows = []
    for p in momenta:
        q = p+t*e+s*f
        r2 = q.dot(q)
        phase = sp.exp(-sp.I*p.dot(tau))
        rows.append([phase*sp.exp(-r2/3)/sp.sqrt(4*sp.pi),
                     phase*(1+r2/5)*sp.exp(-r2/4)/sp.sqrt(4*sp.pi),
                     -sp.I*phase*sp.sqrt(3/(4*sp.pi))*q[2]*sp.exp(-r2/5),
                     sp.I*phase*sp.sqrt(3/(4*sp.pi))*q[0]*sp.exp(-r2/5),
                     sp.I*phase*sp.sqrt(3/(4*sp.pi))*q[1]*sp.exp(-r2/5)])
    beta = sp.Matrix(rows)
    symbolic_h = sp.diag(*[(p+t*e+s*f).dot(p+t*e+s*f) for p in momenta])+beta*d*beta.conjugate().T
    h0 = symbolic_h.subs({t: 0, s: 0})
    l = sp.diff(symbolic_h, t).subs({t: 0, s: 0})
    c = sp.diff(symbolic_h, t, 2).subs({t: 0, s: 0})
    mixed = sp.diff(symbolic_h, t, s).subs({t: 0, s: 0})
    numeric = lambda value: np.array(value.evalf(), dtype=np.complex128)
    real = lambda value: np.array(value, dtype=float).reshape(-1)

    class AnalyticRadial:
        # Explicit analytic finite-NC example, not a stand-in material source.
        upf = SimpleNamespace(angular_l=np.array([0, 0, 1]),
            channel_radial=np.array([0, 1, 2, 2, 2]),
            channel_lm=np.array([0, 0, 1, 2, 3]),
            channel_l=np.array([0, 0, 1, 1, 1]))

        def beta(self, radius):
            r2 = radius**2
            return np.column_stack((np.exp(-r2/3), (1+r2/5)*np.exp(-r2/4),
                                    radius*np.exp(-r2/5)))

    table = AnalyticRadial()
    p = np.array([real(point) for point in momenta])
    fft = PhysicalRhoFFT.create((5, 5, 5), np.array(miller, dtype=int), np.eye(3), 1.)
    placed = place_projectors(centered_projectors(table, p), p, real(tau))
    inputs = SimpleNamespace(sources=(), tables={'X': table},
                             reference=SimpleNamespace(positions_bohr=real(tau)[:, None]))
    background = SimpleNamespace(inputs=inputs, fft=fft, vrs=np.zeros((fft.nnr, 2)))
    factors = PWFactors(fft, np.sum(p*p, axis=1),
                        (SpeciesProjectors('X', (1,), placed, numeric(d)),),
                        real(k), None, ({'kind': 'analytic-NC-oracle'},))
    operator = MaterialOpticalCurrent(MaterialHamiltonian(background, factors, 1, block_size=2))
    identity = np.eye(3, dtype=complex)
    expected_l, expected_c = numeric(l), numeric(c)
    np.testing.assert_allclose(operator.field_action(identity, [0, 0, 0]), numeric(h0), atol=2e-13)
    np.testing.assert_allclose(operator.field_action(identity, .037*real(e)),
                               numeric(symbolic_h.subs({t: sp.Rational(37, 1000), s: 0})), atol=2e-13)
    refinements = []
    for h, linear, contact in operator.refinement_actions(identity, real(e), step_bohr_inv=.004):
        refinements.append(dict(step_bohr_inv=h,
            linear_error=float(np.linalg.norm(linear-expected_l)),
            contact_error=float(np.linalg.norm(contact-expected_c))))
        np.testing.assert_allclose(linear, linear.conj().T, atol=3e-13)
        np.testing.assert_allclose(contact, contact.conj().T, atol=3e-13)
    for name in ('linear_error', 'contact_error'):
        assert refinements[1][name] < .4*refinements[0][name], refinements
        assert refinements[2][name] < .4*refinements[1][name], refinements
    np.testing.assert_allclose(linear, expected_l, atol=3e-7, rtol=3e-7)
    np.testing.assert_allclose(contact, expected_c, atol=3e-7, rtol=3e-7)
    np.testing.assert_allclose(operator.mixed_contact_action(identity, real(e), real(f), step_bohr_inv=.0005),
                               numeric(mixed), atol=5e-7, rtol=5e-7)

    # Scalar derivatives checked independently, not merely matrix agreement.
    vector = sp.Matrix([1, sp.I/2, sp.Rational(2, 3)])
    scalar = (vector.conjugate().T*symbolic_h*vector)[0]
    v = numeric(vector).reshape(-1)
    for order, action in ((1, linear), (2, contact)):
        exact_scalar = complex(sp.diff(scalar, t, order).subs({t: 0, s: 0}).evalf())
        np.testing.assert_allclose(np.vdot(v, action @ v), exact_scalar, atol=1e-6)

    kinetic_l = np.diag(2*(p @ real(e)))
    kinetic_c = 2*float(e.dot(e))*identity
    beta0 = numeric(beta.subs({t: 0, s: 0}))
    beta1 = numeric(sp.diff(beta, t).subs({t: 0, s: 0}))
    beta2 = numeric(sp.diff(beta, t, 2).subs({t: 0, s: 0}))
    diagonal_d = np.diag(np.diag(numeric(d)))
    mutations = {
        'omitted_nonlocal_current': kinetic_l,
        'one_nonlocal_leg': kinetic_l+beta1 @ numeric(d) @ beta0.conj().T,
        'wrong_electron_field_sign': -expected_l,
        'Hartree_Ry_factor_two': expected_l/2,
        'discarded_offdiagonal_channels': kinetic_l+beta1 @ diagonal_d @ beta0.conj().T+beta0 @ diagonal_d @ beta1.conj().T,
    }
    exposed = {name: float(np.linalg.norm(value-expected_l)) for name, value in mutations.items()}
    exposed.update(omitted_nonlocal_contact=float(np.linalg.norm(kinetic_c-expected_c)),
        omitted_cross_contact=float(np.linalg.norm(2*beta1 @ numeric(d) @ beta1.conj().T)),
        omitted_second_projector_legs=float(np.linalg.norm(beta2 @ numeric(d) @ beta0.conj().T+beta0 @ numeric(d) @ beta2.conj().T)),
        omitted_all_contact=float(np.linalg.norm(expected_c)))
    assert min(exposed.values()) > 1e-6, exposed

    # Exact free kinetic/mode Hessian: positivity is NOT inferred by removing A².
    pp, aa = sp.symbols('p a', real=True)
    full_scalar = (pp+aa)**2+aa**2/2
    dropped_scalar = pp**2+2*pp*aa+aa**2/2
    full_hessian = sp.hessian(full_scalar, (pp, aa))
    dropped_hessian = sp.hessian(dropped_scalar, (pp, aa))
    assert full_hessian.det() > 0 and dropped_hessian.det() < 0
    assert eigvalsh(numeric(full_hessian))[0] > 0
    assert eigvalsh(numeric(dropped_hessian))[0] < 0
    print(json.dumps(dict(status='passed', refinement=refinements, exposed_mutations=exposed,
        full_free_field_hessian_det=str(full_hessian.det()),
        missing_A2_hessian_det=str(dropped_hessian.det()),
        scope='independent finite analytic NC example only; not a material result'), indent=2))


if __name__ == '__main__':
    main()
