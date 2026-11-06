#!/usr/bin/env python3
"""Independent exact/literal finite-frame oracles, not production material data.

SymPy checks are independent of implementation algebra. Literal coordinate
waves and four-operation 4^3 enumeration expose sign/fold/phase/coverage errors.
No original WFC is loaded, no Native process is launched.
"""
import argparse
from dataclasses import replace
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from types import SimpleNamespace, MappingProxyType

import numpy as np
import sympy as sp

from substrate_framework.qe_bloch_grid import UNITARY_POLICY, build_bloch_grid
from substrate_framework.qe_real_bloch_grid import (
    REAL_POLICY, RealHamiltonianScope, build_real_bloch_grid,
    validate_real_coverage, transport_coefficients, transform_real_vector,
)
from substrate_framework.material_reference_state import prepare_source
from substrate_framework.material_reference_transport import raw_state_transport
from substrate_framework.zero_temperature_preparation import (
    ModeLevel, FillingPolicy, fill_zero_temperature,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rejects(call, message):
    try:
        call()
    except (ValueError, TypeError):
        return True
    raise AssertionError('mutation accepted: ' + message)


def exact_algebra():
    i = sp.I
    psi = sp.Matrix([[1, 1+i], [i, 2], [1-i, -i]]) / 5
    d = sp.Matrix([[4, i], [-i, 3]]) / 20
    permutation = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    u = sp.diag(1, i, -1) * permutation
    target = u * psi.conjugate()
    gram, rho = psi.H * psi, psi * d * psi.H
    target_gram = target.H * target
    target_rho = target * d.conjugate() * target.H
    zero = lambda m: all(sp.simplify(x) == 0 for x in m)
    require(zero(target_gram - gram.conjugate()), 'exact Gram conjugation')
    require(zero(target_rho - u * rho.conjugate() * u.H), 'exact D conjugation')
    require(not zero(target_gram - gram), 'Gram fixture must detect missing conjugation')
    require(not zero(target * d * target.H - target_rho), 'D fixture must detect missing conjugation')
    require(sp.simplify(sp.trace(d.conjugate() * target_gram) - sp.conjugate(sp.trace(d * gram))) == 0,
            'exact unweighted number invariant')
    h = sp.Matrix([[2, 1, 2], [1, 3, -1], [2, -1, 4]])
    x = sp.diag(0, 1, 3)
    # Real nonlocal hopping: its full commutator current is not just kinetic p.
    current = i * (h*x - x*h)
    force = sp.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 3]])
    force_target = u * force.conjugate() * u.H
    current_target = -u * current.conjugate() * u.H
    fm = sp.simplify(sp.trace(rho * force))
    jm = sp.simplify(sp.trace(rho * current))
    require(fm != 0 and jm != 0, 'force/current sign fixture is sensitive')
    require(sp.simplify(sp.trace(target_rho * force_target) - sp.conjugate(fm)) == 0, 'force even')
    require(sp.simplify(sp.trace(target_rho * current_target) + sp.conjugate(jm)) == 0, 'current odd')
    require(sp.simplify(sp.trace(target_rho * (-current_target)) + sp.conjugate(jm)) != 0,
            'wrong physical-current parity detected')
    # Production symmetric-frame preparation compared to a directly recomputed
    # target Gram/D; this is NOT the independent exact reference above.
    c = np.array(psi, dtype=complex)
    unitary = np.array(u, dtype=complex)
    state = prepare_source(c, [.2, .8], source_identity='literal-source',
                           basis_identity='literal-basis', authority='finite-oracle-only')
    mapped = raw_state_transport(state, antiunitary=True)
    direct = prepare_source(unitary @ c.conj(), [.2, .8], source_identity='literal-source',
                            basis_identity='literal-target', authority='finite-oracle-only')
    require(np.allclose(mapped.gram, direct.gram, atol=1e-14, rtol=0), 'production Gram transport')
    require(np.allclose(mapped.density_coefficients, direct.density_coefficients, atol=1e-12, rtol=0),
            'production D transport')
    require(np.array_equal(mapped.populations, state.populations), 'real populations unchanged')
    return dict(force_mean=str(fm), full_nonlocal_current_mean=str(jm),
        gram_exact=True, density_exact=True, number_exact=True,
        sensitive_mutations=['missing Gram conjugation', 'missing D conjugation', 'wrong current sign'],
        symmetric_frame_d_max_error=float(np.max(np.abs(mapped.density_coefficients-direct.density_coefficients))))


def literal_wave_oracle():
    # Non-self-inverse rotation, nonzero fractional translation and boundary fold.
    s = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    ks = sp.Matrix([-sp.Rational(1, 2), sp.Rational(1, 4), 0])
    kt = sp.Matrix([-sp.Rational(1, 4), 0, -sp.Rational(1, 2)])
    g = sp.Matrix([0, 0, -1])
    f = sp.Matrix([sp.Rational(1, 4), 0, 0])
    m = sp.Matrix([[0, 1, -1], [0, -1, 0], [0, 0, 1]])
    coeff = sp.Matrix([[1, 1+sp.I], [sp.I, 2], [1-sp.I, -sp.I]]) / 5
    require(kt == -s*ks+g, 'independent nonzero reciprocal fold')
    target_m = -s*m - g * sp.ones(1, m.cols)
    x = sp.Matrix(sp.symbols('x y z', real=True))
    # Compare full Bloch exponent at each original row, not just squared norm.
    for row in range(m.cols):
        left = (target_m[:, row]+kt).dot(x) + (m[:, row]+ks).dot(f)
        right = -(m[:, row]+ks).dot(s.T*x-f)
        require(sp.expand(left-right) == 0, 'exact full Bloch exponent')
        require((target_m[:, row]+kt) == -s*(m[:, row]+ks), 'exact physical momentum')
    sn, kn, fn, mn = (np.array(v, dtype=float) for v in (s, ks, f, m))
    kn, fn = kn.ravel(), fn.ravel()
    target_rows = np.array(target_m, dtype=int)
    c = np.array(coeff, dtype=complex)
    slot = SimpleNamespace(orbital_conjugation=True, epsilon=-1,
        operation=SimpleNamespace(rotation=sn, fractional_translation=fn), source_k_fractional=kn)
    actual_c = transport_coefficients(slot, mn, c)
    omega = 8.
    points = np.array([[1/8, 1/4, 3/8], [1/4, 1/8, 1/2], [3/8, 1/2, 1/8]])

    def evaluate(rows, k, columns, coordinates, normalization=np.sqrt(omega)):
        return np.exp(2j*np.pi * coordinates @ (rows + k[:, None])) @ columns / normalization

    # Independent literal pullback and conjugation of FULL real-space waves.
    expected = evaluate(mn, kn, c, points @ sn - fn).conj()
    observed = evaluate(target_rows, np.array(kt, dtype=float).ravel(), actual_c, points)
    error = float(np.max(np.abs(observed-expected)))
    require(error < 2e-14, 'literal transported full Bloch waves')
    mutations = {
        'missing coefficient conjugation': evaluate(target_rows, np.array(kt, dtype=float).ravel(),
             np.exp(2j*np.pi*((mn.T+kn)@fn))[:, None]*c, points),
        'wrong reciprocal fold': evaluate(-sn@mn, np.array(kt, dtype=float).ravel(), actual_c, points),
        'wrong translation phase': evaluate(target_rows, np.array(kt, dtype=float).ravel(),
             np.exp(-2j*np.pi*((mn.T+kn)@fn))[:, None]*c.conj(), points),
        'wrong volume normalization': evaluate(target_rows, np.array(kt, dtype=float).ravel(), actual_c, points, omega),
    }
    gaps = {key: float(np.max(np.abs(value-expected))) for key, value in mutations.items()}
    require(all(value > 1e-3 for value in gaps.values()), 'all wave mutations must be sensitive')
    # Literal uniform Fourier quadrature independently establishes Omega/Ngrid,
    # NO k-star weight and NO spin degeneracy in a single-sector Gram.
    mesh = np.array(list(product(range(8), repeat=3)), dtype=float) / 8
    samples = evaluate(mn, kn, c, mesh)
    quadrature_gram = omega/len(mesh) * samples.conj().T @ samples
    require(np.max(np.abs(quadrature_gram-c.conj().T@c)) < 1e-14, 'volume-normalized Fourier Gram')
    require(np.max(np.abs(2*quadrature_gram-c.conj().T@c)) > .01, 'false spin factor detected')
    return dict(literal_wave_max_error=error, mutation_gaps=gaps,
                fold=list(g), exact_exponents=True, normalization='Omega/Ngrid; 1/sqrt(Omega) orbitals')


def mesh_fixture():
    """Literal dispersed crystal group, not imported production slot records."""
    matrices = [((1,0,0),(0,1,0),(0,0,1)), ((1,0,0),(0,-1,0),(0,0,-1)),
                ((1,0,0),(0,0,-1),(0,-1,0)), ((1,0,0),(0,0,1),(0,1,0))]
    all_indices = list(product(range(4), repeat=3))

    def apply(matrix, h, sign=1):
        return tuple(sign*sum(matrix[a][b]*h[b] for b in range(3)) % 4 for a in range(3))

    # Independent exact integer orbit partition, in QE original mesh order.
    unassigned, representatives, orbits = set(all_indices), [], []
    for h in all_indices:
        if h not in unassigned:
            continue
        orbit = {apply(matrix, h, sign) for matrix in matrices for sign in (1, -1)}
        require(orbit <= unassigned, 'literal crystal orbit partition')
        representatives.append(h)
        orbits.append(orbit)
        unassigned -= orbit
    require(not unassigned and len(representatives) == 21, 'literal dispersed21 source count')
    fold = lambda h: np.array([v/4 if v < 2 else v/4-1 for v in h])
    points = tuple(SimpleNamespace(ik=k+1, npw=3, weight=len(orbit)/64,
                                  k_cart_bohr_inv=fold(h))
                   for k, (h, orbit) in enumerate(zip(representatives, orbits)))
    operations = tuple(SimpleNamespace(index=j+1, rotation=np.array(m),
        fractional_translation=np.zeros(3), equivalent_atoms=np.array([1]), effective_time_reversal=False)
        for j, m in enumerate(matrices))
    ref = SimpleNamespace(source={'path':'literal.xml', 'bytes':'1', 'sha256':'a'*64},
        lsda=True, noncolin=False, spinorbit=False, gamma_only=False, nat=1,
        cell_bohr=np.eye(3), positions_bohr=np.zeros((3,1)), atom_names=('Ni',),
        fft_smooth=(4,4,4), reciprocal_bohr_inv=np.eye(3), nbnd_by_spin=(2,2),
        nks=21, kpoints=points, monkhorst_pack=dict(nk1=4,nk2=4,nk3=4,k1=0,k2=0,k3=0),
        saved_symmetries=SimpleNamespace(nsym=4, operations=operations))
    ref.select = lambda ik, spin: (points[ik-1], np.array([ik+100*spin, ik+100*spin+.5]),
                                  np.array([-.02, 1.02]))
    bindings = {(spin,p.ik):dict(path=f'literal-spin{spin}-k{p.ik}',bytes=1,
                                sha256=sha256(f'{spin}:{p.ik}'.encode()).hexdigest())
                for spin in (1,2) for p in points}
    scope = RealHamiltonianScope(MappingProxyType(ref.source), (), 'literal-oracle-scope-not-material')
    grid = build_real_bloch_grid(ref, bindings, scope=scope, symmetry_policy=REAL_POLICY)
    require(len(grid.slots) == 128, 'all128 explicit mesh/spin slots')
    for slot in grid.slots:
        require(slot.mesh_index in orbits[slot.source_ik-1], 'independent literal source orbit')
        require(slot.weight == 1/64, 'one full-k measure per explicit spin')
        require(slot.source_binding['sha256'] == bindings[slot.ispin,slot.source_ik]['sha256'],
                'same-spin source binding; no physical TR spin swap')
        require(slot.eigenvalues_ry[0] == slot.source_ik+100*slot.ispin, 'same-spin labels')
    selected = next(s for s in grid.slots if s.mesh_index == (3,0,0) and s.ispin == 1)
    require(selected.orbital_conjugation, 'minus-quarter-x requires conjugation')
    spin_swapped = next(s for s in grid.slots if s.mesh_index == (3,0,0) and s.ispin == 2)
    require(np.max(np.abs(selected.eigenvalues_ry-spin_swapped.eigenvalues_ry)) == 100,
            'spin swap negative control must be sensitive')
    rejects(lambda: build_bloch_grid(ref, bindings, symmetry_policy=UNITARY_POLICY), 'unitary-only false coverage')
    rejects(lambda: validate_real_coverage(ref, grid.slots[:-1]), 'missing slot')
    rejects(lambda: validate_real_coverage(ref, grid.slots[:-1]+(grid.slots[0],)), 'duplicate slot')
    rejects(lambda: validate_real_coverage(ref, tuple(replace(s, weight=1/32) for s in grid.slots)), 'double spin measure')
    rejects(lambda: build_real_bloch_grid(ref, dict(list(bindings.items())[:-1]),
            scope=scope, symmetry_policy=REAL_POLICY), 'missing original source')
    rejects(lambda: build_real_bloch_grid(ref, bindings, scope=scope, symmetry_policy=UNITARY_POLICY), 'forged unitary authority')
    # Force/current parity plus a nonidentity spatial rotation are exercised.
    vector_slot = next(s for s in grid.slots if s.orbital_conjugation and s.operation.index != 1)
    v = np.array([1.,2.,3.])
    expected = np.array(matrices[vector_slot.operation.index-1]) @ v
    require(np.array_equal(transform_real_vector(grid, vector_slot, v, kind='current'), -expected), 'odd current rotation')
    require(np.array_equal(transform_real_vector(grid, vector_slot, v[None,:], kind='force'), expected[None,:]), 'even force rotation')
    # Independent neutral filling utility sensitivity, not MV interpretation.
    policy = FillingPolicy('literal saved labels', 'neutral mean N', 'explicit spin slots 1/64',
                           'exact degeneracy only', 'Ry', 'electrons/cell')
    levels = [ModeLevel(f'{s.ispin}:{s.mesh_index}:{b}', s.ispin, float(e), 1/64)
              for s in grid.slots for b,e in enumerate(s.eigenvalues_ry)]
    filling = fill_zero_temperature(levels, 1.5, policy=policy)
    require(abs(filling.weighted_number-1.5) < 1e-14, 'weighted neutral population')
    require(all(0 <= p <= 1 for p in filling.populations), 'Pauli population bound')
    require(sum(filling.populations) != filling.weighted_number, 'weights are measures, not probabilities')
    require(any(p not in (-.02,1.02) for p in filling.populations), 'raw MV deliberately not used')
    return dict(original_sources_per_spin=21, slots=128,
        conjugated_slots=sum(s.orbital_conjugation for s in grid.slots),
        exact_orbit_sizes=[len(o) for o in orbits], source_weight_sum=sum(p.weight for p in points),
        weighted_number=filling.weighted_number,
        rejected=['unitary-only coverage', 'missing slot', 'duplicate slot', 'double spin weight',
                  'missing source binding', 'unitary authority substitution'], spin_swap_gap=100)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = dict(schema='control-preparation-independent-oracle-v1',
        exact_finite_frame=exact_algebra(), literal_bloch_wave=literal_wave_oracle(),
        literal_mesh_population=mesh_fixture(), status='pass',
        scope='independent exact and literal finite constructions; not material data/convergence')
    # SymPy integers are serialized as ordinary integers, never invented floats.
    report['literal_bloch_wave']['fold'] = [int(x) for x in report['literal_bloch_wave']['fold']]
    here = Path(__file__).resolve()
    registration = here.with_name('control-preparation-registration.json')
    report['sources'] = [dict(path=str(p), bytes=p.stat().st_size,
        sha256=sha256(p.read_bytes()).hexdigest()) for p in (here, registration)]
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=report['status'], output=str(args.output.resolve()))))


if __name__ == '__main__':
    main()
