#!/usr/bin/env python3
"""Independent finite UPF/Fourier matrix oracle; explicitly run with --output.

Uses actual physical_finiteq, PhysicalFieldVertices and MaterialHamiltonian
objects built on a tiny synthetic finite geometry. References use literal
spherical/radial quadrature and explicit observable/density matrices, not
production projector jets/actions. PW-square matrices exist ONLY here.
No propagation, actual Ni result, or full-material qualification is claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / 'src'), str(HERE), str(HERE.parent / '0027'),
                str(HERE.parent / '0025'), str(HERE.parent / '0024')]

import numpy as np
from mixed_frame_response import DensityBlock, DensityDyad, density_dyads
from mixed_material_moments import _one_body, response_blocks
from mixed_physical_field_dyads import physical_field_dyad_moments
from run_field_derivative_oracle import literal_integral
from run_physical_spatial_oracle import model, momenta, dense_ionic
from substrate_framework.material_physical_operators import (
    PhysicalFieldVertices, physical_hamiltonian, physical_finiteq)
from substrate_framework.qe_bloch_edge import BlochEdge


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def source_bindings():
    paths = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.is_relative_to(ROOT):
                paths.add(path)
    return [binding(p) for p in sorted(paths)]


def compare(report, name, actual, expected, tolerance=3e-10):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape:
        raise AssertionError((name, actual.shape, expected.shape))
    defect = float(np.max(np.abs(actual - expected), initial=0.))
    scale = max(1., float(np.max(np.abs(expected), initial=0.)))
    report[name] = defect / scale
    if not np.isfinite(defect) or defect > tolerance * scale:
        raise AssertionError((name, defect, scale))


def expose(report, name, mutant, expected, floor=1e-8):
    difference = float(np.linalg.norm(np.asarray(mutant) - np.asarray(expected)))
    scale = max(1., float(np.linalg.norm(expected)))
    report[name] = difference / scale
    if not np.isfinite(difference) or difference <= floor * scale:
        raise AssertionError(('unexposed mutation', name, difference, scale))


def refuse(report, name, operation):
    try:
        operation()
    except (ValueError, TypeError):
        report[name] = 'refused'
    else:
        raise AssertionError('invalid boundary admitted: ' + name)


def literal_endpoint(inputs, background, basis):
    p = momenta(basis)
    table, ref = inputs.tables['X'], inputs.reference
    x = np.indices(ref.fft_smooth).reshape(3, -1, order='F').T / np.asarray(ref.fft_smooth)
    plane = np.exp(2j * np.pi * x @ basis.miller)
    h0 = np.diag(np.sum(p * p, axis=1)).astype(complex)
    h0 += plane.conj().T @ (background.vrs[:, basis.slot.ispin-1, None] * plane) / len(x)
    field = np.array([np.diag(2 * p[:, j]).astype(complex) for j in range(3)])
    kinetic = field.copy()
    center = literal_integral(table, p, np.eye(3)[0], np.eye(3)[0])[0]
    left_terms, right_terms = np.zeros_like(field), np.zeros_like(field)
    diagonal_channel = np.zeros_like(field)
    for tau in ref.positions_bohr.T:
        phase = np.exp(-1j * p @ tau)[:, None]
        beta = center * phase
        d = table.upf.expanded_d
        h0 += beta @ d @ beta.conj().T
        for axis in range(3):
            first = literal_integral(table, p, np.eye(3)[axis], np.eye(3)[axis])[1] * phase
            left_terms[axis] += first @ d @ beta.conj().T
            right_terms[axis] += beta @ d @ first.conj().T
            dd = np.diag(np.diag(d))
            diagonal_channel[axis] += first @ dd @ beta.conj().T + beta @ dd @ first.conj().T
    field += left_terms + right_terms
    return NS(p=p, center=center, h0=h0, field=field, kinetic=kinetic,
              left_terms=left_terms, right_terms=right_terms,
              diagonal_channel=kinetic + diagonal_channel)


def dense_block(inputs, edge, block, endpoints, diagonal):
    source, target = endpoints[block.source_endpoint], endpoints[block.target_endpoint]
    rho = sum(d.weight * d.left @ d.D @ d.right.conj().T for d in block.dyads)
    ref = inputs.reference
    x = np.indices(ref.fft_smooth).reshape(3, -1, order='F').T / np.asarray(ref.fft_smooth)
    ps = np.exp(2j * np.pi * x @ edge.source_basis.miller)
    pt = np.exp(2j * np.pi * x @ (edge.target_basis.miller + edge.reciprocal_fold[:, None]))
    density = np.zeros((len(x), 2), complex)
    omega = abs(float(np.linalg.det(ref.cell_bohr)))
    density[:, edge.ispin-1] = np.sum((pt @ rho) * ps.conj(), axis=1) / omega
    gradient = []
    for atom in range(ref.nat):
        for axis in range(3):
            unit = np.zeros((ref.nat, 3)); unit[atom, axis] = 1
            derivative = dense_ionic(inputs, edge, unit, (source.center, target.center))
            gradient.append(np.trace(derivative.conj().T @ rho))
    return dict(density=density, gradient=np.asarray(gradient),
                field_gradient_electronic=np.asarray([np.trace(a @ rho) for a in source.field])
                    if diagonal else np.zeros(3, complex),
                one_body=np.trace(source.h0 @ rho) if diagonal else 0j)


def block_checks(inputs, grid, background, report, mutations, boundaries, rng):
    frames = [grid.load_frame(next(s for s in grid.slots
              if s.ispin == 1 and s.mesh_index == (0, 0, z))) for z in range(4)]
    hamiltonians = [physical_hamiltonian(background, f.basis, 1, block_size=2) for f in frames]
    endpoints = [literal_endpoint(inputs, background, f.basis) for f in frames]
    ranks, weights = (2, 3, 4, 2), (.13, .29, .37, .21)
    random = lambda shape: rng.normal(size=shape) + 1j*rng.normal(size=shape)
    psi = tuple(random((len(e.p), rank)) for e, rank in zip(endpoints, ranks))
    target_plus, target_minus = (1, 2, 3, 0), (3, 0, 1, 2)
    plus = tuple(random((len(endpoints[target_plus[i]].p), ranks[i])) for i in range(4))
    minus = tuple(random((len(endpoints[target_minus[i]].p), ranks[i])) for i in range(4))
    mixed = tuple(random(a.shape) for a in psi)
    vectors = tuple(random((rank, 1)) for rank in ranks)
    coefficients = tuple(v@v.conj().T for v in vectors)
    blocks = density_dyads(psi, plus, minus, mixed, coefficients,
                         target_plus=target_plus, target_minus=target_minus,
                         source_weights=weights)
    captured = None
    for component, collection in blocks.items():
        diagonal = component in ('mean', 'mixed')
        q = (0., 0., 0. if diagonal else (.25 if component == 'plus' else -.25))
        for i, block in enumerate(collection):
            s, t = block.source_endpoint, block.target_endpoint
            edge = BlochEdge.create(inputs.reference, frames[s].basis, frames[t].basis, q)
            coupling = physical_finiteq(inputs, edge, hamiltonians[s].factors,
                                       hamiltonians[t].factors, block_size=2)
            vertex = PhysicalFieldVertices(hamiltonians[s], projector_block_size=2)
            actual = response_blocks({component: block}, coupling=coupling,
                                     field_vertex=vertex if diagonal else None,
                                     diagonal=diagonal)[component]
            # Exercise the actual one-body primitive used by material energy,
            # independently of response-only factor reuse.
            actual['one_body'] = (sum(_one_body(d, hamiltonians[s].apply, 2)
                                      for d in block.dyads) if diagonal else 0j)
            expected = dense_block(inputs, edge, block, endpoints, diagonal)
            for key in expected:
                compare(report, f'{component}/{i}/{key}', actual[key], expected[key])
            if diagonal and component == 'mixed' and i == 1:
                captured = block, coupling, vertex, hamiltonians[s], expected
            if not diagonal:
                sizes = [len(e.p) for e in endpoints]
                offsets = np.cumsum([0]+sizes)
                full = np.zeros((sum(sizes), sum(sizes)), complex)
                rho = sum(d.weight*d.left@d.D@d.right.conj().T for d in block.dyads)
                full[offsets[t]:offsets[t+1], offsets[s]:offsets[s+1]] = rho
                hglobal = np.zeros_like(full)
                for k, endpoint in enumerate(endpoints):
                    hglobal[offsets[k]:offsets[k+1], offsets[k]:offsets[k+1]] = endpoint.h0
                compare(report, f'{component}/{i}/global_geometry',
                        np.trace(hglobal@full), actual['one_body'])
                global_field = []
                for axis in range(3):
                    observable = np.zeros_like(full)
                    for k, endpoint in enumerate(endpoints):
                        observable[offsets[k]:offsets[k+1], offsets[k]:offsets[k+1]] = endpoint.field[axis]
                    global_field.append(np.trace(observable@full))
                compare(report, f'{component}/{i}/global_field_geometry',
                        global_field, actual['field_gradient_electronic'])
                refuse(boundaries, f'{component}/{i}/forbidden_field_vertex',
                    lambda: response_blocks({'x': block}, coupling=coupling,
                                            field_vertex=vertex, diagonal=False))
                refuse(boundaries, f'{component}/{i}/false_diagonal',
                    lambda: response_blocks({'x': block}, coupling=coupling,
                                            field_vertex=None, diagonal=True))
    block, coupling, vertex, h0, expected = captured
    kwargs = dict(coupling=coupling, field_vertex=vertex, diagonal=True)
    refuse(boundaries, 'diagonal_missing_field', lambda: response_blocks(
        {'x': block}, **{**kwargs, 'field_vertex': None}))
    refuse(boundaries, 'diagonal_wrong_geometry', lambda: response_blocks(
        {'x': block}, **{**kwargs, 'diagonal': False}))
    refuse(boundaries, 'nonboolean_geometry', lambda: response_blocks(
        {'x': block}, **{**kwargs, 'diagonal': 1}))
    refuse(boundaries, 'empty_density_block', lambda: response_blocks(
        {'x': DensityBlock(1, 1, ())}, **kwargs))
    rho = sum(d.weight*d.left@d.D@d.right.conj().T for d in block.dyads)
    expose(mutations, 'discard_H0_offdiagonal_action',
           np.trace(np.diag(np.diag(endpoints[1].h0))@rho), expected['one_body'])
    for omitted in range(len(block.dyads)):
        reduced = DensityBlock(block.target_endpoint, block.source_endpoint,
                               tuple(d for j, d in enumerate(block.dyads) if j != omitted))
        mutant = response_blocks({'x': reduced}, **kwargs)['x']
        expose(mutations, 'omit_mixed_dyad_'+str(omitted), mutant['density'], expected['density'])
    for name, transform in (
            ('transpose_original_D', lambda d: DensityDyad(d.left, d.D.T, d.right, d.source_index, d.weight)),
            ('population_diagonal_D', lambda d: DensityDyad(d.left, np.diag(np.diag(d.D)), d.right, d.source_index, d.weight)),
            ('weight_twice', lambda d: DensityDyad(d.left, d.D, d.right, d.source_index, d.weight**2)),
            ('weight_omitted', lambda d: DensityDyad(d.left, d.D, d.right, d.source_index, 1.))):
        changed = DensityBlock(block.target_endpoint, block.source_endpoint,
                               tuple(transform(d) for d in block.dyads))
        expose(mutations, name, response_blocks({'x': changed}, **kwargs)['x']['density'],
               expected['density'])
    gamma = density_dyads((psi[1],), (random(psi[1].shape),),
                          (random(psi[1].shape),), (mixed[1],),
                          (coefficients[1],), target_plus=(0,), target_minus=(0,),
                          source_weights=(weights[1],))
    gamma_edge = BlochEdge.create(inputs.reference, frames[1].basis, frames[1].basis, (0., 0., 0.))
    grouped = {name: collection[0] for name, collection in gamma.items()}
    # All four independent groups share physical factors in ONE invocation.
    actual_groups = response_blocks(grouped, **kwargs)
    for name, density_block in grouped.items():
        reference = dense_block(inputs, gamma_edge, density_block, [endpoints[1]], True)
        for key in ('density', 'gradient', 'field_gradient_electronic'):
            compare(report, f'Gamma/grouped/{name}/{key}', actual_groups[name][key], reference[key])
    rectangular = {}
    for index, (nl, nr, weight) in enumerate(((2, 3, .23), (3, 2, .61))):
        dyad = DensityDyad(random((len(endpoints[1].p), nl)), random((nl, nr)),
                           random((len(endpoints[1].p), nr)), 0, weight)
        rectangular[str(index)] = DensityBlock(0, 0, (dyad,))
    rect_values = response_blocks(rectangular, **kwargs)
    for name, density_block in rectangular.items():
        reference = dense_block(inputs, gamma_edge, density_block, [endpoints[1]], True)
        for key in ('density', 'gradient', 'field_gradient_electronic'):
            compare(report, f'Gamma/rectangular/{name}/{key}', rect_values[name][key], reference[key])
    for density_block in blocks['mixed']:
        for d in density_block.dyads:
            if d.D is not coefficients[d.source_index] or d.weight != weights[d.source_index]:
                raise AssertionError('borrowed original coefficient/weight lost')
    return frames[1], hamiltonians[1], endpoints[1], random


def field_checks(inputs, frame, hamiltonian, endpoint, random, report, mutations, boundaries):
    vertex = PhysicalFieldVertices(hamiltonian, projector_block_size=2)
    left, right = random((len(endpoint.p), 3)), random((len(endpoint.p), 2))
    coefficient = random((3, 2))
    weight = .173
    rho = weight * left @ coefficient @ right.conj().T
    trace = lambda operators: np.array([np.trace(a @ rho) for a in operators])
    actual = physical_field_dyad_moments(vertex, left, right, coefficient, weight=weight)
    expected = trace(endpoint.field)
    compare(report, 'rectangular_nonHermitian_coefficient', actual, expected)
    for name, operators in (
            ('omit_left_projector_jet', endpoint.field - endpoint.left_terms),
            ('omit_right_projector_jet', endpoint.field - endpoint.right_terms),
            ('omit_UPF_offdiagonal_channels', endpoint.diagonal_channel),
            ('kinetic_only', endpoint.kinetic)):
        expose(mutations, name, trace(operators), expected)
    expose(mutations, 'real_projection', expected.real, expected)
    expose(mutations, 'conjugate_trace', expected.conj(), expected)
    expose(mutations, 'conjugate_coefficient', physical_field_dyad_moments(
        vertex, left, right, coefficient.conj(), weight=weight), expected)
    # Both coefficient absorption branches and non-divisible block widths.
    compare(report, 'rectangular_absorb_opposite_leg', physical_field_dyad_moments(
        vertex, right, left, coefficient.conj().T, weight=weight), expected.conj())
    zero = physical_field_dyad_moments(vertex, left, right, np.zeros_like(coefficient), weight=weight)
    compare(report, 'singular_zero_coefficient', zero, np.zeros(3))
    compare(report, 'zero_weight', physical_field_dyad_moments(
        vertex, left, right, coefficient, weight=0.), np.zeros(3))
    same = physical_field_dyad_moments(vertex, left, left, np.eye(3), weight=weight)
    compare(report, 'same_borrowed_frame', same, np.array([
        weight * np.trace(left.conj().T @ a @ left) for a in endpoint.field]))
    refuse(boundaries, 'wrong_field_vertex', lambda: physical_field_dyad_moments(
        hamiltonian, left, right, coefficient, weight=weight))
    refuse(boundaries, 'wrong_coefficient_shape', lambda: physical_field_dyad_moments(
        vertex, left, right, coefficient.T, weight=weight))
    refuse(boundaries, 'wrong_PW_rows', lambda: physical_field_dyad_moments(
        vertex, left[:-1], right, coefficient, weight=weight))
    refuse(boundaries, 'complex_source_weight', lambda: physical_field_dyad_moments(
        vertex, left, right, coefficient, weight=weight+1j))
    refuse(boundaries, 'nonfinite_coefficient', lambda: physical_field_dyad_moments(
        vertex, left, right, coefficient * np.nan, weight=weight))
    # Independent finite differences of literal shifted physical projectors.
    direction = np.array([.23, -.31, .19])
    table, p = inputs.tables['X'], endpoint.p
    def scalar(t):
        shifted = p + t * direction
        center = literal_integral(table, shifted, direction, direction)[0]
        operator = np.diag(np.sum(shifted**2, axis=1)).astype(complex)
        for tau in inputs.reference.positions_bohr.T:
            beta = center * np.exp(-1j * p @ tau)[:, None]
            operator += beta @ table.upf.expanded_d @ beta.conj().T
        return np.trace(operator @ rho)
    errors = []
    for step in (.008, .004, .002):
        value = (scalar(step) - scalar(-step)) / (2 * step)
        errors.append(float(abs(value - direction @ actual)))
    if not errors[-1] < errors[0] / 8 or errors[-1] > 2e-6:
        raise AssertionError(('literal field finite-difference refinement', errors))
    report['literal_field_difference_errors'] = errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = dict(schema='P254-mixed-material-moments-finite-oracle-v1',
                  scope='finite synthetic independent UPF/Fourier matrices only; no material propagation',
                  production_APIs={
                      'mixed_material_moments.response_blocks': '(groups, *, coupling, field_vertex, diagonal)',
                      'mixed_physical_field_dyads.physical_field_dyad_moments': '(vertex, left, right, coefficient, *, weight)'},
                  checks={}, mutations={}, boundaries={})
    rng = np.random.default_rng(254280503)
    with TemporaryDirectory(prefix='p254-mixed-material-moments-') as directory:
        inputs, grid, background = model(directory)
        frame, hamiltonian, endpoint, random = block_checks(
            inputs, grid, background, report['checks'], report['mutations'], report['boundaries'], rng)
        field_checks(inputs, frame, hamiltonian, endpoint, random,
                     report['checks'], report['mutations'], report['boundaries'])
    report['implementation_and_reference_sources'] = source_bindings()
    report['status'] = 'passed'
    text = json.dumps(report, indent=2, allow_nan=False) + '\n'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(text)
    print(json.dumps(dict(status='passed', output=binding(args.output))))


if __name__ == '__main__':
    main()
