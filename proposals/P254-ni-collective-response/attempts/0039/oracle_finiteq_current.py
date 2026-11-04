#!/usr/bin/env python3
"""A39 focused oracle: literal real-space gauged NC scalar, not a material run.

Reference changes the load-bearing representation: tensor spherical/radial
nodes directly integrate the gauged two-point kernel. Its straight-line
integral is exp(i q.midpoint)*sinc(q.separation/2), not production line nodes,
shifted projectors, solid harmonics derivatives, or continuity reconstruction.
Old radial/uniform proof is carried as a bound receipt, never rerun here.
Root executes; this source itself asserts no numerical success.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0028'), str(HERE.parent/'0027')]

import numpy as np
from numpy.polynomial.legendre import leggauss
from run_physical_spatial_oracle import model, momenta
from substrate_framework.material_physical_operators import physical_hamiltonian, PhysicalFieldVertices
from substrate_framework.material_finiteq_current import MaterialFiniteQCurrent, MODEL
from substrate_framework.qe_bloch_edge import BlochEdge


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=str(len(data)), sha256=hashlib.sha256(data).hexdigest())


def implementations():
    names = ('material_finiteq_current.py', 'material_finiteq_projectors.py',
             'material_physical_operators.py', 'qe_material_field_derivatives.py',
             'qe_bloch_finiteq.py', 'qe_bloch_edge.py', 'qe_bloch_finiteq_ionic.py',
             'qe_material_hamiltonian.py', 'qe_nc_projectors.py')
    return [binding(HERE/'oracle_finiteq_current.py'),
            binding(HERE/'export_finiteq_current_material.py'),
            *[binding(ROOT/'src/substrate_framework'/name) for name in names],
            binding(HERE.parent/'0028/run_physical_spatial_oracle.py'),
            binding(HERE.parent/'0028/run_field_derivative_oracle.py'),
            binding(HERE.parent/'0027/bloch_oracle_fixture.py')]


def compare(records, name, actual, expected, tolerance=3e-10):
    value = float(np.linalg.norm(actual-expected)/max(1., np.linalg.norm(expected)))
    records[name] = value
    if not np.isfinite(value) or value > tolerance:
        raise AssertionError((name, value, tolerance))


def expose(records, name, mutant, expected, floor=1e-7):
    value = float(np.linalg.norm(mutant-expected)/max(1., np.linalg.norm(expected)))
    records[name] = value
    if not np.isfinite(value) or value <= floor:
        raise AssertionError(('unexposed mutation', name, value))


def realspace_nodes(table):
    """Literal retained Simpson radial measure and independent sphere harmonics."""
    mu, w = leggauss(14)
    phi = np.arange(28)*2*np.pi/28
    z, az = np.meshgrid(mu, phi, indexing='ij')
    n = np.column_stack(((np.sqrt(1-z*z)*np.cos(az)).ravel(),
                         (np.sqrt(1-z*z)*np.sin(az)).ravel(), z.ravel()))
    x, y, z = n.T
    c0, c1, c2 = 1/np.sqrt(4*np.pi), np.sqrt(3/(4*np.pi)), np.sqrt(5/(4*np.pi))
    angular = np.column_stack((np.full(len(n), c0), c1*z, -c1*x, -c1*y,
        c2*(3*z*z-1)/2, -c2*np.sqrt(3)*x*z, -c2*np.sqrt(3)*y*z,
        c2*np.sqrt(3)*(x*x-y*y)/2, c2*np.sqrt(3)*x*y))
    aw = np.repeat(w, 28)*2*np.pi/28
    upf = table.upf
    rw = upf.rab*np.array([1., 4., 2., 4., 1.])/3
    points, amplitudes = [], []
    for j, radius in enumerate(upf.r):
        points.append(radius*n)
        amplitudes.append(angular[:, upf.channel_lm]*aw[:, None]
            *(rw[j]*radius*upf.beta_r[j, upf.channel_radial]/np.sqrt(table.omega_bohr3)))
    return np.concatenate(points), np.concatenate(amplitudes)


def literal_nonlocal(inputs, edge, a, *, epsilon=None, via_atom=False):
    """Dense two-point kernel, bounded real-space row blocks, for small oracle only.

    epsilon=None differentiates the actual scalar. Otherwise evaluate the
    scalar exp(-i epsilon integral A.dl), allowing an independent finite
    difference. via_atom is the distinct two-segment path mutation.
    """
    table = inputs.tables['X']
    points, amplitude = realspace_nodes(table)
    ps, pt = momenta(edge.source_basis), momenta(edge.target_basis)
    q = edge.q_cart_bohr_inv
    result = np.zeros((len(pt), len(ps)), complex)
    for tau in inputs.reference.positions_bohr.T:
        r = points+tau
        source_plane = np.exp(1j*r @ ps.T)
        target_plane = np.exp(-1j*pt @ r.T)
        for start in range(0, len(r), 64):
            sl = slice(start, start+64)
            separation = r[sl, None, :]-r[None, :, :]
            midpoint = (r[sl, None, :]+r[None, :, :])/2
            line = ((separation @ a)*np.exp(1j*midpoint @ q)
                    *np.sinc((separation @ q)/(2*np.pi)))
            if via_atom:
                # r' -> tau -> r, each segment integrated analytically.
                left = r[sl]-tau
                right = tau-r
                line = ((left @ a)[:, None]*np.exp(1j*(tau+left/2) @ q)[:, None]
                        *np.sinc((left @ q)/(2*np.pi))[:, None]
                        +(right @ a)[None, :]*np.exp(1j*(r+right/2) @ q)[None, :]
                        *np.sinc((right @ q)/(2*np.pi))[None, :])
            kernel = amplitude[sl] @ table.upf.expanded_d @ amplitude.conj().T
            kernel *= -1j*line if epsilon is None else np.exp(-1j*epsilon*line)
            result += target_plane[:, sl] @ (kernel @ source_plane)
    return result


def literal_kinetic(edge, a, *, fold=None):
    ps, pt = momenta(edge.source_basis), momenta(edge.target_basis)
    fold = edge.reciprocal_fold if fold is None else fold
    matrix = np.zeros((len(pt), len(ps)), complex)
    for ti, mt in enumerate(edge.target_basis.miller.T):
        for si, ms in enumerate(edge.source_basis.miller.T):
            if np.array_equal(mt+fold, ms):
                matrix[ti, si] = (pt[ti]+ps[si]) @ a
    return matrix


def run():
    residuals, mutants = {}, {}
    with TemporaryDirectory(prefix='p254-a39-oracle-') as directory:
        inputs, grid, background = model(directory)
        frames = [grid.load_frame(next(s for s in grid.slots if s.ispin == 1
                  and s.mesh_index == mesh)) for mesh in ((0, 0, 1), (0, 0, 2))]
        source, target = frames
        edge = BlochEdge.create(inputs.reference, source.basis, target.basis, (0., 0., .25))
        if edge.source_dimension == edge.target_dimension or not np.any(edge.reciprocal_fold):
            raise AssertionError('unequal PW grids and nonzero reciprocal fold required')
        sh, th = [physical_hamiltonian(background, f.basis, 1, block_size=2) for f in frames]
        owner = MaterialFiniteQCurrent(edge, sh, th, line_order=16, block_size=2)
        # q is along z: this genuinely transverse complex direction has no z leg.
        a = np.array([.23+.11j, -.17+.19j, 0.])
        action = owner.prepare(a)
        kinetic = literal_kinetic(edge, a)
        nonlocal_matrix = literal_nonlocal(inputs, edge, a)
        dense = kinetic+nonlocal_matrix
        rng = np.random.default_rng(2543901)
        x = rng.normal(size=(edge.source_dimension, 3))+1j*rng.normal(size=(edge.source_dimension, 3))
        y = rng.normal(size=(edge.target_dimension, 3))+1j*rng.normal(size=(edge.target_dimension, 3))
        compare(residuals, 'transverse_gauged_scalar_forward', action.apply(x), dense @ x)
        compare(residuals, 'independent_true_adjoint', action.adjoint_apply(y), dense.conj().T @ y)
        compare(residuals, 'hilbert_pair', y.conj().T @ action.apply(x), action.adjoint_apply(y).conj().T @ x)
        compare(residuals, 'current_Ry_to_Ha_sign', action.current_apply(x), -.5*dense @ x)
        compare(residuals, 'current_adjoint_units', action.current_adjoint_apply(y), -.5*dense.conj().T @ y)
        compare(residuals, 'column_partition', action.apply(x[:, :1]), (dense @ x)[:, :1])
        compare(residuals, 'line_order_refinement', MaterialFiniteQCurrent(edge, sh, th,
            line_order=8).prepare(a).apply(x), dense @ x)
        # Exact gauge identity for the *unprojected* scalar at the retained endpoint
        # rows. Do not commute two endpoint-truncated H matrices: that drops
        # intermediate momenta and is not this identity.
        q = edge.q_cart_bohr_inv
        longitudinal = owner.prepare(q).apply(x)
        # For nonlocal commutator, shifted momenta must be evaluated even outside
        # endpoint support; the direct kernel implements V exp(iq.r)-exp(iq.r) V.
        table = inputs.tables['X']
        points, amplitudes = realspace_nodes(table)
        ps, pt = momenta(source.basis), momenta(target.basis)
        commutator = np.zeros_like(dense)
        for tau in inputs.reference.positions_bohr.T:
            r = points+tau
            bt = np.exp(-1j*pt @ r.T) @ amplitudes
            bs = np.exp(-1j*ps @ r.T) @ amplitudes
            btq = np.exp(-1j*(pt-q) @ r.T) @ amplitudes
            bsq = np.exp(-1j*(ps+q) @ r.T) @ amplitudes
            commutator += bt @ table.upf.expanded_d @ bsq.conj().T-btq @ table.upf.expanded_d @ bs.conj().T
        # q.(d_pt+d_ps)V integrates to V(pt,ps+q)-V(pt-q,ps).
        expected_long = literal_kinetic(edge, q)+commutator
        compare(residuals, 'longitudinal_continuity', longitudinal, expected_long @ x)
        diagonal = BlochEdge.create(inputs.reference, source.basis, source.basis, (0., 0., 0.))
        uniform = MaterialFiniteQCurrent(diagonal, sh, sh).prepare(a)
        vertices = PhysicalFieldVertices(sh, projector_block_size=2)
        old = vertices.prepare(a.real).apply(x)+1j*vertices.prepare(a.imag).apply(x)
        compare(residuals, 'analytic_uniform_existing_vertex', uniform.apply(x), old)
        uniform_dense = literal_kinetic(diagonal, a)+literal_nonlocal(inputs, diagonal, a)
        compare(residuals, 'uniform_gauged_scalar_reference', uniform.apply(x), uniform_dense @ x)
        expose(mutants, 'omit_both_nonlocal_legs', kinetic @ x, dense @ x)
        # Each omitted leg is independently integrated from real-space projector
        # amplitudes; these are not arbitrary fractions of the total vertex.
        for label in ('target', 'source'):
            # A leg-only production mutation is constructed from independently
            # integrated shifted projector factors below, not from half the sum.
            from scipy.integrate import quad_vec
            def leg(t):
                u, v = pt-t*q, ps+(1-t)*q
                total = np.zeros_like(dense)
                for tau in inputs.reference.positions_bohr.T:
                    btu = np.exp(-1j*u @ points.T) @ amplitudes
                    bsv = np.exp(-1j*v @ points.T) @ amplitudes
                    du = np.exp(-1j*u @ points.T) @ ((-1j*(points @ a))[:, None]*amplitudes)
                    dv = np.exp(-1j*v @ points.T) @ ((-1j*(points @ a.conj()))[:, None]*amplitudes)
                    phase = np.exp(-1j*(pt[:, None, :]-ps[None, :, :]-q) @ tau)
                    total += phase*(du @ table.upf.expanded_d @ bsv.conj().T if label == 'target'
                                    else btu @ table.upf.expanded_d @ dv.conj().T)
                return total
            one_leg = quad_vec(leg, 0., 1., epsabs=1e-11, epsrel=1e-11)[0]
            expose(mutants, 'omit_'+label+'_leg', (dense-one_leg) @ x, dense @ x)
        expose(mutants, 'reverse_current_sign', .5*dense @ x, -.5*dense @ x)
        expose(mutants, 'conjugate_polarization', owner.prepare(a.conj()).apply(x), dense @ x)
        expose(mutants, 'transpose_instead_of_dagger', dense.T @ y, dense.conj().T @ y)
        expose(mutants, 'wrong_fold', (literal_kinetic(edge, a, fold=-edge.reciprocal_fold)+nonlocal_matrix) @ x, dense @ x)
        expose(mutants, 'wrong_source_target_phase', (kinetic+nonlocal_matrix.conj()) @ x, dense @ x)
        expose(mutants, 'two_segment_transverse_path', (kinetic+literal_nonlocal(inputs, edge, a, via_atom=True)) @ x, dense @ x)
        finite_errors = []
        for step in (.002, .001):
            finite = (literal_nonlocal(inputs, edge, a, epsilon=step)
                      -literal_nonlocal(inputs, edge, a, epsilon=-step))/(2*step)
            finite_errors.append(float(np.linalg.norm(finite-nonlocal_matrix)))
        if finite_errors[-1] > 2e-7 or finite_errors[-1] >= finite_errors[0]/3:
            raise AssertionError(('gauged scalar finite difference refinement', finite_errors))
        return dict(residuals=residuals, mutations=mutants, finite_difference_errors=finite_errors,
                    dimensions=[edge.source_dimension, edge.target_dimension], fold=edge.reciprocal_fold.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    proof = HERE.parent/'0028/physical-spatial-oracle-003.json'
    old = json.loads(proof.read_text())
    if old.get('status') != 'passed':
        raise ValueError('existing admitted physical radial/field proof required')
    result = dict(schema='P254-A39-finiteq-current-oracle-v1', model=MODEL,
        implementation_sources=implementations(), carried_radial_uniform_proof=binding(proof),
        evidence_scope='new finite-q vertex only; no material response or continuum accuracy', **run(), status='passed')
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'], result=binding(args.output))))


if __name__ == '__main__':
    main()
