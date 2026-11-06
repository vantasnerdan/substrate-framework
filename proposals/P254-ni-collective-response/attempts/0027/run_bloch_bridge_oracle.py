"""NEW finite A27 oracle source; run only under Root's later source-bound grant.

Literal PW sums, a QE-convention FFT pullback, and finite differences of a
separately assembled full-space NC scalar operator are independent comparators.
No material decoding, historical replay, or physical-output claim occurs here.
"""
from dataclasses import replace
from itertools import product
import json
import tempfile
import numpy as np
from bloch_oracle_fixture import fixture, source_rows
from substrate_framework.qe_bloch_grid import build_bloch_grid, iter_bloch_frames, qe_fold, UNITARY_POLICY
from substrate_framework.qe_bloch_frame import transport_bloch_frame, admit_source
from substrate_framework.qe_bloch_qzero import BlochQZeroBridge
from substrate_framework.qe_bloch_ionic_action import BlochIonicAction
from substrate_framework.qe_material_hamiltonian import build_pw_factors_for_basis, MaterialHamiltonian
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.qe_saved_symmetry import parse_saved_symmetries


def close(a, b, tolerance=2e-10):
    actual, expected = np.asarray(a), np.asarray(b)
    if actual.shape != expected.shape or not (np.isfinite(actual).all() and np.isfinite(expected).all()):
        raise AssertionError('comparison requires equal shapes and finite operands')
    error = float(np.max(np.abs(actual-expected), initial=0))
    if error > tolerance:
        raise AssertionError(f'independent residual {error} exceeds {tolerance}')
    return error


def sensitive(correct, wrong, minimum=1e-7):
    expected, mutant = np.asarray(correct), np.asarray(wrong)
    if expected.shape != mutant.shape or not (np.isfinite(expected).all() and np.isfinite(mutant).all()):
        raise AssertionError('mutation requires equal shapes and finite operands')
    error = float(np.max(np.abs(expected-mutant), initial=0))
    if error <= minimum:
        raise AssertionError(f'mutation not exposed: {error}')
    return error


def rejected(call):
    try:
        call()
    except ValueError:
        return
    raise AssertionError('malformed convention was accepted')


def coordinates(shape):
    return np.array(list(product(*(range(n) for n in shape))))/shape


def literal_wave(miller, coefficients, x):
    return np.exp(2j*np.pi*(x @ miller)) @ coefficients


def source_fft_pullback(source, slot, shape):
    """Separate scatter/ifft/grid permutation/fold/fft path, no row transport."""
    cube = np.zeros((*shape, source.nbnd), dtype=complex)
    cube[tuple(source.miller[i] % shape[i] for i in range(3))] = source.coefficients
    waves = np.fft.ifftn(cube, axes=(0,1,2))*np.prod(shape)
    x = coordinates(shape)
    arg = x @ slot.operation.rotation - slot.operation.fractional_translation
    indices = np.rint(arg*shape).astype(int) % shape
    pulled = waves[tuple(indices[:, i] for i in range(3))]
    pulled *= np.exp(-2j*np.pi*(x @ slot.reciprocal_shift))[:, None]
    pulled *= np.exp(-2j*np.pi*(slot.source_k_fractional @ slot.operation.fractional_translation))
    return pulled


def grid_oracle(inputs, bindings, root, rotations):
    ref = inputs.reference
    grid = build_bloch_grid(ref, bindings, symmetry_policy=UNITARY_POLICY)
    assert len(grid.slots) == 128
    assert len({(s.ispin, s.mesh_index) for s in grid.iter_slots()}) == 128
    for operation, rotation in zip(ref.saved_symmetries.operations, rotations):
        close(operation.rotation, rotation, 0)
        assert operation.raw_time_reversal is None and operation.effective_time_reversal is False
    info = root.find('output/symmetries/symmetry/info')
    info.set('time_reversal', 'false')
    explicit = parse_saved_symmetries(root, ref.source, ref.nat, noncolin=False)
    assert explicit.operations[0].raw_time_reversal is False
    assert explicit.operations[0].effective_policy == 'explicit XML attribute'
    del info.attrib['time_reversal']
    wrong_weight = replace(ref.kpoints[0], weight=2*ref.kpoints[0].weight)
    rejected(lambda: build_bloch_grid(replace(ref, kpoints=(wrong_weight, *ref.kpoints[1:])),
                                     bindings, symmetry_policy=UNITARY_POLICY))
    incomplete = dict(bindings); del incomplete[2, 10]
    rejected(lambda: build_bloch_grid(ref, incomplete, symmetry_policy=UNITARY_POLICY))
    asymmetric = next(s for s in ref.saved_symmetries.operations if not np.array_equal(s.rotation, s.rotation.T))
    sensitive(asymmetric.rotation, asymmetric.rotation.T)
    close(qe_fold([.5,-.5,1.5,-1.5]), [-.5,.5,-.5,.5], 0)
    sensitive(qe_fold([.5]), np.array([.5])-np.rint([.5]))
    slot = next(s for s in grid.slots if np.any(s.reciprocal_shift) and s.ispin == 2)
    source = load_wavefunctions(slot.source_binding)
    frame = grid.load_frame(slot)
    x = coordinates(ref.fft_smooth)
    target = literal_wave(frame.basis.miller, frame.coefficients, x)
    pulled = literal_wave(source.miller, source.coefficients, x @ slot.operation.rotation)
    pulled *= np.exp(-2j*np.pi*(x @ slot.reciprocal_shift))[:, None]
    residual = close(target, pulled)
    close(target, source_fft_pullback(source, slot, ref.fft_smooth))
    pulled_cube = source_fft_pullback(source, slot, ref.fft_smooth).reshape((*ref.fft_smooth, source.nbnd))
    transformed = np.fft.fftn(pulled_cube, axes=(0,1,2))/np.prod(ref.fft_smooth)
    close(transformed[tuple(frame.basis.miller[i] % ref.fft_smooth[i] for i in range(3))],
          frame.coefficients)
    sensitive(target, literal_wave(slot.operation.rotation @ source.miller, source.coefficients, x))
    assert frame.coefficients is source.coefficients or np.array_equal(frame.coefficients, source.coefficients)
    assert {tuple(m) for m in frame.basis.miller.T} == {tuple(m) for m in source_rows(slot.k_fractional, 2.3).T}
    for spin in (1,2):
        own = next(s for s in grid.slots if s.source_ik == 2 and s.ispin == spin and s.operation.index == 1)
        f = grid.load_frame(own)
        original = load_wavefunctions(own.source_binding)
        close(f.basis.miller, original.miller, 0)
        assert f.gauge_identity.endswith(':raw-saved-columns')
    wrong_spin = replace(source, ispin=1)
    rejected(lambda: admit_source(ref, slot, wrong_spin))
    missing = replace(source, miller=source.miller[:, 1:], coefficients=source.coefficients[1:], igwx=source.igwx-1)
    rejected(lambda: admit_source(ref, slot, missing))
    # Same row count but an interior hole replaced by an outside row.
    corrupted = source.miller.copy(); corrupted[:, 0] = [7, 0, 0]
    rejected(lambda: admit_source(ref, slot, replace(source, miller=corrupted)))
    rejected(lambda: frame.basis.__class__(**{**frame.basis.__dict__, 'miller': frame.basis.miller+1}).validate(ref))
    # Original atom periodic images must be kept even though all translations vanish.
    inversion = next(s for s in grid.slots if np.any(s.atom_image_offsets))
    close(inversion.operation.rotation.T @ np.linalg.solve(ref.cell_bohr, ref.positions_bohr)
          - inversion.operation.fractional_translation[:,None],
          np.linalg.solve(ref.cell_bohr, ref.positions_bohr)[:, inversion.operation.equivalent_atoms-1]
          + inversion.atom_image_offsets)
    bad_op = replace(ref.saved_symmetries.operations[1], equivalent_atoms=np.array([2,1]))
    bad_saved = replace(ref.saved_symmetries, operations=(ref.saved_symmetries.operations[0], bad_op, *ref.saved_symmetries.operations[2:]))
    rejected(lambda: build_bloch_grid(replace(ref, saved_symmetries=bad_saved), bindings, symmetry_policy=UNITARY_POLICY))
    # Translated coordinate origin makes all nontrivial operation translations
    # nonzero but commensurate; compare both PW and source-FFT phases.
    shift = np.array([1,2,3])/9
    operations = tuple(replace(op, fractional_translation=op.rotation.T @ shift-shift)
                       for op in ref.saved_symmetries.operations)
    moved = replace(ref, positions_bohr=ref.positions_bohr+ref.cell_bohr @ shift[:,None],
                    saved_symmetries=replace(ref.saved_symmetries, operations=operations))
    moved_grid = build_bloch_grid(moved, bindings, symmetry_policy=UNITARY_POLICY)
    moved_slot = next(s for s in moved_grid.slots if s.operation.index > 1 and np.any(s.operation.fractional_translation))
    moved_source = load_wavefunctions(moved_slot.source_binding)
    moved_frame = moved_grid.load_frame(moved_slot)
    close(literal_wave(moved_frame.basis.miller, moved_frame.coefficients, x),
          source_fft_pullback(moved_source, moved_slot, ref.fft_smooth))
    sensitive(moved_frame.coefficients, moved_source.coefficients)
    # A raw nonsingular band-gauge change must transform D contragrediently.
    psi = frame.coefficients
    d = np.array([[.3,.04j],[-.04j,.2]])
    w = np.array([[1.1,.2j],[.13,.7]])
    inv = np.linalg.inv(w)
    changed = inv @ d @ inv.conj().T
    rho = psi @ d @ psi.conj().T
    close(rho, (psi @ w) @ changed @ (psi @ w).conj().T)
    sensitive(rho, (psi @ w) @ d @ (psi @ w).conj().T)
    # Iterate without accumulating frames; spin-down row order remains independent.
    seen = 0
    for streamed in iter_bloch_frames(grid):
        assert streamed.basis.miller.shape[1] == streamed.basis.slot.expected_rows
        seen += 1
        del streamed
    assert seen == 128
    return grid, frame, dict(literal_fft_residual=residual, full_slots=seen)


def literal_operator(m, k, positions, local_miller):
    """Independent finite PW local plus rank-one Gaussian NC matrix."""
    p = m.T + k
    out = np.zeros((len(p),len(p)), complex)
    support = {tuple(g) for g in local_miller.T}
    for i, row in enumerate(m.T):
        for j, col in enumerate(m.T):
            g = row-col
            if tuple(g) in support:
                out[i,j] += sum(-.7*np.exp(-.17*(g@g))*np.exp(-1j*g@tau) for tau in positions.T)
    for tau in positions.T:
        beta = np.exp(-.2*np.sum(p*p,axis=1))/np.sqrt(4*np.pi)*np.exp(-1j*(p@tau))
        out += 1.3*beta[:,None]*beta.conj()[None,:]
    return out


def action_oracle(inputs, frame):
    basis, ref = frame.basis, inputs.reference
    bridge = BlochQZeroBridge(inputs, basis)
    factors = build_pw_factors_for_basis(inputs, basis)
    action = BlochIonicAction(inputs, basis, factors=factors)
    m, k = basis.miller, basis.k_cart_bohr_inv
    n = m.shape[1]
    eye = np.eye(n, dtype=complex)  # deliberately small full-space oracle only
    direction = np.array([[.2+.13j,-.31,.17j],[-.19,.27j,.43]])
    prepared = action.prepare(direction)
    def full_apply(op):
        return np.column_stack([op.apply(eye[:, start:start+8]) for start in range(0,n,8)])
    actual = full_apply(prepared)
    eps = 2e-6
    expected = np.zeros_like(actual)
    for part, coefficient in ((direction.real,1), (direction.imag,1j)):
        plus = literal_operator(m,k,ref.positions_bohr+eps*part.T,inputs.density.miller)
        minus = literal_operator(m,k,ref.positions_bohr-eps*part.T,inputs.density.miller)
        expected += coefficient*(plus-minus)/(2*eps)
    residual = close(actual, expected, 2e-8)
    adjoint = np.column_stack([prepared.adjoint_apply(eye[:, start:start+8]) for start in range(0,n,8)])
    close(adjoint, actual.conj().T)
    close(full_apply(action.prepare(1j*direction)), 1j*actual)
    sensitive(actual, full_apply(action.prepare(direction.conj())))
    # Mutation: p=G changes the radial NC factors, not just an algebraically
    # cancelling constant in the two momentum-derivative sides.
    gamma_expected = np.zeros_like(actual)
    for part, coefficient in ((direction.real,1),(direction.imag,1j)):
        gamma_expected += coefficient*(literal_operator(m,np.zeros(3),ref.positions_bohr+eps*part.T,inputs.density.miller)
                          - literal_operator(m,np.zeros(3),ref.positions_bohr-eps*part.T,inputs.density.miller))/(2*eps)
    p_mutation = sensitive(actual, gamma_expected)
    local, bra, ket = prepared.components(eye[:, :8])
    sensitive(actual[:, :8], local+bra)
    sensitive(actual[:, :8], local+ket)
    x = coordinates(ref.fft_smooth)
    # coordinates() is C-order; explicitly convert to physical F-order fields.
    profile_c = np.cos(2*np.pi*x[:,0]) + .3j*np.sin(2*np.pi*x[:,1]) + .6
    field = profile_c.reshape(ref.fft_smooth).ravel(order='F')
    waves = np.exp(2j*np.pi*(x @ m))
    literal_local = waves.conj().T @ (profile_c[:,None]*waves)/len(x)
    close(bridge.local_apply(field, eye), literal_local)
    close(bridge.local_adjoint_apply(field, eye), literal_local.conj().T)
    left, right = eye[:, :3]+.2j*eye[:, -3:], eye[:, 2:5]+.17*eye[:, -4:-1]
    c, weight = np.array([.2,.4j,-.3]), .125
    density = bridge.density_from_dyads(left,right,c,weight=weight)
    literal_density_c = np.sum((waves @ left)*(waves @ right).conj()*c,axis=1)*weight/np.linalg.det(ref.cell_bohr)
    expected_density = literal_density_c.reshape(ref.fft_smooth).ravel(order='F')
    close(density[:, frame.ispin-1], expected_density)
    close(density[:, 2-frame.ispin], np.zeros(density.shape[0], dtype=complex))
    sensitive(density[:, frame.ispin-1], expected_density*np.linalg.det(ref.cell_bohr))
    sensitive(density[:, frame.ispin-1], density[:, 2-frame.ispin])
    pairing = ref.cell_bohr[0,0]**3/len(x)*np.sum(density[:,frame.ispin-1]*field)
    close(pairing, weight*np.sum(c*np.diag(right.conj().T @ literal_local @ left)))
    # Gradient comparator projects literal density onto the declared density G
    # support before differentiation, keeping orbital support out of this step.
    dg = np.exp(2j*np.pi*(x @ inputs.density.miller))
    dc = dg.conj().T @ literal_density_c/len(x)
    grad = bridge.density_gradient(density)
    for axis in range(3):
        expected_g = (dg @ (1j*inputs.density.miller[axis]*dc)).reshape(ref.fft_grid).ravel(order='F')
        close(grad[axis,:,frame.ispin-1], expected_g)
    sensitive(grad[:,:,frame.ispin-1], grad[:,:,frame.ispin-1]+1j*k[:,None]*expected_density)
    background = type('Background', (), {})()
    background.fft, background.inputs = bridge.density_fft, inputs
    background.vrs = np.column_stack([field.real, .7*field.real+.2])
    h0 = MaterialHamiltonian(background, factors, frame.ispin)
    p = m.T+k
    nonlocal_only = literal_operator(m,k,ref.positions_bohr,np.empty((3,0),int))
    v_c = background.vrs[:,frame.ispin-1].reshape(ref.fft_grid,order='F').ravel()
    h_expected = np.diag(np.sum(p*p,axis=1)) + waves.conj().T @ (v_c[:,None]*waves)/len(x) + nonlocal_only
    close(h0.apply(eye), h_expected)
    return dict(nonGamma_full_space_derivative_residual=residual, p_equals_G_mutation=p_mutation,
                full_pw_dimension=n, saved_columns=frame.coefficients.shape[1])


def main():
    with tempfile.TemporaryDirectory(prefix='a27-bloch-oracle-') as directory:
        inputs, bindings, root, rotations = fixture(directory)
        _, frame, grid_report = grid_oracle(inputs,bindings,root,rotations)
        actions = action_oracle(inputs,frame)
    print(json.dumps(dict(scope='NEW finite analytic model; not material execution',
                          grid=grid_report, actions=actions), indent=2))


if __name__ == '__main__':
    main()
