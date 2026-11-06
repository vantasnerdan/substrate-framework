"""One named two-band folded nonGamma q=0 consumer of accepted A25 fields."""
from dataclasses import asdict
import json
import numpy as np
from .qe_saved_native import load_wavefunctions
from .qe_bloch_frame import transport_bloch_frame
from .qe_material_hamiltonian import build_pw_factors_for_basis, MaterialHamiltonian
from .qe_bloch_ionic_action import BlochIonicAction
from .qe_bloch_qzero import BlochQZeroBridge


def selected_slot(grid):
    return next(slot for slot in grid.slots if slot.ispin == 1 and np.any(slot.reciprocal_shift)
                and np.any(slot.k_cart_bohr_inv))


def probe_arrays(background, grid, descriptor, budget):
    inputs, slot = background.inputs, selected_slot(grid)
    source = load_wavefunctions(slot.source_binding)
    frame = transport_bloch_frame(inputs.reference,slot,source)
    if frame.basis.identity != descriptor.basis_identity:
        raise ValueError('selected raw frame differs from completed coverage')
    raw = frame.coefficients[:, :2]
    budget()
    factors = build_pw_factors_for_basis(inputs,frame.basis)
    h0 = MaterialHamiltonian(background,factors,slot.ispin,block_size=2).apply(raw)
    budget()
    ionic = BlochIonicAction(inputs,frame.basis,factors=factors)
    direction = np.zeros((37,3),complex)
    direction[31,0], direction[0,1] = 1., 1j
    local, bra, ket = ionic.prepare(direction).components(raw)
    budget()
    bridge = ionic.bridge
    shape = inputs.reference.fft_grid
    # Physical F-order grid, exact named cos(b1)+i sin(b2) profile.
    x = np.arange(shape[0])/shape[0]
    y = np.arange(shape[1])/shape[1]
    profile = np.broadcast_to(np.cos(2*np.pi*x)[:,None,None]
                 + 1j*np.sin(2*np.pi*y)[None,:,None],shape).ravel(order='F')
    number = bridge.local_apply(profile,raw)
    coefficient, weight = np.array([1+.25j]), slot.weight
    density = bridge.density_from_dyads(raw[:, :1],raw[:,1:2],coefficient,weight=weight)
    gradient = bridge.density_gradient(density)
    omega, nnr = bridge.orbital_fft.omega, bridge.orbital_fft.nnr
    pairs = np.array([[omega/nnr*np.sum(density[:,0]), weight*coefficient[0]*np.vdot(raw[:,1],raw[:,0])],
                      [omega/nnr*np.sum(density[:,0]*profile),weight*coefficient[0]*np.vdot(raw[:,1],number[:,0])]])
    budget()
    return dict(probe_slot_json=np.asarray(json.dumps(asdict(descriptor),sort_keys=True)),
                bands_one_based=np.array([1,2],dtype=np.int64), source_block=source.coefficients[:,:2],
                raw_block=raw, target_miller=frame.basis.miller, h0_action_ry=h0,
                saved_label_residual_ry=h0-raw*frame.eigenvalues_ry[None,:2],
                ionic_local=local,ionic_bra=bra,ionic_ket=ket,ionic_action=local+bra+ket,
                direction=direction,number_profile=profile,number_action=number,
                density=density,density_gradient=gradient,dyad_coefficients=coefficient,
                mode_weight=np.asarray(weight),comparison_pairs=pairs)


def validate_probe_arrays(arrays, grid, records):
    slot = selected_slot(grid)
    descriptor = next(r for r in records if (r.ispin,r.mesh_slot) == (slot.ispin,slot.mesh_slot))
    npw, nnr = descriptor.n_pw, int(np.prod(grid.reference.fft_grid))
    shapes = dict(probe_slot_json=(),bands_one_based=(2,),source_block=(npw,2),raw_block=(npw,2),
                  target_miller=(3,npw),h0_action_ry=(npw,2),saved_label_residual_ry=(npw,2),
                  ionic_local=(npw,2),ionic_bra=(npw,2),ionic_ket=(npw,2),ionic_action=(npw,2),
                  direction=(37,3),number_profile=(nnr,),number_action=(npw,2),density=(nnr,2),
                  density_gradient=(3,nnr,2),dyad_coefficients=(1,),mode_weight=(),comparison_pairs=(2,2))
    if set(arrays) != set(shapes):
        raise ValueError('complete selected probe archive required')
    for key,shape in shapes.items():
        a = arrays[key]
        if a.shape != shape or (a.dtype.kind in 'fc' and not np.isfinite(a).all()):
            raise ValueError('invalid selected probe array: '+key)
    if json.loads(str(arrays['probe_slot_json'])) != json.loads(json.dumps(asdict(descriptor))):
        raise ValueError('selected probe does not belong to complete coverage')
    if (not np.array_equal(arrays['bands_one_based'],[1,2]) or arrays['mode_weight'] != slot.weight
            or not np.array_equal(arrays['dyad_coefficients'],[1+.25j])):
        raise ValueError('selected raw bands/dyad/weight differ')
    from material_pilot_io import array_record
    if array_record(arrays['target_miller']) != descriptor.target_miller_record:
        raise ValueError('selected target PW row identity differs')
    # All actual saved translations are zero; retain explicit transport formula.
    source_miller = np.linalg.solve(slot.operation.rotation,
                       arrays['target_miller']+slot.reciprocal_shift[:,None])
    phase = np.exp(-2j*np.pi*((source_miller.T+slot.source_k_fractional) @ slot.operation.fractional_translation))
    expected = phase[:,None]*arrays['source_block']
    if not np.array_equal(expected,arrays['raw_block']):
        raise ValueError('selected raw columns were changed beyond declared unitary transport')
    for name in ('integrated_dyad_trace','number_density_bilinear'):
        index = ('integrated_dyad_trace','number_density_bilinear').index(name)
        left,right = arrays['comparison_pairs'][index]
        if abs(left-right) > 1e-12+1e-10*abs(right):
            raise ValueError('selected bridge comparison failed: '+name)
    if not np.array_equal(arrays['ionic_action'],arrays['ionic_local']+arrays['ionic_bra']+arrays['ionic_ket']):
        raise ValueError('both ionic derivative sides are not retained')
    return {name: {'left':[float(pair[0].real),float(pair[0].imag)],
                   'right':[float(pair[1].real),float(pair[1].imag)],
                   'absolute_error':float(abs(pair[0]-pair[1]))}
            for name,pair in zip(('integrated_dyad_trace','number_density_bilinear'),arrays['comparison_pairs'])}
