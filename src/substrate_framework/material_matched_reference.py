"""Complete conditional controls, evaluated at each ORIGINAL rounded geometry.

This is a matched frozen Taylor coefficient model, not a response ranking.
No covariance of moments, positive-time antiunitary evolution, SCF relaxation,
Native per-case Hessian, or equivalence to the older loaded remainder is assumed.
"""
import numpy as np

from .material_reference_state import prepare_source, sector_density_moments
from .material_frozen_state_curvature import (
    frozen_onsite_curvature, local_density_curvature, OnsiteCurvature)
from .material_physical_operators import build_physical_factors, physical_coupling
from .material_physical_field_moments import full_column_moments
from .material_mixed_field_ionic import full_column_mixed
from .qe_material_hamiltonian import MaterialHamiltonian
from .qe_real_bloch_grid import transport_real_state

MODEL = 'matched-direct-physical-onsite-plus-finite-ionic-taylor-v1'
KEYS = ('density_bohr3', 'positive_gradient_ry_bohr', 'local_onsite_ry_bohr2',
        'nonlocal_onsite_ry_bohr2', 'field_gradient', 'current_hartree',
        'contact', 'mixed_ionic_field', 'number')


def slot_identity(slot):
    return dict(spin=int(slot.ispin), mesh_index=list(slot.mesh_index),
                source_ik=int(slot.source_ik), source=dict(slot.source_binding),
                weight=float(slot.weight), operation=int(slot.operation.index),
                orbital_conjugation=bool(slot.orbital_conjugation),
                reciprocal_shift=slot.reciprocal_shift.tolist())


def ledger_slots(grid, ledger):
    """Authenticate original ragged population ledger, without re-filling."""
    n = len(grid.slots)
    fields = ('spin', 'mesh_index', 'source_ik', 'source_sha256', 'measures',
              'orbital_conjugation', 'reciprocal_shift', 'basis_identity', 'band_slices')
    if any(len(ledger[k]) != n for k in fields):
        raise ValueError('complete original slot ledger required')
    cursor, result = 0, []
    populations = np.asarray(ledger['populations'])
    for i, slot in enumerate(grid.slots):
        start, stop = map(int, ledger['band_slices'][i])
        if (start != cursor or stop-start != len(slot.eigenvalues_ry)
                or int(ledger['spin'][i]) != slot.ispin
                or tuple(ledger['mesh_index'][i]) != slot.mesh_index
                or int(ledger['source_ik'][i]) != slot.source_ik
                or str(ledger['source_sha256'][i]) != slot.source_binding['sha256']
                or float(ledger['measures'][i]) != slot.weight
                or bool(ledger['orbital_conjugation'][i]) != slot.orbital_conjugation
                or not np.array_equal(ledger['reciprocal_shift'][i], slot.reciprocal_shift)
                or not np.array_equal(ledger['energies_ry'][start:stop], slot.eigenvalues_ry)):
            raise ValueError('original slot/population/source measure mismatch')
        pop = populations[start:stop]
        if not np.isfinite(pop).all() or np.any((pop < 0) | (pop > 1)):
            raise ValueError('inadmissible original populations')
        result.append((i, slot, pop))
        cursor = stop
    if cursor != len(populations):
        raise ValueError('unaccounted population entries')
    return tuple(result)


def source_state(source, admission, selected, *, authority):
    if not selected:
        raise ValueError('nonempty original source orbit required')
    pop = selected[0][2]
    if any(not np.array_equal(pop, row[2]) for row in selected):
        raise ValueError('source orbit has different populations')
    return prepare_source(source.coefficients, pop, source_identity=source.source['sha256'],
                          basis_identity=admission['source_rows_sha256'], authority=authority)


def direct_slot(control, background, slot, source, state, admission, *, row_block=256,
                cached_field=None):
    """All rows/columns/atoms; build physical factors ONCE, inject into H and g.

    Real basis support is owned by the shared physical factory. Its own real
    scope validation must run; this module never casts it to a unitary basis.
    """
    frame, target = transport_real_state(control.grid, slot, source, state,
                                          source_admission=admission)
    factors = build_physical_factors(background.inputs, frame.basis)
    h = MaterialHamiltonian(background, factors, slot.ispin, 8)
    g = physical_coupling(background.inputs, frame.basis, factors, block_size=8)
    if cached_field is None:
        spatial = sector_density_moments(frame, target, g)
        curvature = frozen_onsite_curvature(frame.coefficients, target.density_coefficients, g,
                                            density_bohr3=spatial['density_bohr3'])
        optical = full_column_moments(h, frame.coefficients, target.density_coefficients,
                                      row_block=row_block)
        mixed = full_column_mixed(h, frame.coefficients, target.density_coefficients,
                                  row_block=row_block)
        field_reuse = None
    else:
        optical, mixed, field_reuse, legacy = cached_field
        spatial = dict(density_bohr3=legacy['density_bohr3'],
            number_pair=legacy['number_pair'],
            positive_gradient_ry_bohr=legacy['positive_gradient_ry_bohr'] +
                optical['delta_nonlocal_positive_gradient_ry_bohr'].reshape(-1))
        # The CACHE is A25: convert its nonlocal terms ONCE. Never add these
        # corrections to freshly calculated physical-factor coefficients.
        curvature = OnsiteCurvature(local_density_curvature(spatial['density_bohr3'],g),
            legacy['conditional_nonlocal_ry_bohr2']+optical['delta_nonlocal_onsite_ry_bohr2'])
    arrays = dict(density_bohr3=spatial['density_bohr3'],
        positive_gradient_ry_bohr=spatial['positive_gradient_ry_bohr'],
        local_onsite_ry_bohr2=curvature.local, nonlocal_onsite_ry_bohr2=curvature.nonlocal_term,
        field_gradient=optical['gradient'], current_hartree=optical['current_hartree'],
        contact=optical['contact'], mixed_ionic_field=mixed['mixed_ionic_field'],
        number=np.asarray(optical['number']), number_pair=spatial['number_pair'])
    if cached_field is None:
        arrays.update(
            physical_minus_factory_energy=np.asarray(optical['delta_nonlocal_energy_ry']),
            physical_minus_factory_gradient=optical['delta_nonlocal_positive_gradient_ry_bohr'],
            physical_minus_factory_onsite=optical['delta_nonlocal_onsite_ry_bohr2'])
    # Fresh deltas diagnose roundoff only; cached legacy conversion is above.
    metadata = dict(slot=slot_identity(slot), basis_identity=frame.basis.identity,
        frame_identity=frame.gauge_identity, rows=frame.coefficients.shape[0],
        columns=frame.coefficients.shape[1], optical_work=optical['work'], mixed_work=mixed['work'],
        scalar_channel_residuals=optical['rotational_scalar_residuals'],
        field_reuse=field_reuse,
        geometry_policy='direct-original-target', normalization='unweighted-one-explicit-spin-slot',
        orbital_conjugation_scope='initial Gram AND D only; not positive-time propagation')
    return arrays, metadata


def aggregate_complete(grid, records):
    """Deterministic original slot weights; rejects duplicates, omissions, extras.

    records yields (slot identity, unweighted arrays). No multiplicity shortcuts.
    """
    expected = {(s.ispin, tuple(s.mesh_index)): s for s in grid.slots}
    seen, totals, numbers = set(), None, np.zeros(2, complex)
    iterator = iter(records)
    # Caller supplies canonical order; enforce it so accumulation is repeatable.
    ordered = sorted(expected)
    for key in ordered:
        try:
            identity, arrays = next(iterator)
        except StopIteration as error:
            raise ValueError('incomplete full target coverage') from error
        actual = (identity['spin'], tuple(identity['mesh_index']))
        if actual in seen or actual != key or identity != slot_identity(expected[key]):
            raise ValueError('noncanonical, duplicated or mismatched original slot')
        seen.add(actual)
        if not all(k in arrays and np.isfinite(arrays[k]).all() for k in KEYS):
            raise ValueError('missing/nonfinite complete slot coefficients')
        weight = expected[key].weight
        if totals is None:
            totals = {k: np.zeros_like(arrays[k], dtype=complex) for k in KEYS if k != 'density_bohr3'}
            totals['density_bohr3'] = np.zeros((len(arrays['density_bohr3']), 2), complex)
        for name in KEYS:
            value = np.asarray(arrays[name])
            if name == 'density_bohr3':
                if value.shape != totals[name][:, key[0]-1].shape:
                    raise ValueError('physical density dimensions differ')
                totals[name][:, key[0]-1] += weight*value
            else:
                if value.shape != totals[name].shape:
                    raise ValueError('original coefficient dimensions differ: '+name)
                totals[name] += weight*value
        numbers[key[0]-1] += weight*arrays['number']
    if next(iterator, None) is not None:
        raise ValueError('extra target record')
    if totals is None:
        raise ValueError('empty material grid')
    totals['weighted_spin_number'] = numbers
    return totals


def assemble_frozen(aggregate, ionic_force, ionic_q):
    """Taylor gradient g_el-F_II; C(q)=C_II(q)+absolute onsite C_el."""
    onsite = OnsiteCurvature(aggregate['local_onsite_ry_bohr2'],
                            aggregate['nonlocal_onsite_ry_bohr2']).matrix()
    nat = len(aggregate['local_onsite_ry_bohr2'])
    ionic_q = np.asarray(ionic_q)
    if ionic_q.shape[-4:] != (nat, 3, nat, 3):
        raise ValueError('full all-atom ionic q grid required')
    force = np.asarray(ionic_force)
    if force.shape != (3*nat,):
        raise ValueError('all-atom bare ionic force required')
    from .material_force_lattice import q_grid_to_force_lattice
    q = ionic_q + onsite.reshape(nat, 3, nat, 3)
    return dict(electronic_onsite_matrix_ry_bohr2=onsite,
                ionic_force_ry_bohr=force,
                force_at_origin_ry_bohr=force-aggregate['positive_gradient_ry_bohr'],
                positive_total_gradient_ry_bohr=aggregate['positive_gradient_ry_bohr']-force,
                common_q_grid_ry_bohr2=q,
                common_force_lattice_ry_bohr2=q_grid_to_force_lattice(q))


def finite_ionic_reference(inputs, contract):
    """Sourced finite ionic linear/q-quadratic terms, NOT invented Native data.

    Linear scalar uses its explicit spherical G cutoff. The quadratic uses
    saved signed density G and d2ionq's shifted-G/head/rgen convention. Together
    they DEFINE the frozen Taylor polynomial; equality to one finite Ewald
    derivative or convergence is not asserted. Same prescription for every case.
    """
    from .qe_ionic_ewald import EwaldSumContract, ionic_linear_term
    from .qe_ionic_curvature import QEIonicCurvature
    ref = inputs.reference
    linear = ionic_linear_term(inputs, contract=EwaldSumContract(**contract),
                               coordinate_identity=ref.source['sha256']+':original-atom-order')
    curvature = QEIonicCurvature(cell_bohr=ref.cell_bohr, positions_bohr=ref.positions_bohr,
        valence_charges=linear.valence_charges,
        g_cart_bohr_inv=inputs.density.reciprocal_bohr_inv @ inputs.density.miller,
        alat_bohr=ref.alat_bohr, alpha_bohr_inv2=contract['alpha_bohr_inv2'],
        real_cutoff_bohr=contract['real_cutoff_bohr'], max_real_vectors=10000)
    shape = tuple(ref.monkhorst_pack['nk'+str(i)] for i in (1, 2, 3))
    if shape != (4, 4, 4):
        raise ValueError('matched full4x4x4 contract required')
    q = np.empty(shape+(ref.nat, 3, ref.nat, 3), complex)
    for index in np.ndindex(shape):
        fractional = np.asarray(index)/np.asarray(shape)
        fractional -= np.floor(fractional+0.5)
        q[index] = curvature.matrix(ref.reciprocal_bohr_inv @ fractional).reshape(ref.nat, 3, ref.nat, 3)
    masses = {s.name: s.mass_amu for s in ref.species}
    return linear, q, np.asarray([masses[n] for n in ref.atom_names])
