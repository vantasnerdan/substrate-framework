"""Join existing matched aggregate and preparations to ClampedArray, once only."""
import numpy as np

from matched_joint_case import load_case
from matched_clamped_isotopes import isotope_masses
from substrate_framework.material_clamped_array import ClampedArray
from substrate_framework.ionic_gaussian_preparation import prepare_local_trap_vacuum
from substrate_framework.material_joint_field_scalar import conditional_field_gaussian
from substrate_framework.physical_field_environment import UniformMaxwellField
from substrate_framework.material_modes import AMU_RY


def build_clamped(arrays, atomic_masses, *, authority, control, prepared=None,
                  reality_tolerance=1e-10):
    """No contractions: use full frozen electronic force/current from this case.

    independent-regions builds ONE local block, replicated by IndependentRegions.
    Intracell offdiagonal ionic couplings survive; only intercell C is removed.
    """
    if control not in ('periodic', 'independent-regions'):
        raise ValueError('explicit periodic or independent-regions control required')
    lattice = np.asarray(arrays['common_force_lattice_ry_bohr2'])
    shape = lattice.shape[:3]
    n = int(np.prod(shape))
    masses = np.repeat(np.asarray(atomic_masses)*AMU_RY, 3)
    charges = arrays['valence_charges']
    d = len(masses)
    ionic = prepared[0] if prepared else prepare_local_trap_vacuum(lattice, atomic_masses,
        state_authority=authority, reality_relative_tolerance=reality_tolerance)
    local = control == 'independent-regions'
    volume = float(abs(np.linalg.det(arrays['cell_bohr'])))
    field = prepared[1] if prepared and not local else UniformMaxwellField(
        volume*(1 if local else n), np.eye(3), authority+'; explicit '+control+' Maxwell volume')
    photon = prepared[2] if prepared and not local else conditional_field_gaussian(
        field, arrays['contact'], masses, charges, 1 if local else n,
        state_authority=authority, reality_relative_tolerance=reality_tolerance)
    force = (np.asarray(arrays['ionic_force_ry_bohr']).reshape(d)-
             np.asarray(arrays['positive_gradient_ry_bohr']).reshape(d))
    model = ClampedArray(lattice[:1,:1,:1] if local else lattice, masses, charges,
        np.zeros((d,d)), arrays['contact'], np.asarray(arrays['mixed_ionic_field']).reshape(d,3),
        force, arrays['field_gradient'], field.polarizations, field.mass,
        ionic.local.ordered, photon['ordered_covariance'], ionic.omega_ry)
    if local:
        from matched_clamped_regions import IndependentRegions
        model = IndependentRegions(model, shape)
    preparation = dict(ionic_ordered=ionic.local.ordered, trap_omega_ry=ionic.omega_ry,
        field_ordered=photon['ordered_covariance'], field_mass=field.mass,
        field_preparation=photon, force_per_cell_ry_bohr=force,
        field_gradient_per_cell=arrays['field_gradient'], mass_ry=masses)
    return model, preparation


def load_clamped(case_name, variant, control):
    case = load_case(case_name)
    masses, mass_record = isotope_masses(case.arrays, variant)
    prepared = ((case.ionic_preparation, case.field, case.field_preparation)
                if variant == 'original' else None)
    model, preparation = build_clamped(case.arrays, masses,
        authority=case.identity+'; '+variant+'; clamped-electron '+control,
        control=control, prepared=prepared, reality_tolerance=case.reality_tolerance)
    identity = dict(case=case_name, matched_joint_identity=case.identity,
        matched_joint_specification=case.specification, aggregate=case.aggregate_binding,
        archive=case.aggregate['archive'], mass=mass_record, control=control,
        atom_indices=case.arrays['atom_indices'], atom_names=case.arrays['atom_names'],
        positions_bohr=case.arrays['positions_bohr'], cell_bohr=case.arrays['cell_bohr'],
        valence_charges=case.charges, cell_indices=model.indices,
        preparation=preparation, phase_space_order='u(original atoms,xyz),p(original atoms,xyz)',
        field_order='a_xyz,Pi_xyz', electrons_held=True, extra_delta_curvature='exact zero')
    return model, identity
