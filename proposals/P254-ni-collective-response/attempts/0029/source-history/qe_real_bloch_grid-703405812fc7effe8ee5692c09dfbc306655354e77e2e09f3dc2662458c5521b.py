"""Within-spin orbital conjugation for bound real collinear NC Hamiltonians.

This is NOT physical spin-1/2 time reversal.  It neither swaps spin nor
promotes lattice-only rotations to crystal operations.  See A28's bound
control-preparation-registration.json for original QE source blocks.
The old unitary BlochBasis/production authorities deliberately remain distinct.
"""
from dataclasses import dataclass
from hashlib import sha256
from itertools import product
import json
from types import MappingProxyType
import xml.etree.ElementTree as ET

import numpy as np

from .qe_bloch_grid import BlochSlot, _operations, qe_fold, readonly
from .qe_bloch_frame import admit_source, _digest
from .qe_material_grid import PhysicalRhoFFT
from .qe_saved_native import read_bound_bytes, load_wavefunctions

REAL_POLICY = 'qe7.5-real-collinear-orbital-conjugation-identity-unitary-first-v1'


def _identity(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


@dataclass(frozen=True)
class RealHamiltonianScope:
    reference_source: object
    upf_sources: tuple
    identity: str


def admit_real_hamiltonian(reference, upfs):
    """Check original XML plus decoded bound ordinary-NC UPFs, not a user flag.

    Reality here is the operator model's property, not a claim that every
    rounded saved background is exactly invariant under a crystal operation.
    """
    raw, source = read_bound_bytes(reference.source)
    root = ET.fromstring(raw)
    if (not reference.lsda or reference.noncolin or reference.spinorbit
            or reference.gamma_only or len(reference.nbnd_by_spin) != 2):
        raise ValueError('separate full-complex collinear spin Hamiltonians required')
    allowed = {'control_variables', 'atomic_species', 'atomic_structure', 'dft',
               'spin', 'bands', 'basis', 'electron_control', 'k_points_IBZ',
               'ion_control', 'cell_control', 'symmetry_flags', 'free_positions', 'twoch_'}
    inp = root.find('input')
    if inp is None or any(node.tag not in allowed for node in inp):
        raise ValueError('unadmitted XML field/model extension')
    for prefix in ('input', 'output'):
        dft = root.find(prefix + '/dft')
        if dft is None or [n.tag for n in dft] != ['functional'] or dft.findtext('functional', '').strip() != 'PBE':
            raise ValueError('only explicit pure PBE without Hubbard/hybrid extensions admitted')
    for path, expected in {
        'input/spin/lsda': 'true', 'input/spin/noncolin': 'false',
        'input/spin/spinorbit': 'false', 'output/algorithmic_info/uspp': 'false',
        'output/algorithmic_info/paw': 'false', 'input/symmetry_flags/noinv': 'false',
        'input/symmetry_flags/no_t_rev': 'false', 'input/control_variables/fcp': 'false',
        'input/control_variables/rism': 'false', 'input/twoch_/twochem': 'false',
    }.items():
        if root.findtext(path, '').strip() != expected:
            raise ValueError('unadmitted real-Hamiltonian flag: ' + path)
    spin = root.find('input/spin')
    if {n.tag for n in spin} != {'lsda', 'noncolin', 'spinorbit'}:
        raise ValueError('external/constrained spin extensions not admitted')
    if set(upfs) != {s.name for s in reference.species}:
        raise ValueError('exact original species UPF coverage required')
    sources = []
    for species in reference.species:
        upf = upfs[species.name]
        # Re-read header from the already bound bytes retained by NCUPF.
        if sha256(upf.original_bytes).hexdigest() != upf.sha256 or upf.element != species.name:
            raise ValueError('UPF original bytes/species differ')
        h = ET.fromstring(upf.original_bytes).find('PP_HEADER').attrib
        if (h.get('pseudo_type') != 'NC' or h.get('relativistic') != 'scalar'
                or h.get('functional', '').strip() != 'PBE'
                or any(h.get(k, '').lower() not in ('f', 'false', '.false.')
                       for k in ('has_so', 'is_ultrasoft', 'is_paw', 'is_coulomb', 'core_correction'))
                or not np.array_equal(upf.raw_dij, upf.raw_dij.T)
                or any(np.iscomplexobj(a) or not np.isfinite(a).all()
                       for a in (upf.raw_dij, upf.beta_r, upf.local_ry))):
            raise ValueError('real scalar ordinary-NC UPF required')
        sources.append(dict(path=upf.path, bytes=str(len(upf.original_bytes)), sha256=upf.sha256))
    identity = _identity(dict(policy=REAL_POLICY, xml=dict(source), upfs=sources))
    return RealHamiltonianScope(MappingProxyType(dict(source)),
                                tuple(MappingProxyType(s) for s in sources), identity)


@dataclass(frozen=True)
class RealBlochSlot(BlochSlot):
    orbital_conjugation: bool
    scope_identity: str

    @property
    def epsilon(self):
        return -1 if self.orbital_conjugation else 1


@dataclass(frozen=True)
class RealBlochGrid:
    reference: object
    slots: tuple
    scope: RealHamiltonianScope
    symmetry_policy: str = REAL_POLICY

    def require_slot(self, slot):
        if not any(slot is s for s in self.slots):
            raise ValueError('slot does not belong to this real-Hamiltonian grid')

    def iter_slots(self):
        return iter(self.slots)

    def load_frame(self, slot):
        self.require_slot(slot)
        return transport_real_frame(self, slot, load_wavefunctions(slot.source_binding))


def validate_real_coverage(reference, slots):
    """Independent cardinality/source/weight gate; no metadata-only admission claim."""
    mesh = tuple(reference.monkhorst_pack['nk' + str(i)] for i in (1, 2, 3))
    count = int(np.prod(mesh))
    expected = {(spin, index) for spin in (1, 2) for index in product(*(range(n) for n in mesh))}
    if len(slots) != len(expected) or {(s.ispin, s.mesh_index) for s in slots} != expected:
        raise ValueError('incomplete or duplicate explicit spin/mesh coverage')
    for spin in (1, 2):
        for point in reference.kpoints:
            orbit = [s for s in slots if s.ispin == spin and s.source_ik == point.ik]
            if (not orbit or any(s.weight != 1 / count for s in orbit)
                    or not np.isclose(len(orbit) / count, point.weight, atol=2e-12, rtol=0)):
                raise ValueError('source orbit/weight differs from original XML')
    if {s.source_ik for s in slots} != {p.ik for p in reference.kpoints}:
        raise ValueError('extra or missing original source')


def build_real_bloch_grid(reference, wfc_bindings, *, scope, symmetry_policy):
    if (symmetry_policy != REAL_POLICY or not isinstance(scope, RealHamiltonianScope)
            or dict(scope.reference_source) != dict(reference.source)):
        raise ValueError('explicit bound orbital-conjugation policy/scope required')
    mesh = tuple(reference.monkhorst_pack['nk' + str(i)] for i in (1, 2, 3))
    shifts = np.array([reference.monkhorst_pack['k' + str(i)] for i in (1, 2, 3)])
    if min(mesh) < 1 or np.any(shifts):
        raise ValueError('complete unshifted mesh required by this preparation policy')
    keys = {(spin, p.ik) for spin in (1, 2) for p in reference.kpoints}
    if set(wfc_bindings) != keys:
        raise ValueError('exact complete explicit spin/source binding coverage required')
    bindings = {key: MappingProxyType({f: str(value[f]) for f in ('path', 'bytes', 'sha256')})
                for key, value in wfc_bindings.items()}
    operations = _operations(reference)  # admitted crystal operations ONLY
    sources = [(p, np.linalg.solve(reference.reciprocal_bohr_inv, p.k_cart_bohr_inv))
               for p in reference.kpoints]
    slots = []
    for number, index in enumerate(product(*(range(n) for n in mesh))):
        k = qe_fold(np.array(index) / mesh)
        choices = []
        for conjugate in (False, True):
            epsilon = -1 if conjugate else 1
            for op, images in operations:
                identity = np.array_equal(op.rotation, np.eye(3)) and not np.any(op.fractional_translation)
                for point, ks in sources:
                    delta = k - epsilon * op.rotation @ ks
                    g = np.rint(delta).astype(np.int64)
                    if np.allclose(delta, g, atol=2e-10, rtol=0):
                        preferred = not conjugate and identity and np.allclose(k, ks, atol=2e-10, rtol=0)
                        choices.append((not preferred, conjugate, op.index, point.ik, point, ks, op, images, g))
        if not choices:
            raise ValueError('real crystal group does not cover mesh slot ' + str(index))
        _, conjugate, _, _, point, ks, op, images, g = min(choices, key=lambda c: c[:4])
        for spin in (1, 2):
            _, energies, raw = reference.select(point.ik, spin)
            slots.append(RealBlochSlot(spin, index, number, point.ik, bindings[spin, point.ik], op,
                readonly(g), readonly(ks), readonly(k), readonly(reference.reciprocal_bohr_inv @ k),
                energies, raw, point.npw, 1 / int(np.prod(mesh)), REAL_POLICY, images,
                conjugate, scope.identity))
    validate_real_coverage(reference, slots)
    return RealBlochGrid(reference, tuple(slots), scope)


def real_basis_identity(miller, reciprocal, slot):
    return _identity(dict(arrays=_digest(miller, reciprocal, slot.k_cart_bohr_inv,
        slot.operation.rotation, slot.operation.fractional_translation, slot.reciprocal_shift),
        source=dict(slot.source_binding), policy=slot.gauge_policy, scope=slot.scope_identity,
        spin=slot.ispin, mesh=slot.mesh_index, ik=slot.source_ik,
        operation=slot.operation.index, orbital_conjugation=slot.orbital_conjugation))


@dataclass(frozen=True)
class RealBlochBasis:
    miller: np.ndarray
    k_cart_bohr_inv: np.ndarray
    reciprocal_bohr_inv: np.ndarray
    row_map: np.ndarray
    source_miller: np.ndarray
    slot: RealBlochSlot
    reference_source: object
    source_admission: object
    identity: str
    scope: RealHamiltonianScope

    def validate(self, reference):
        s, audit = self.slot, self.source_admission
        if (dict(self.reference_source) != dict(reference.source)
                or dict(self.scope.reference_source) != dict(reference.source)
                or s.scope_identity != self.scope.identity or s.gauge_policy != REAL_POLICY
                or type(s.orbital_conjugation) is not bool):
            raise ValueError('different real-Hamiltonian basis authority')
        point, energies, raw = reference.select(s.source_ik, s.ispin)
        operations = {op.index: (op, image) for op, image in _operations(reference)}
        if s.operation.index not in operations:
            raise ValueError('operation is not an admitted crystal symmetry')
        op, images = operations[s.operation.index]
        if (any(not np.array_equal(getattr(op, name), getattr(s.operation, name)) for name in
                ('rotation', 'fractional_translation', 'equivalent_atoms'))
                or not np.array_equal(images, s.atom_image_offsets)
                or audit['source_sha256'] != s.source_binding['sha256']
                or audit['ispin'] != s.ispin or audit['ik'] != s.source_ik
                or audit['source_rows_sha256'] != _digest(self.source_miller)
                or self.source_miller.shape != (3, point.npw) or s.expected_rows != point.npw
                or not np.array_equal(s.eigenvalues_ry, energies) or not np.array_equal(s.raw_occupations, raw)):
            raise ValueError('source rows/spin/labels/operation admission mismatch')
        mesh = np.array([reference.monkhorst_pack['nk' + str(i)] for i in (1, 2, 3)])
        if (not np.allclose(reference.reciprocal_bohr_inv @ s.source_k_fractional,
                            point.k_cart_bohr_inv, atol=2e-13, rtol=0)
                or not np.allclose(s.k_fractional, qe_fold(np.array(s.mesh_index) / mesh), atol=2e-10, rtol=0)
                or not np.allclose(s.k_fractional, s.epsilon * op.rotation @ s.source_k_fractional
                                   + s.reciprocal_shift, atol=2e-10, rtol=0)
                or not np.array_equal(s.reciprocal_shift, np.rint(s.reciprocal_shift))
                or not np.array_equal(self.miller, s.epsilon * op.rotation @ self.source_miller
                                      - s.reciprocal_shift[:, None])
                or not np.array_equal(self.row_map, np.arange(point.npw))
                or not np.array_equal(self.k_cart_bohr_inv, s.k_cart_bohr_inv)
                or not np.allclose(self.k_cart_bohr_inv, reference.reciprocal_bohr_inv @ s.k_fractional, atol=2e-13, rtol=0)
                or not np.allclose(self.reciprocal_bohr_inv, reference.reciprocal_bohr_inv, atol=2e-13, rtol=0)
                or self.identity != real_basis_identity(self.miller, self.reciprocal_bohr_inv, s)):
            raise ValueError('real orbital row/fold/phase identity differs')
        PhysicalRhoFFT.create(reference.fft_smooth, self.miller, self.reciprocal_bohr_inv,
                              abs(np.linalg.det(reference.cell_bohr)))
        return self


def real_basis(grid, slot, source, admission):
    grid.require_slot(slot)
    if (dict(source.source) != dict(slot.source_binding) or source.ispin != slot.ispin
            or source.ik != slot.source_ik or admission['source_rows_sha256'] != _digest(source.miller)):
        raise ValueError('cached admission does not belong to actual source')
    m = readonly(slot.epsilon * slot.operation.rotation @ source.miller - slot.reciprocal_shift[:, None])
    return RealBlochBasis(m, slot.k_cart_bohr_inv, source.reciprocal_bohr_inv,
        readonly(np.arange(source.igwx)), source.miller, slot, grid.reference.source,
        admission, real_basis_identity(m, source.reciprocal_bohr_inv, slot), grid.scope).validate(grid.reference)


def transport_coefficients(slot, miller, coefficients):
    """Inherited-order c_t; exact -k/-G fold and conjugated translation phase."""
    result = coefficients.conj() if slot.orbital_conjugation else coefficients
    f = slot.operation.fractional_translation
    if np.any(f):
        phase = np.exp(-slot.epsilon * 2j * np.pi * ((miller.T + slot.source_k_fractional) @ f))
        result = phase[:, None] * result
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class RealBlochFrame:
    basis: RealBlochBasis
    coefficients: np.ndarray
    eigenvalues_ry: np.ndarray
    raw_occupations: np.ndarray
    gauge_identity: str

    @property
    def ispin(self):
        return self.basis.slot.ispin


def transport_real_frame(grid, slot, source, *, source_admission=None):
    grid.require_slot(slot)
    audit = admit_source(grid.reference, slot, source) if source_admission is None else source_admission
    basis = real_basis(grid, slot, source, audit)
    c = transport_coefficients(slot, source.miller, source.coefficients)
    return RealBlochFrame(basis, c, slot.eigenvalues_ry, slot.raw_occupations,
                          basis.identity + ':raw-saved-columns')


def transport_real_state(grid, slot, source, state, *, source_admission=None):
    from .material_reference_transport import raw_state_transport
    if state.source_identity != source.source['sha256'] or state.gram.shape != (source.nbnd, source.nbnd):
        raise ValueError('state is not the original full-column source state')
    frame = transport_real_frame(grid, slot, source, source_admission=source_admission)
    return frame, raw_state_transport(state, antiunitary=slot.orbital_conjugation)


def build_real_pw_factors(inputs, basis):
    """Explicit new policy adapter; NEVER pretend this is a unitary BlochBasis."""
    from .qe_material_hamiltonian import _build_basis_factors
    if not isinstance(basis, RealBlochBasis):
        raise ValueError('checked real-Hamiltonian Bloch basis required')
    basis.validate(inputs.reference)
    actual = {(name, table.upf.sha256) for name, table in inputs.tables.items()}
    expected = {(s.name, next(x['sha256'] for x in basis.scope.upf_sources
                 if x['path'].endswith('/' + s.pseudo_file))) for s in inputs.reference.species}
    if actual != expected:
        raise ValueError('factor UPFs differ from the admitted real-Hamiltonian scope')
    s = basis.slot
    provenance = dict(kind='real-collinear-orbital-conjugation-basis', identity=basis.identity,
        source=dict(s.source_binding), xml_source=dict(basis.reference_source), scope=basis.scope.identity,
        ispin=s.ispin, mesh_index=s.mesh_index, source_ik=s.source_ik, operation=s.operation.index,
        orbital_conjugation=s.orbital_conjugation, gauge_policy=REAL_POLICY)
    return _build_basis_factors(inputs, basis.miller, basis.k_cart_bohr_inv,
                               basis.reciprocal_bohr_inv, None, provenance)


def transform_real_vector(grid, slot, value, *, kind):
    """Covariance at transformed geometry, not saved-background invariance.

    Full atom-major force/positive-gradient: +S and irt. Integrated physical
    current: epsilon*S, including NC nonlocal velocity. Complex raw roundoff
    is conjugated for the antilinear route, never discarded or repaired.
    """
    grid.require_slot(slot)
    value = np.asarray(value)
    transformed = value.conj() if slot.orbital_conjugation else value
    s = slot.operation.rotation
    if kind in ('force', 'positive-gradient'):
        if value.shape != (grid.reference.nat, 3):
            raise ValueError('all original atom-major Cartesian vectors required')
        return transformed[slot.operation.equivalent_atoms - 1] @ s.T
    if kind == 'current':
        if value.shape != (3,):
            raise ValueError('integrated Cartesian current vector required')
        return slot.epsilon * s @ transformed
    raise ValueError('explicit force/positive-gradient/current parity required')
