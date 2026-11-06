"""Explicit spin/mesh slots for the XML-order UNITARY canonical extension.

Metadata construction does not decode WFCs or assert a physical preparation.
Identity is preferred on saved representatives; otherwise the first XML
operation wins. This policy is not a little-group invariance theorem.
"""
from dataclasses import dataclass
from itertools import product
from types import MappingProxyType
import numpy as np

UNITARY_POLICY = 'qe7.5-xml-order-unitary-identity-first-v1'


def readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def qe_fold(value):
    """QE NINT ties away from zero; positive mesh half-boundary becomes -1/2."""
    value = np.asarray(value, dtype=float)
    return value - np.copysign(np.floor(np.abs(value) + .5), value)


@dataclass(frozen=True)
class BlochSlot:
    ispin: int
    mesh_index: tuple[int, int, int]
    mesh_slot: int                       # zero based, h3 fastest; NOT source ik
    source_ik: int
    source_binding: object
    operation: object
    reciprocal_shift: np.ndarray         # kt = S ks + g
    source_k_fractional: np.ndarray
    k_fractional: np.ndarray
    k_cart_bohr_inv: np.ndarray
    eigenvalues_ry: np.ndarray
    raw_occupations: np.ndarray
    expected_rows: int
    weight: float                        # full-grid per-cell mean quadrature
    gauge_policy: str
    atom_image_offsets: np.ndarray       # S.T x_a - f = x_irt(a) + integer


@dataclass(frozen=True)
class BlochGrid:
    reference: object
    slots: tuple[BlochSlot, ...]
    symmetry_policy: str

    def iter_slots(self):
        return iter(self.slots)

    def load_frame(self, slot):
        """Decode one bound source and return its checked transported frame."""
        from .qe_bloch_frame import load_bloch_frame
        self.require_slot(slot)
        return load_bloch_frame(self.reference, slot)

    def require_slot(self, slot):
        if not any(slot is existing for existing in self.slots):
            raise ValueError('slot does not belong to this grid')


def _operations(reference):
    saved = reference.saved_symmetries
    if saved is None or reference.noncolin or reference.spinorbit or reference.gamma_only:
        raise ValueError('saved scalar full-complex crystal symmetries required')
    # Deliberately actual-geometry scope; no unsourced general-cell transpose rule.
    cell = reference.cell_bohr
    if (not np.allclose(cell, np.eye(3)*cell[0, 0], atol=2e-12, rtol=0)
            or cell[0, 0] <= 0 or len(set(reference.fft_smooth)) != 1):
        raise ValueError('unitary route requires diagonal cubic cell/equal FFT dimensions')
    x = np.linalg.solve(cell, reference.positions_bohr)
    accepted = []
    for op in saved.operations[:saved.nsym]:
        s, f, irt = op.rotation, op.fractional_translation, op.equivalent_atoms
        if op.effective_time_reversal is not False:
            continue
        if f is None or irt is None or not np.array_equal(s.T @ s, np.eye(3, dtype=int)):
            raise ValueError('complete signed-permutation crystal operation required')
        if not np.allclose(f * reference.fft_smooth, np.rint(f * reference.fft_smooth), atol=2e-10, rtol=0):
            raise ValueError('translation is not commensurate with physical FFT grid')
        delta = s.T @ x - f[:, None] - x[:, irt-1]
        images = np.rint(delta).astype(np.int64)
        if not np.allclose(delta, images, atol=2e-7, rtol=0):
            raise ValueError('saved atom permutation/pullback disagrees with actual geometry')
        if any(reference.atom_names[a] != reference.atom_names[b-1] for a, b in enumerate(irt)):
            raise ValueError('unitary operation changes chemical species')
        accepted.append((op, readonly(images)))
    if not accepted:
        raise ValueError('no admitted unitary crystal operations')
    return accepted


def build_bloch_grid(reference, wfc_bindings, *, symmetry_policy):
    """Build complete metadata; bindings map (one-based ispin, source ik) to bytes.

    Every spin header/row set is admitted independently on load, never inferred
    from another spin. No all-frame cache or invented XML k record is created.
    """
    if symmetry_policy != UNITARY_POLICY:
        raise ValueError('explicit XML-order unitary policy required')
    mesh = tuple(reference.monkhorst_pack[f'nk{i}'] for i in (1, 2, 3))
    shifts = np.array([reference.monkhorst_pack[f'k{i}'] for i in (1, 2, 3)])
    if min(mesh) < 1 or not np.isin(shifts, (0, 1)).all():
        raise ValueError('invalid Monkhorst-Pack mesh')
    keys = {(spin, p.ik) for spin in range(1, len(reference.nbnd_by_spin)+1) for p in reference.kpoints}
    if set(wfc_bindings) != keys:
        raise ValueError('exact complete explicit spin/source binding coverage required')
    bindings = {key: MappingProxyType({field: str(value[field]) for field in ('path', 'bytes', 'sha256')})
                for key, value in wfc_bindings.items()}
    operations = _operations(reference)
    sources = [(p, np.linalg.solve(reference.reciprocal_bohr_inv, p.k_cart_bohr_inv)) for p in reference.kpoints]
    slots, counts = [], {p.ik: 0 for p in reference.kpoints}
    for number, index in enumerate(product(*(range(n) for n in mesh))):
        k = qe_fold((np.array(index) + shifts/2) / mesh)
        choices = []
        for op, images in operations:
            identity = np.array_equal(op.rotation, np.eye(3, dtype=int)) and np.array_equal(op.fractional_translation, np.zeros(3))
            for point, ks in sources:
                delta = k - op.rotation @ ks
                g = np.rint(delta).astype(np.int64)
                if np.allclose(delta, g, atol=2e-10, rtol=0):
                    preferred = identity and np.allclose(k, ks, atol=2e-10, rtol=0)
                    choices.append((not preferred, op.index, point.ik, point, ks, op, images, g))
        if not choices:
            raise ValueError(f'unitary operations do not cover full mesh slot {index}')
        _, _, _, point, ks, op, images, g = min(choices, key=lambda c: c[:3])
        counts[point.ik] += 1
        for spin in range(1, len(reference.nbnd_by_spin)+1):
            _, energies, raw = reference.select(point.ik, spin)
            slots.append(BlochSlot(spin, index, number, point.ik, bindings[spin, point.ik], op,
                         readonly(g), readonly(ks), readonly(k), readonly(reference.reciprocal_bohr_inv @ k),
                         energies, raw, point.npw, 1/np.prod(mesh), symmetry_policy, images))
    for point in reference.kpoints:
        if not np.isclose(counts[point.ik]/np.prod(mesh), point.weight, atol=2e-12, rtol=0):
            raise ValueError('unitary orbit counts disagree with original XML quadrature')
    return BlochGrid(reference, tuple(slots), symmetry_policy)


def iter_bloch_frames(grid, requested_slots=None):
    """Yield source-grouped frames, retaining only one source and one target.

    Consumer must release each yielded frame before advancing to retain this
    memory bound. Retaining outputs is explicitly the caller's allocation.
    """
    from .qe_saved_native import load_wavefunctions
    from .qe_bloch_frame import admit_source, transport_bloch_frame
    slots = list(grid.slots if requested_slots is None else requested_slots)
    for slot in slots:
        grid.require_slot(slot)
    slots.sort(key=lambda s: (s.ispin, s.source_ik, s.mesh_slot))
    source, key = None, None
    for slot in slots:
        next_key = (slot.ispin, slot.source_ik)
        if next_key != key:
            source = None
            source = load_wavefunctions(slot.source_binding)
            key = next_key
            audit = admit_source(grid.reference, slot, source)
        yield transport_bloch_frame(grid.reference, slot, source, source_admission=audit)
