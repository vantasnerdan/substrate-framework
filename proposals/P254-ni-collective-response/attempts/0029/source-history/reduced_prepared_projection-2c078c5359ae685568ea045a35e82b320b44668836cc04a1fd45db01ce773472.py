"""Persisted per-slot preparation geometry and small NC overlaps.

Prepared edges require no native WFC, Gram eigensolve or full beta bank. A later
representative physical-action witness may independently rebuild full factors.
"""
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np
from substrate_framework.qe_real_bloch_grid import RealBlochBasis
from substrate_framework.qe_bloch_finiteq import BlochFiniteQBridge
from substrate_framework.qe_nc_radial import build_radial_tables
from reduced_paired_local_support import create_support
from reduced_local_gradient import local_gradient_matrix


@dataclass(frozen=True)
class SmallSpecies:
    name: str
    atom_indices: tuple
    active_d: np.ndarray


@dataclass(frozen=True)
class PreparedIonicEdge:
    inputs: object
    edge: object
    bridge: object
    nat: int
    ndof: int
    gq: np.ndarray
    vloc: dict
    species: tuple
    source_overlaps: dict
    target_overlaps: dict
    support_metadata: dict
    local_gradient: np.ndarray | None = None


def basis_arrays(basis):
    return {name: np.asarray(getattr(basis, name)) for name in
            ('miller', 'k_cart_bohr_inv', 'reciprocal_bohr_inv', 'row_map', 'source_miller')}


def restore_basis(arrays, metadata, slot, grid):
    basis = RealBlochBasis(**{name: arrays[name] for name in
        ('miller', 'k_cart_bohr_inv', 'reciprocal_bohr_inv', 'row_map', 'source_miller')},
        slot=slot, reference_source=MappingProxyType(dict(grid.reference.source)),
        source_admission=MappingProxyType(metadata['source_admission']),
        identity=metadata['basis_identity'], scope=grid.scope)
    return basis.validate(grid.reference)


def overlap_arrays(overlaps):
    out = {}
    for i, (name, (base, derivative)) in enumerate(sorted(overlaps.items())):
        out['overlap_base_%d' % i] = base
        out['overlap_derivative_%d' % i] = derivative
    return out


def restore_overlaps(arrays, names):
    return {name: (arrays['overlap_base_%d' % i], arrays['overlap_derivative_%d' % i])
            for i, name in enumerate(names)}


def prepare_ionic_edge(inputs, edge, source_overlaps, target_overlaps, minus_fractional,
                       *, local_cache=None):
    bridge = BlochFiniteQBridge(inputs, edge)
    cache_key = tuple(edge.q_fractional)
    cached = None if local_cache is None else local_cache.get(cache_key)
    cache_hit = cached is not None
    if cached is None:
        support = create_support(bridge.density_fft, edge.q_fractional, minus_fractional)
        gq = support.wavevectors.T
        moduli = np.linalg.norm(gq, axis=1)
        required = float(np.max(moduli, initial=0.))
        vloc = {}
        for name in dict.fromkeys(inputs.reference.atom_names):
            table = inputs.tables[name]
            if table.qmax_bohr_inv < required:
                table = build_radial_tables(table.upf, table.omega_bohr3, required)
            vloc[name] = table.vloc(moduli)
        gradient, gradient_metadata = local_gradient_matrix(inputs, support, vloc)
        cached = support, gq, vloc, gradient, gradient_metadata
        if local_cache is not None:
            local_cache[cache_key] = cached
            while len(local_cache) > 2:
                del local_cache[next(iter(local_cache))]
    support, gq, vloc, gradient, gradient_metadata = cached
    bridge.density_fft = support
    groups = []
    for name in dict.fromkeys(inputs.reference.atom_names):
        atoms = tuple(i+1 for i, n in enumerate(inputs.reference.atom_names) if n == name)
        species = SmallSpecies(name, atoms, inputs.tables[name].upf.expanded_d)
        groups.append((species, species))
    metadata = dict(original_count=support.original_count, paired_count=len(support.nl),
        aliased_additive_modes=int(len(support.nl)-len(np.unique(support.nl))),
        support_rule='original union (-original-q-qminus); original UPF radial extension',
        local_gradient_cache_hit=cache_hit, **gradient_metadata)
    return PreparedIonicEdge(inputs, edge, bridge, inputs.reference.nat, 3*inputs.reference.nat,
                            gq, vloc, tuple(groups), source_overlaps, target_overlaps, metadata, gradient)
