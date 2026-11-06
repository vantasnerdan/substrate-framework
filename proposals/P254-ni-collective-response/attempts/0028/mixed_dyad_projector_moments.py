"""Physical field projector traces with original ordered-dyad accumulation.

Original-frame overlaps and placed radial jets are shared without identifying
independent terms. Coefficients remain separate from the original UPF channel
matrix. Each term adds species to its own kinetic field before its original
weight is applied, preserving cancellation-sensitive material arithmetic.
"""
import numpy as np

from mixed_physical_field_dyads import _overlap
from substrate_framework.material_mixed_field_ionic import first_cartesian_jets
from substrate_framework.material_physical_operators import (
    PhysicalFieldVertices, _weighted_overlap)
from substrate_framework.paired_field_moments import _place_first
from substrate_framework.qe_material_field_derivatives import UPFRadial


def _cached_overlap(cache, bank, frame, block_size):
    key = id(bank), id(frame)
    if key not in cache:
        cache[key] = _overlap(bank, frame, block_size)
    return cache[key]


def _with_coefficient(cache, species, overlap, coefficient):
    key = id(species.active_d), id(overlap), id(coefficient)
    if key not in cache:
        cache[key] = _weighted_overlap(species, overlap, overlap.shape[1])@coefficient
    return cache[key]


def field_projector_moments(vertex, terms, kinetic):
    """Return grouped Tr(L X), retaining both radial derivatives and weight.

    ``kinetic`` contains each original term's unweighted kinetic field trace.
    Only small original-column overlaps and one species/axis placed workspace
    are retained; no full-PW drives, density matrix or H0 action is formed.
    """
    if not isinstance(vertex, PhysicalFieldVertices):
        raise TypeError('original PhysicalFieldVertices required')
    terms, kinetic = tuple(terms), tuple(kinetic)
    if len(terms) != len(kinetic):
        raise ValueError('one kinetic field per original ordered term required')
    fields = [np.array(value, dtype=complex, copy=True) for value in kinetic]
    for species in vertex.hamiltonian.factors.species:
        inputs = vertex.hamiltonian.background.inputs
        jets = first_cartesian_jets(UPFRadial(inputs.tables[species.name]), vertex.momenta)
        derivative = np.empty_like(species.beta, order='F')
        base_overlaps, base_coefficients = {}, {}
        for axis, center in enumerate(jets[1:]):
            _place_first(vertex, species, center, derivative)
            # The placed bank is overwritten on the next axis.
            first_overlaps, first_coefficients = {}, {}
            for index, term in enumerate(terms):
                bl = _cached_overlap(base_overlaps, species.beta, term.left, vertex.block_size)
                br = _cached_overlap(base_overlaps, species.beta, term.right, vertex.block_size)
                dl = _cached_overlap(first_overlaps, derivative, term.left, vertex.block_size)
                dr = _cached_overlap(first_overlaps, derivative, term.right, vertex.block_size)
                weighted_base = _with_coefficient(base_coefficients, species, bl, term.D)
                weighted_first = _with_coefficient(first_coefficients, species, dl, term.D)
                fields[index][axis] += (np.vdot(dr, weighted_base)
                                        + np.vdot(br, weighted_first))
    result = {}
    for term, field in zip(terms, fields):
        if term.group not in result:
            result[term.group] = np.zeros(3, complex)
        result[term.group] += term.weight*field
    return result
