"""Restore a completed authenticated population specification without refilling."""
from itertools import groupby
from math import fsum

from .material_preparation import FullGridSector, MaterialPopulationSpecification, SourceBinding
from .material_preparation_inputs import preparation_policy
from .zero_temperature_preparation import (ModeLevel, FillingPolicy, DegenerateBlock,
                                           ZeroTemperatureFilling)


def restore_population_specification(coverage, arrays, metadata):
    """Caller authenticates receipt and validates arrays against coverage first."""
    reference = SourceBinding.from_mapping(coverage.grid.reference.source)
    if SourceBinding.from_mapping(metadata['reference_binding']) != reference:
        raise ValueError('saved population XML source differs')
    charge_bindings = tuple((name, SourceBinding.from_mapping(record))
                           for name, record in metadata['charge_bindings'])
    expected = tuple((name, SourceBinding.from_mapping(record))
                     for name, record in coverage.recipe['inputs']['upfs'].items())
    if dict(charge_bindings) != dict(expected):
        raise ValueError('saved population valence source differs from complete coverage')
    sectors = tuple(FullGridSector(r.ispin, r.mesh_index, r.source_ik,
        SourceBinding.from_mapping(r.source_binding), r.gauge_identity, r.basis_identity,
        r.n_pw, coverage.coverage_identity + ':' + r.gauge_identity) for r in coverage.slot_records)
    levels = tuple(ModeLevel(f'{sector.key}:band={band + 1}', sector.spin,
                    float(arrays['energies_ry'][i, band]), float(arrays['measures'][i]))
                   for i, sector in enumerate(sectors) for band in range(346))
    populations = tuple(float(p) for p in arrays['populations'].ravel())
    blocks = []
    for energy, group in groupby(sorted(range(len(levels)), key=lambda i: levels[i].energy),
                                 key=lambda i: levels[i].energy):
        indices = tuple(group)
        blocks.append(DegenerateBlock(energy, indices, fsum(levels[i].measure for i in indices),
                                      populations[indices[0]]))
    policy = preparation_policy(metadata['settings'], coverage.coverage_identity)
    filling_policy = FillingPolicy(policy.symmetric_frame_authority, policy.neutral_number_authority,
        'explicit independent full-k slots; 1/Nk per primitive cell; ' + policy.coverage_authority,
        policy.exact_degeneracy_authority, 'Ry', 'electrons per primitive cell')
    filling = ZeroTemperatureFilling(levels, populations, tuple(blocks), metadata['neutral_target'],
        metadata['weighted_number'], tuple(tuple(p) for p in metadata['spin_numbers']),
        tuple(metadata['chemical_potential_interval']), metadata['boundary_log_fugacity'], filling_policy)
    return MaterialPopulationSpecification(reference, charge_bindings,
        tuple(tuple(p) for p in metadata['valence_charges']), (4, 4, 4), sectors,
        tuple((i * 346, (i + 1) * 346) for i in range(128)), filling, policy)
