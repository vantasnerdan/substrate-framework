"""Conditional fixed-mean-number FD limit of explicitly labelled mode energies.

Exact supplied float equality, not a near-degeneracy threshold, defines blocks.
The single common multiplier leaves longitudinal magnetization free. Measures
count independent modes per declared number unit; they are not probabilities
or fluctuation coefficients. No temperature, canonical-N or equilibrium claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import groupby
from math import fsum, isfinite, log


@dataclass(frozen=True)
class ModeLevel:
    identity: str
    spin: int
    energy: float
    measure: float


@dataclass(frozen=True)
class FillingPolicy:
    """Source labels for the chosen operator, mean-number and mode measure."""
    operator_authority: str
    number_authority: str
    measure_authority: str
    degeneracy_authority: str
    energy_units: str
    number_units: str

    def __post_init__(self):
        if any(not isinstance(value, str) or not value.strip()
               for value in vars(self).values()):
            raise ValueError('all filling policy source labels are required')


@dataclass(frozen=True)
class DegenerateBlock:
    energy: float
    indices: tuple[int, ...]
    measure: float
    population: float


@dataclass(frozen=True)
class ZeroTemperatureFilling:
    levels: tuple[ModeLevel, ...]
    populations: tuple[float, ...]
    blocks: tuple[DegenerateBlock, ...]
    target_mean_number: float
    weighted_number: float
    spin_numbers: tuple[tuple[int, float], ...]
    chemical_potential_interval: tuple[float | None, float | None]
    boundary_log_fugacity: float | None
    policy: FillingPolicy


def fill_zero_temperature(levels, target_mean_number, *, policy: FillingPolicy):
    """Return scalar block populations in input order, with no spin factor.

    At a fractional block e*, mu(theta)=e*+theta*log(p/(1-p))+o(theta).
    For a gapped integer fill the returned interval bounds possible limiting
    multipliers, not a uniquely selected mu. None denotes an infinite bound.
    At an empty/full retained space the multiplier can lie below/above its
    spectrum. The initially empty external complement is a separate authority.
    This function only defines the mean-number limit, not a finite-T solver.
    """
    if not isinstance(policy, FillingPolicy):
        raise TypeError('explicit FillingPolicy required')
    levels = tuple(levels)
    if not levels or len({level.identity for level in levels}) != len(levels):
        raise ValueError('nonempty distinct mode identities required')
    for level in levels:
        if (not isinstance(level.identity, str) or not level.identity.strip()
                or type(level.spin) is not int or level.spin not in (1, 2)
                or not isfinite(level.energy) or not isfinite(level.measure)
                or level.measure <= 0):
            raise ValueError('finite energy, positive measure and explicit spin 1/2 required')
    target = float(target_mean_number)
    capacity = fsum(level.measure for level in levels)
    if not isfinite(target) or not 0 <= target <= capacity:
        raise ValueError('mean-number target outside retained-mode capacity')
    populations = [0.0] * len(levels)
    blocks, filled_measures = [], []
    interval, log_fugacity = (None, None), None
    ordered = sorted(range(len(levels)), key=lambda i: levels[i].energy)
    for energy, members in groupby(ordered, key=lambda i: levels[i].energy):
        indices = tuple(members)
        measure = fsum(levels[i].measure for i in indices)
        below = fsum(filled_measures)
        through = fsum((*filled_measures, measure))
        if target >= through:
            population = 1.0
        elif target <= below:
            population = 0.0
        else:
            population = (target - below) / measure
        blocks.append(DegenerateBlock(energy, indices, measure, population))
        for i in indices:
            populations[i] = population
        filled_measures.append(measure)
    fractional = [block for block in blocks if 0 < block.population < 1]
    if fractional:
        block = fractional[0]
        interval = (block.energy, block.energy)
        log_fugacity = log(block.population) - log(1 - block.population)
    else:
        occupied = [block.energy for block in blocks if block.population == 1]
        empty = [block.energy for block in blocks if block.population == 0]
        interval = (occupied[-1] if occupied else None, empty[0] if empty else None)
    spin_numbers = tuple((spin, fsum(level.measure * populations[i]
                                  for i, level in enumerate(levels) if level.spin == spin))
                         for spin in sorted({level.spin for level in levels}))
    return ZeroTemperatureFilling(
        levels, tuple(populations), tuple(blocks), target,
        fsum(level.measure * populations[i] for i, level in enumerate(levels)),
        spin_numbers, interval, log_fugacity, policy)
