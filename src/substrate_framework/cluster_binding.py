"""Matched bulk-interstitial electronic binding of sparse cluster occupancies."""
from __future__ import annotations

from math import isfinite
from typing import Mapping


def occupancy_bindings_eV(cluster_energies: Mapping[int, float],
                          bulk_energy_eV: float, bulk_D_energy_eV: float) -> dict:
    """Return per-added-D block means without inventing uncomputed occupancies.

    Positive binding favors transfer from the matched bulk interstitial to the
    vacancy. These electronic energies are not ionic finite-temperature free energies.
    """
    occupancies = sorted(cluster_energies)
    if (not occupancies or occupancies[0] != 0 or
            any(not isinstance(n, int) or n < 0 for n in occupancies) or
            not all(isfinite(value) for value in
                    [bulk_energy_eV, bulk_D_energy_eV, *cluster_energies.values()])):
        raise ValueError('Require finite matched energies and integer occupancies starting at0')
    insertion = bulk_D_energy_eV-bulk_energy_eV
    blocks = []
    for lower, upper in zip(occupancies[:-1], occupancies[1:]):
        added = upper-lower
        binding = insertion+(cluster_energies[lower]-cluster_energies[upper])/added
        blocks.append(dict(n_from=lower, n_to=upper, added_D=added,
                           binding_eV_per_added_D=binding,
                           interpretation='single increment' if added == 1 else 'block mean'))
    excess = {n: energy-cluster_energies[0]-n*insertion
              for n, energy in cluster_energies.items()}
    return dict(bulk_interstitial_insertion_eV=insertion, blocks=blocks,
                excess_cluster_energies_eV=excess,
                reference='matched_bulk_interstitial_D')
