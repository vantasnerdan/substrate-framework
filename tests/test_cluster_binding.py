import pytest

from substrate_framework.cluster_binding import occupancy_bindings_eV


def test_sparse_binding_uses_added_particle_count_and_conserves_species_gauge():
    energies = {0: -90.0, 1: -94.2, 2: -98.6, 4: -107.4, 6: -116.0}
    result = occupancy_bindings_eV(energies, -100.0, -104.0)
    assert [r['binding_eV_per_added_D'] for r in result['blocks']] == pytest.approx([.2, .4, .4, .3])
    # Arbitrary species zero shifts must cancel in a balanced trapping reaction.
    metal_zero, D_zero = 7.3, -2.1
    shifted = {n: energy+7*metal_zero+n*D_zero for n, energy in energies.items()}
    alternate = occupancy_bindings_eV(shifted, -100+8*metal_zero, -104+8*metal_zero+D_zero)
    assert [r['binding_eV_per_added_D'] for r in alternate['blocks']] == pytest.approx(
        [r['binding_eV_per_added_D'] for r in result['blocks']], abs=1e-12)
    assert alternate['excess_cluster_energies_eV'] == pytest.approx(result['excess_cluster_energies_eV'])
