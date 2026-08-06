from __future__ import annotations

import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.local_horizon_gravity import (
    LocalHorizonEinsteinLedger,
    local_horizon_einstein_ledger,
)


def test_public_package_exports_local_horizon_gravity_api() -> None:
    assert framework.LocalHorizonEinsteinLedger is LocalHorizonEinsteinLedger
    assert framework.local_horizon_einstein_ledger is local_horizon_einstein_ledger


def test_clausius_and_raychaudhuri_fix_the_full_einstein_coefficient() -> None:
    entropy_density, action, speed, cosmological = sp.symbols(
        "eta_A J c Lambda", positive=True
    )
    ledger = local_horizon_einstein_ledger(
        entropy_density,
        action,
        speed,
        cosmological_constant=cosmological,
    )
    assert ledger.unruh_energy_per_inverse_length == action * speed / (2 * sp.pi)
    assert ledger.clausius_stress_coefficient == 2 * sp.pi / (action * speed)
    assert ledger.einstein_stress_coupling == 2 * sp.pi / (
        action * speed * entropy_density
    )
    assert ledger.newton_constant == speed**3 / (4 * action * entropy_density)
    assert ledger.einstein_stress_coupling == (
        8 * sp.pi * ledger.newton_constant / speed**4
    )
    assert ledger.cosmological_constant == cosmological
    assert ledger.radiative_tensor_polarizations == 2


def test_action_and_area_entropy_are_load_bearing_not_independent_aliases() -> None:
    baseline = local_horizon_einstein_ledger(3, 5, 7)
    doubled_action = local_horizon_einstein_ledger(3, 10, 7)
    doubled_entropy = local_horizon_einstein_ledger(6, 5, 7)
    assert doubled_action.newton_constant == baseline.newton_constant / 2
    assert doubled_entropy.newton_constant == baseline.newton_constant / 2


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ((1.0, 1, 1), "entropy_area_density"),
        ((1, 0, 1), "action_scale"),
        ((1, 1, -1), "signal_speed"),
    ],
)
def test_local_horizon_inputs_must_be_exact_and_positive(arguments, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        local_horizon_einstein_ledger(*arguments)
