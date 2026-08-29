"""Exact gauge/material, Newtonian-sign, and shear checks for P248."""

from __future__ import annotations

import sympy as sp

from substrate_framework.optical_gothic import (
    newtonian_optical_gradient_ledger,
    optical_continuity_compatibility,
    optical_shear_decomposition,
)
from substrate_framework.verification import CheckLedger


def run() -> int:
    ledger = CheckLedger("P248-COMPATIBILITY")
    mean = sp.symbols("nbar", positive=True)
    dt, adv, div = sp.symbols("dt adv div", real=True)

    material_sign = optical_continuity_compatibility(
        mean,
        dt,
        adv,
        div,
        flow_orientation=-1,
    )
    ledger.check(
        "material-sign exact residual relation",
        sp.simplify(
            material_sign.harmonic_residual
            - (2 * mean * material_sign.material_residual - mean**2 * div)
        )
        == 0,
    )
    on_material = optical_continuity_compatibility(
        mean,
        -adv - mean * div,
        adv,
        div,
        flow_orientation=-1,
    )
    ledger.check(
        "joint material and harmonic iff incompressible",
        sp.simplify(on_material.harmonic_residual + mean**2 * div) == 0,
    )
    static_inhomogeneous = optical_continuity_compatibility(mean, 0, 0, 0)
    ledger.check(
        "inhomogeneous static counterexample to homogeneity claim",
        static_inhomogeneous.material_residual == 0
        and static_inhomogeneous.harmonic_residual == 0,
    )

    paper_sign = optical_continuity_compatibility(
        mean,
        -adv - mean * div,
        adv,
        div,
        flow_orientation=1,
    )
    ledger.check(
        "printed shift sign gives different compatibility",
        sp.simplify(
            paper_sign.harmonic_residual
            + mean * (4 * adv + 3 * mean * div)
        )
        == 0,
    )

    gradient_squared = sp.symbols("gradient_squared", positive=True)
    energy = newtonian_optical_gradient_ledger(3, 5, gradient_squared)
    ledger.check(
        "printed gradient energy has opposite Newtonian sign",
        sp.simplify(
            energy.paper_energy_density + energy.required_newtonian_energy_density
        )
        == 0,
    )
    ledger.mutation_sensitive(
        "Newtonian sign oracle",
        lambda sign: sp.simplify(
            sign * energy.paper_energy_density
            - energy.required_newtonian_energy_density
        )
        == 0,
        -1,
        [1],
    )

    shear = optical_shear_decomposition([4, 2, 1])
    ledger.check(
        "additive determinant-mean strain is not trace-free",
        shear.additive_trace == 1,
    )
    ledger.check(
        "log shear trace vanishes exactly",
        shear.logarithmic_shear_trace == 0,
    )
    ledger.check(
        "log shear normalized determinant is one",
        shear.normalized_determinant == 1,
    )
    isotropic = optical_shear_decomposition([3, 3, 3])
    ledger.check(
        "log shear exact isotropic reduction",
        isotropic.logarithmic_shear == (0, 0, 0),
    )

    # The claimed on-shell definition cancels the field equation if promoted
    # to an arbitrary-field closure: box_g = -kappa*(T - box_g/kappa).
    box_g, coupling, matter = sp.symbols("box_g coupling matter", nonzero=True)
    substituted_residual = sp.simplify(
        box_g + coupling * (matter - box_g / coupling)
    )
    ledger.check(
        "off-shell self-source substitution is a matter-only tautology",
        substituted_residual == coupling * matter,
    )

    deformation = sp.diag(2, 1, 1)
    jacobian = deformation.det()
    cauchy = sp.eye(3)
    first_piola = jacobian * cauchy * deformation.inv().T
    volume_only = jacobian * cauchy
    ledger.check(
        "volume density alone is not first Piola stress",
        first_piola != volume_only,
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(run())
