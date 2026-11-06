"""Expose and verify the Gamma projection/signed-spectrum composition boundary."""

import json
from pathlib import Path
import sys
import traceback

import numpy as np
import sympy as sp

from substrate_framework import material_modes as modes
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0003"))
from oracle_artifacts import sha256, save_receipt


def exact_fixture(ledger):
    spring = sp.Matrix([[7, 2, -1], [2, 5, 3], [-1, 3, 11]])
    force = spring.row_join(-spring).col_join((-spring).row_join(spring))
    inverse_mass_root = sp.diag(*([sp.Rational(1, 2)] * 3 + [sp.Rational(1, 3)] * 3))
    raw = inverse_mass_root * force * inverse_mass_root
    raw[0, 4] += 1 + 2 * sp.I
    hermitian = (raw + raw.H) / 2
    translations = (2 * sp.eye(3)).col_join(3 * sp.eye(3))
    projector = sp.eye(6) - translations * (translations.H * translations).inv() * translations.H
    ledger.check("exact acoustic projector", projector.H == projector and projector**2 == projector)
    ledger.check("exact translation annihilation", projector * translations == sp.zeros(6, 3))
    projected = projector * hermitian * projector
    ledger.check("exact projected Hermiticity", projected == projected.H)
    roots = projected.charpoly().as_poly().all_roots(radicals=False)
    ledger.check("fixture exact null and positive spectrum",
                 sum(root.is_zero is True for root in roots) == 3
                 and sum(root.is_positive is True for root in roots) == 3)
    return projected, roots


def numeric_fixture():
    masses = np.array([4., 9.])
    spring = np.array([[7., 2., -1.], [2., 5., 3.], [-1., 3., 11.]])
    force = np.block([[spring, -spring], [-spring, spring]]).reshape(2, 3, 2, 3)
    scale = modes.RYDBERG_SI / modes.BOHR_RADIUS_SI**2 / modes.AMU_SI
    operator = modes.mass_weighted_operator(force, masses, units="si")
    operator[0, 4] += (1 + 2j) * scale
    return modes.project_hermitian(operator), masses, scale


def composed_spectra(ledger, exact, roots, receipt):
    operator, masses, scale = numeric_fixture()
    original = operator.copy()
    expected = np.sort(np.array([float(root.evalf(40)) for root in roots])) * scale
    tolerance = 256 * len(operator) * np.finfo(float).eps * float(exact.norm().evalf(40)) * scale
    observations = {"engineering_absolute_operator_allowance_si": tolerance,
                    "exact_dimensionless_eigenvalues": [str(root) for root in roots]}
    receipt["observations"] = observations
    for sign in (1, -1):
        projected = modes.project_gamma_translations(sign * operator, masses, q_fractional=[0, 0, 0])
        residuals = modes.operator_residuals(projected, masses)
        observations[str(sign)] = {"raw_projected_residuals_si": residuals}
        spectrum = modes.signed_modes(projected, units="si")
        reference = np.sort(sign * expected)
        ledger.check(f"composed exact-root spectrum sign={sign}",
                     np.max(np.abs(spectrum.eigenvalues - reference)) <= tolerance)
        ledger.check(f"composed translation and eigenpair residual sign={sign}",
                     residuals["translation_right_frobenius"] <= tolerance
                     and spectrum.eigenpair_residual <= tolerance)
        observations[str(sign)]["actual_eigenvalues_si"] = spectrum.eigenvalues.tolist()
        observations[str(sign)]["signed_frequencies_thz"] = spectrum.frequency_thz.tolist()
        if sign < 0:
            ledger.check("three finite instabilities survive projection and spectrum",
                         np.count_nonzero(spectrum.eigenvalues < -tolerance) == 3
                         and np.all(spectrum.frequency_thz[:3] < 0))
    ledger.check("raw Hermitian input unchanged", np.array_equal(operator, original))
    asymmetric = operator.copy()
    asymmetric[0, 4] += 1j * scale
    projected = modes.project_gamma_translations(asymmetric, masses, q_fractional=[0, 0, 0])
    ledger.check("genuine projected asymmetry retained",
                 np.linalg.norm(projected - projected.conj().T) > tolerance)
    try:
        modes.signed_modes(projected, units="si")
    except ValueError:
        rejected = True
    else:
        rejected = False
    ledger.check("genuine asymmetric spectrum remains rejected", rejected)


def main():
    sources = [Path(__file__), HERE / "record.json", Path(modes.__file__),
               HERE.parent / "0003/oracle_artifacts.py"]
    receipt = {"status": "running", "source_sha256": {str(path): sha256(path) for path in sources},
               "versions": {"python": sys.version, "numpy": np.__version__, "sympy": sp.__version__},
               "scope": "Synthetic exact projection/spectrum API regression, not material evidence."}
    ledger = CheckLedger("P254-Gamma-projection-composition")
    try:
        exact, roots = exact_fixture(ledger)
        composed_spectra(ledger, exact, roots, receipt)
        ledger.finish()
        receipt["status"] = "complete"
    except Exception:
        receipt["status"] = "failed"
        receipt["raw_exception"] = traceback.format_exc()
        print("P254_GAMMA_PROJECTION_COMPOSITION_FAILED", flush=True)
        raise
    finally:
        receipt["passed_checks"] = ledger.passed
        save_receipt(HERE / "receipt.json", receipt)
    print("P254_GAMMA_PROJECTION_COMPOSITION_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
