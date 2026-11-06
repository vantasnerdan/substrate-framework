#!/usr/bin/env python3
"""Run independent finite mathematical energy oracles; never actual Ni propagation.

Usage: python run_mixed_excitation_energy_oracle.py --output /path/report.json
The report contains actual source bytes/SHA/raw text, exact differentiation,
dense finite-periodic energy/refinement, normalization and targeted mutations.
Only the caller executes this file; importing it runs no checks.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as sy

from joint_mean_energy_components import components as authoritative_components
from mixed_excitation_energy import mixed_components
from mixed_excitation_energy_oracle_support import (
    FIELDS, KEYS, exact_derivatives, fixture, local_impulse_normalization,
    periodic_model, periodic_scalar, physical_laplacian,
)


def bindings():
    root = Path(__file__).resolve().parent
    names = ("mixed_excitation_energy.py", "run_mixed_excitation_energy_oracle.py",
             "mixed_excitation_energy_oracle_support.py", "joint_mean_energy_components.py",
             "moving-coherent-energy-construction-registration-001.json",
             "propagated-excitation-energy-source-grounding-001.json")
    result = []
    for name in names:
        path = root/name
        raw = path.read_bytes()
        result.append(dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(),
                           raw_source=raw.decode("utf-8")))
    return result


def number(value):
    z = complex(value)
    return dict(real=z.real, imag=z.imag)


def maximum_defect(actual, expected):
    return max(abs(complex(actual[key]-expected[key])) for key in expected)


def exact_check():
    args = fixture(exact=True)
    expected, source_defects = exact_derivatives(args, authoritative_components)
    actual = mixed_components(**args)
    defects = {key: sy.simplify(actual[key]-value) for key, value in expected.items()}
    return dict(pass_check=all(value == 0 for value in defects.values())
                           and all(value == 0 for value in source_defects.values()),
                exact_derivative_defects={key: str(value) for key, value in defects.items()},
                literal_vs_authoritative_defects={key: str(value) for key, value in source_defects.items()},
                components={key: number(sy.N(value)) for key, value in actual.items()},
                complex_defect_preserved=any(sy.im(value) != 0 for value in actual.values()))


def dense_check(gamma, *, matched=False):
    args = fixture(gamma=gamma)
    phase, cp, cm, dense = periodic_model(args, mode=0 if gamma else 1)
    if matched:
        args["coefficients"].delta_curvature = np.zeros_like(args["coefficients"].delta_curvature)
    args.update(curvature_plus=cp, curvature_minus=cm)
    expected = mixed_components(**args)
    energy = lambda lam, eta: periodic_scalar(args, phase, dense, lam, eta)
    refinements = []
    numerical = []
    for h in (0.04, 0.02, 0.01, 0.005):
        sample = physical_laplacian(energy, h)
        numerical.append(sample)
        refinements.append(dict(h=h, max_component_error=maximum_defect(sample, expected),
                                component_errors={key: abs(complex(sample[key]-expected[key]))
                                                  for key in expected}))
    extrapolated = {key: (4*numerical[-1][key]-numerical[-2][key])/3 for key in expected}
    scale = max(1, max(abs(complex(value)) for value in expected.values()))
    # Independently pair full dense Fourier vectors, with no guessed N factor.
    ss = args["stages"]
    mean = np.tile(ss["mean"].u, len(phase))
    mixed = np.tile(ss["mixed"].u, len(phase))
    plus = np.concatenate([v*ss["plus"].u for v in phase])
    minus = np.concatenate([v.conjugate()*ss["minus"].u for v in phase])
    dense_mixed = (mean@dense@mixed+mixed@dense@mean+plus@dense@minus+minus@dense@plus)/2
    pair_defect = abs(complex(dense_mixed-expected["ionic_C"]))
    physical = energy(0.02+0.01j, 0.02-0.01j)
    physical_imaginary = max(abs(complex(value).imag) for value in physical.values())
    final_error = maximum_defect(extrapolated, expected)
    improving = refinements[-1]["max_component_error"] < refinements[0]["max_component_error"]/40
    return dict(pass_check=final_error < 2e-7*scale and pair_defect < 1e-10*scale
                           and physical_imaginary < 1e-10*scale and improving,
                gamma=gamma, matched_delta_C_zero=matched, finite_model_cells=len(phase), refinements=refinements,
                richardson_max_component_error=final_error, tolerance=2e-7*scale,
                dense_fourier_curvature_pair_defect=pair_defect,
                physical_eta_equals_conjugate_lambda_imaginary_defect=physical_imaginary,
                independent_curvature_difference=float(np.max(np.abs(cp-cm))),
                nonstationary_mean_gradient_norm=float(np.linalg.norm(ss["mean"].gradient)),
                nonzero_mixed_density_norm=float(np.linalg.norm(ss["mixed"].density)),
                nonzero_mixed_one_body=number(args["one_body"]["mixed"]),
                components={key: number(value) for key, value in expected.items()}), args, expected

def mutation_checks(args, expected):
    """Each mutation is evaluated against the independently validated target."""
    mutations = {}
    def record(name, wrong):
        gap = maximum_defect(wrong, expected)
        mutations[name] = dict(detected=gap > 1e-7, max_component_gap=gap,
                               changed_components=[key for key in expected
                                                   if abs(complex(wrong[key]-expected[key])) > 1e-7])
    omitted = copy.deepcopy(args)
    for name in FIELDS:
        setattr(omitted["stages"]["mixed"], name,
                np.zeros_like(getattr(omitted["stages"]["mixed"], name)))
    omitted["one_body"]["mixed"] = 0
    record("omit_all_second_mean_data", mixed_components(**omitted))
    record("omit_second_electronic_mean", {**expected, "bare_electronic": 0})
    record("half_mixed_derivative", {key: value/2 for key, value in expected.items()})
    record("multiply_global_maxwell_by_N", {**expected, "maxwell": expected["maxwell"]*args["reference"].ncells})
    wrong_q = {**args, "curvature_plus": args["curvature_minus"],
               "curvature_minus": args["curvature_plus"]}
    record("swap_q_minusq_orientation", mixed_components(**wrong_q))
    r, ss = args["reference"], args["stages"]
    # Removing only canonical p--Za cross energy, keeping p² and (Za)².
    z = np.repeat(args["charges"], 3)
    cross = 0j
    for left, right in (("mean", "mixed"), ("mixed", "mean"),
                        ("plus", "minus"), ("minus", "plus")):
        cross += np.sum(ss[left].p*z*np.tile(ss[right].shift, len(args["charges"]))/r.masses)
    record("omit_minimal_kinetic_p_Za_cross", {**expected,
           "ionic_kinetic_minimal": expected["ionic_kinetic_minimal"]+r.ncells*cross})
    record("double_matched_physical_C", {**expected, "ionic_C": 2*expected["ionic_C"]})
    return mutations


def input_checks():
    cases = {}
    def rejected(name, args):
        try:
            mixed_components(**args)
        except (ValueError, TypeError) as error:
            cases[name] = dict(rejected=True, error=str(error))
        else:
            cases[name] = dict(rejected=False)
    args = fixture()
    del args["stages"]["mixed"]
    rejected("missing_mixed_stage", args)
    args = fixture()
    args["one_body"]["mixed"] = None
    rejected("missing_mixed_electronic_moment", args)
    args = fixture()
    args["stages"]["mixed"].density = None
    rejected("missing_mixed_density", args)
    args = fixture(gamma=False)
    args["stages"]["plus"].shift[0] = 1
    rejected("nonGamma_first_Maxwell_field", args)
    return dict(pass_check=all(case["rejected"] for case in cases.values()), cases=cases)


def run():
    exact = exact_check()
    gamma, ga, ge = dense_check(True)
    non_gamma, qa, qe = dense_check(False)
    matched, ma, me = dense_check(False, matched=True)
    mutations = mutation_checks(ma, me)
    # A Gamma fixture separately exercises BOTH first uniform-field directions.
    mutations["Gamma_omit_minimal_cross"] = mutation_checks(ga, ge)["omit_minimal_kinetic_p_Za_cross"]
    normalization = local_impulse_normalization(mixed_components)
    normalization["pass_check"] = all(value < 1e-12 for key, value in normalization.items()
                                            if key.endswith("defect"))
    normalization["mutations_detected"] = (normalization["wrong_collective_units_dose_gap"] > 1e-4
                                            and normalization["missing_energy_half_gap"] > 1e-4)
    inputs = input_checks()
    passed = (exact["pass_check"] and exact["complex_defect_preserved"]
              and gamma["pass_check"] and non_gamma["pass_check"] and matched["pass_check"]
              and me["delta_C"] == 0 and inputs["pass_check"]
              and normalization["pass_check"] and normalization["mutations_detected"]
              and all(case["detected"] for case in mutations.values()))
    return dict(pass_check=passed, exact_sympy=exact, dense_periodic_gamma=gamma,
                dense_periodic_nonGamma=non_gamma, dense_periodic_matched=matched,
                local_impulse_normalization=normalization,
                mutations=mutations, required_inputs=inputs)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    options = parser.parse_args()
    report = dict(schema="P254-mixed-excitation-energy-oracle-v1", source_bindings=bindings(),
                  scope="Finite mathematical fixtures only; no Ni propagation, finite-dose licence, heat or ranking.",
                  derivative="d2/(dlambda deta), no half; bilinear independent q/-q; physical eta=lambda*",
                  units="Physical per-cell matter amplitudes; GLOBAL Maxwell coordinate/momentum and mass")
    try:
        report.update(run())
    except Exception as error:
        report.update(pass_check=False, execution_error=dict(type=type(error).__name__, message=str(error)))
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(output=str(options.output.resolve()), pass_check=report["pass_check"])))
    return 0 if report["pass_check"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
