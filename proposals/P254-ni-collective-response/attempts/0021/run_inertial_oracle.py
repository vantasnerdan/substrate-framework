#!/usr/bin/env python3
"""Bounded finite coupled-action oracle; no Native/PW execution or material claim.

Run with --output NEW.json. Output is opened exclusively before numerical work.
All failures are retained; a failed section does not suppress other sections.
"""
from __future__ import annotations

import os
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                  "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_variable] = "1"

import argparse
import hashlib
import json
from pathlib import Path
import signal
import sys
import time
import traceback

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, expm, solve
import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from substrate_framework.inertial_fermion_dynamics import (
    FermionChannel, InertialFermionDynamics, ionic_masses_ry,
)
from substrate_framework.screened_force_response import screened_response

ALLOWANCE = 1e-10
PREREGISTRATION = "2cb8638a"
HERE = Path(__file__).resolve().parent


class Checks:
    def __init__(self):
        self.records = []

    def close(self, name, actual, expected):
        actual, expected = np.asarray(actual), np.asarray(expected)
        error = float(np.linalg.norm(actual - expected) / max(1.0, np.linalg.norm(expected)))
        self.records.append({"name": name, "scaled_error": error if np.isfinite(error) else None,
                             "allowance": ALLOWANCE, "passed": bool(np.isfinite(error) and error <= ALLOWANCE)})

    def truth(self, name, value, **details):
        self.records.append({"name": name, "passed": bool(value), **details})

    def exposed(self, name, wrong, right):
        separation = float(np.linalg.norm(np.asarray(wrong)-np.asarray(right))
                           / max(1.0, np.linalg.norm(right)))
        self.truth(name, np.isfinite(separation) and separation > 100*ALLOWANCE,
                   scaled_separation=separation if np.isfinite(separation) else None,
                   required_separation=100*ALLOWANCE)

    def section(self, name, operation):
        try:
            operation()
        except TimeoutError:
            raise
        except Exception as exc:
            self.truth(name + "/exception", False, exception=repr(exc), traceback=traceback.format_exc())


def binding(path):
    path = Path(path)
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def source_bindings(checks):
    inventory = json.loads((HERE / "source-inventory.json").read_text())
    names = {"variational_fermion_static.py", "screened_force_response.py",
             "finite_fermion_response.py", "constants.f90", "dynamics_module.f90",
             "dyndia.f90", "constructive-dynamics-source.original.txt",
             "constructive-dynamics-source-correction.json"}
    rows = []
    for registered in inventory["files"]:
        if Path(registered["path"]).name in names:
            current = binding(registered["path"])
            checks.truth("source/" + Path(registered["path"]).name,
                         current["sha256"] == registered["sha256"]
                         and current["bytes"] == int(registered["bytes"]))
            rows.append(current)
    rows.extend(binding(path) for path in (HERE / "registration.json", HERE / "source-inventory.json",
        Path(__file__), ROOT / "src/substrate_framework/inertial_fermion_dynamics.py"))
    constants = Path(next(row["path"] for row in rows if Path(row["path"]).name == "constants.f90")).read_text()
    checks.truth("source/QE_Ry_mass_is_half_Hartree_mass", "AMU_RY           = AMU_AU / 2.0_DP" in constants)
    registration = json.loads((HERE / "registration.json").read_text())
    checks.truth("source/preregistered_inventory", binding(HERE / "source-inventory.json")["sha256"]
                 == registration["source_inventory"]["sha256"])
    return rows


def symbolic_checks(checks):
    mass, frozen, gap, population, coupling, kernel, z = sp.symbols("M C d f g K z", nonzero=True, real=True)
    generator = sp.Matrix([[0, 1/mass, 0, 0], [-frozen, 0, -2*coupling, 0],
        [0, 0, 0, -gap], [population*coupling, 0, gap+2*population*kernel, 0]])
    metric = sp.Matrix([[frozen, 0, 2*coupling, 0], [0, 1/mass, 0, 0],
        [2*coupling, 0, 2*gap/population+4*kernel, 0], [0, 0, 0, 2*gap/population]])
    checks.truth("symbolic/same_action_energy", sp.simplify(generator.T*metric+metric*generator) == sp.zeros(4))
    x_over_u = -gap*population*coupling/(gap*(gap+2*population*kernel)-z*z)
    electronic = -sp.I*z*sp.eye(2)-generator[2:, 2:]
    solved = electronic.inv()*generator[2:, :1]
    checks.truth("symbolic/electronic_retarded_sign", sp.simplify(solved[0]-x_over_u) == 0)
    ionic_inverse = frozen + 2*coupling*x_over_u - z*z*mass
    equations = (-sp.I*z*sp.eye(4)-generator)*sp.Matrix([1, -sp.I*z*mass, solved[0], solved[1]])
    checks.truth("symbolic/inertial_frequency_inverse",
                 all(sp.simplify(value) == 0 for value in (equations-sp.Matrix([0, ionic_inverse, 0, 0]))))
    wrong = generator.copy()
    wrong[1, 2] *= -1
    checks.truth("symbolic/wrong_force_sign_exposed", sp.simplify(wrong.T*metric+metric*wrong) != sp.zeros(4))
    amu, electron = sp.symbols("amu electron", positive=True)
    checks.truth("symbolic/Ry_mass_twofactor", sp.simplify((amu/(2*electron))/(amu/electron)) == sp.Rational(1, 2))


def fixture():
    rng = np.random.default_rng(2250021)
    channels, raw = [], []
    for spin, weight, energies, occupations in (
        (1, .7, [-.8, .2, 1.3], [.82, .41, .08]),
        (-1, 1.3, [-.6, .4, 1.7], [.73, .28, .04]),
    ):
        frame, _ = np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))
        def operators(count, scale):
            a = rng.normal(size=(count, 3, 3))+1j*rng.normal(size=(count, 3, 3))
            return scale*(a+a.conj().transpose(0, 2, 1))/2
        g = operators(2, .23)
        local = operators(1, .35)[0]
        densities = np.zeros((2, 3, 3), dtype=complex)
        densities[0 if spin == 1 else 1] = local
        h = (frame*np.asarray(energies)) @ frame.conj().T
        rho = (frame*np.asarray(occupations)) @ frame.conj().T
        channels.append(FermionChannel(h, occupations, frame, g, densities, weight, spin))
        raw.append(dict(h=h, rho=rho, g=g, d=densities, weight=weight, spin=spin,
                        frame=frame, e=np.asarray(energies), f=np.asarray(occupations)))
    # Actual observed species masses/AMU_RY, but two scalar test coordinates.
    # This finite electronic fixture is not the actual material electronic model.
    mass = ionic_masses_ry([58.6934, 1.008], [1, 2], 911.4442431086565, cartesian=1)
    kernel = np.array([[.31, -.09], [-.09, .27]])
    frozen = np.array([[1.8, -.24], [-.24, 1.1]])  # independently supplied, not screened fit
    return InertialFermionDynamics(channels, kernel, frozen, mass), raw


def original_rhs(raw, kernel, frozen, masses, state, force=None):
    """Independent original-basis full EOM, no implementation observables/pack."""
    u, p = state[:2], state[2:4]
    rhos = [state[4+9*a:13+9*a].reshape(3, 3) for a in range(len(raw))]
    density = np.array([sum(c["weight"]*np.trace(c["d"][s] @ rho)
                           for c, rho in zip(raw, rhos)) for s in range(2)])
    gradient = np.array([sum(c["weight"]*np.trace(c["g"][a] @ rho)
                            for c, rho in zip(raw, rhos)) for a in range(2)])
    derivative = []
    for c, rho in zip(raw, rhos):
        v = sum(u[a]*c["g"][a] for a in range(2)) + sum((kernel @ density)[s]*c["d"][s] for s in range(2))
        derivative.append((-1j*(c["h"] @ rho-rho @ c["h"]+v @ c["rho"]-c["rho"] @ v)).ravel())
    return np.concatenate([p/masses, -frozen @ u-gradient+(0 if force is None else force), *derivative])


def original_generator(raw, kernel, frozen, masses):
    """Independent Kronecker/outer-product assembly; row-major vectorization."""
    sizes = [len(c["e"])**2 for c in raw]
    electronic = block_diag(*[-1j*(np.kron(c["h"], np.eye(3))-np.kron(np.eye(3), c["h"].T)) for c in raw])
    density_map = np.concatenate([c["weight"]*c["d"].transpose(0, 2, 1).reshape(2, 9) for c in raw], axis=1)
    gradient_map = np.concatenate([c["weight"]*c["g"].transpose(0, 2, 1).reshape(2, 9) for c in raw], axis=1)
    drive_g = np.concatenate([np.column_stack([(-1j*(op @ c["rho"]-c["rho"] @ op)).ravel() for op in c["g"]]) for c in raw])
    drive_d = np.concatenate([np.column_stack([(-1j*(op @ c["rho"]-c["rho"] @ op)).ravel() for op in c["d"]]) for c in raw])
    electronic += drive_d @ kernel @ density_map
    total = sum(sizes)
    return np.block([[np.zeros((2, 2)), np.diag(1/masses), np.zeros((2, total))],
                     [-frozen, np.zeros((2, 2)), -gradient_map],
                     [drive_g, np.zeros((total, 2)), electronic]])


def independent_energy(raw, kernel, frozen, masses, state):
    u, p = state[:2].real, state[2:4].real
    density, gradient, bare = np.zeros(2, complex), np.zeros(2, complex), 0.0
    for a, c in enumerate(raw):
        rho = state[4+9*a:13+9*a].reshape(3, 3)
        r = c["frame"].conj().T @ rho @ c["frame"]
        for i in range(3):
            for j in range(i+1, 3):
                bare += c["weight"]*(c["e"][j]-c["e"][i])/(c["f"][i]-c["f"][j])*abs(r[i, j])**2
        density += c["weight"]*np.array([np.trace(d @ rho) for d in c["d"]])
        gradient += c["weight"]*np.array([np.trace(g @ rho) for g in c["g"]])
    return dict(ionic_kinetic=float(.5*np.sum(p*p/masses)), ionic_frozen=float(.5*u @ frozen @ u),
                electronic_orbit=float(bare), electron_ion=float((u @ gradient).real),
                density_interaction=float((.5*density @ kernel @ density).real))


def susceptibility(raw, z):
    """Settled +i susceptibility to -B*u; physical +g curvature is minus screened chi."""
    result = np.zeros((4, 4), complex)
    for c in raw:
        operators = np.array([c["frame"].conj().T @ op @ c["frame"] for op in np.concatenate([c["g"], c["d"]])])
        for i in range(3):
            for j in range(3):
                if i != j:
                    result += c["weight"]*(c["f"][i]-c["f"][j])/(c["e"][j]-c["e"][i]-z)*np.outer(operators[:, i, j], operators[:, j, i])
    return result


def numerical_checks(checks, details):
    model, raw = fixture()
    kernel, frozen, masses = model.interaction, model.frozen_curvature, model.masses
    independent = original_generator(raw, kernel, frozen, masses)
    generator = model.constant_generator()
    checks.close("EOM/independent_original_basis_blocks", generator, independent)
    basis = model.tangent_basis()
    rng = np.random.default_rng(22521)
    coordinates = rng.normal(size=basis.shape[1])*.07
    coordinates[2:4] *= np.sqrt(masses)
    initial = basis @ coordinates
    checks.close("EOM/direct_full_rhs", model.derivative(initial), original_rhs(raw, kernel, frozen, masses, initial))
    checks.close("mass/unequal_source_values", masses, np.array([58.6934, 1.008])*911.4442431086565)
    checks.truth("mass/unequal", masses[0] != masses[1])
    checks.close("mass/one_based_species_and_cartesian", ionic_masses_ry([58.6934, 1.008], [2, 1], 911.4442431086565),
                 np.repeat(np.array([1.008, 58.6934])*911.4442431086565, 3))
    _, gradient = model.observables(initial)
    checks.close("force/physical_negative_gradient", model.derivative(initial)[2:4], -frozen @ initial[:2]-gradient)
    checks.exposed("mutation/force_sign", -frozen @ initial[:2]+gradient, model.derivative(initial)[2:4])
    for name, wrong_mass in (("extra_factor_two", 2*masses), ("uniform_mass", np.repeat(masses[0], 2))):
        wrong = original_generator(raw, kernel, frozen, wrong_mass)
        checks.exposed("mutation/"+name, wrong @ initial, generator @ initial)
    wrong_sign = generator.copy()
    wrong_sign[2:4, 4:] *= -1
    no_coupling = generator.copy()
    no_coupling[4:, :2] = 0
    no_cross = original_generator(raw, np.diag(np.diag(kernel)), frozen, masses)
    checks.exposed("mutation/dropped_electronic_ionic_drive", no_coupling @ initial, generator @ initial)
    checks.exposed("mutation/dropped_cross_spin_kernel", no_cross @ initial, generator @ initial)
    unit_weights = [dict(c, weight=1.0) for c in raw]
    checks.exposed("mutation/implicit_spin_weights", original_generator(unit_weights, kernel, frozen, masses) @ initial, generator @ initial)

    times = np.linspace(0, 1.7, 6)
    trajectory = solve_ivp(lambda t, y: original_rhs(raw, kernel, frozen, masses, y),
                           (0, times[-1]), initial, t_eval=times, method="DOP853", rtol=2e-13, atol=2e-14)
    checks.truth("time/independent_scipy_solver_success", trajectory.success, message=trajectory.message)
    checks.truth("time/all_samples", trajectory.y.shape[1] == len(times))
    energy_initial = model.energy_terms(initial)
    for index, tau in enumerate(times[:trajectory.y.shape[1]]):
        state = model.propagate(initial, tau)
        checks.close(f"time/direct_EOM/{index}", state, trajectory.y[:, index])
        checks.close(f"time/independent_exponential/{index}", state, expm(tau*independent) @ initial)
        checks.close(f"time/populations/{index}", model.conserved_populations(state)[0], model.conserved_populations(initial)[0])
        checks.close(f"time/particle_spin/{index}", model.conserved_populations(state)[1:], model.conserved_populations(initial)[1:])
        checks.close(f"energy/conserved/{index}", model.energy_terms(state)["total"], energy_initial["total"])
        terms = independent_energy(raw, kernel, frozen, masses, state)
        for name, expected in terms.items():
            checks.close(f"energy/source_term/{index}/{name}", model.energy_terms(state)[name], expected)
    # Nonzero populations are allowed by EOM, but not silently added to orbit energy.
    population_state = initial.copy()
    population_state[4:13] += .04*np.eye(3).ravel()
    checks.close("population/nonzero_trace_stationary", model.conserved_populations(model.derivative(population_state))[0], [0, 0])
    checks.close("population/weighted_number_spin_definition", model.conserved_populations(population_state)[1:], [.7*.12, .7*.12])

    # Quadratic metric is independently assembled from source-term evaluations;
    # no trajectory or work-integral fitting enters its definition.
    def energy(x):
        return sum(independent_energy(raw, kernel, frozen, masses, basis @ x).values())
    dimension = basis.shape[1]
    identity = np.eye(dimension)
    metric = np.empty((dimension, dimension))
    for i in range(dimension):
        metric[i, i] = 2*energy(identity[i])
        for j in range(i):
            metric[i, j] = metric[j, i] = energy(identity[i]+identity[j])-energy(identity[i])-energy(identity[j])
    tangent_generator = basis.conj().T @ generator @ basis
    checks.close("energy/tangent_closed", generator @ basis, basis @ tangent_generator)
    checks.close("energy/real_tangent_generator", tangent_generator.imag, np.zeros_like(tangent_generator.real))
    checks.close("energy/exact_metric_generator_identity", tangent_generator.T @ metric+metric @ tangent_generator, np.zeros_like(metric))
    for name, wrong in (("wrong_force_sign", wrong_sign), ("dropped_drive", no_coupling)):
        projected = basis.conj().T @ wrong @ basis
        checks.exposed("energy/mutation/"+name, projected.T @ metric+metric @ projected, np.zeros_like(metric))
    velocity = initial[2:4].real/masses
    force = np.array([.4, -.3])
    delta = model.derivative(initial, force)-model.derivative(initial)
    checks.close("energy/external_power", coordinates @ metric @ (basis.conj().T @ delta), velocity @ force)

    for index, z in enumerate((.37+.19j, 1.21+.27j)):
        cel = -screened_response(susceptibility(raw, z), 2, kernel)
        expected_inverse = frozen+cel-z*z*np.diag(masses)
        checks.close(f"frequency/susceptibility_curvature/{index}", model.electronic_curvature(z), cel)
        checks.close(f"frequency/inertia_inverse/{index}", model.frequency_inverse(z), expected_inverse)
        drive = np.array([.2+.1j, -.3+.07j])
        state = model.driven_state(z, drive)
        reference_source = np.r_[np.zeros(2), drive, np.zeros(model.state_size-4)]
        reference = solve(-1j*z*np.eye(model.state_size)-independent, reference_source)
        checks.close(f"frequency/full_original_basis/{index}", state, reference)
        checks.close(f"frequency/composed_displacement/{index}", state[:2], solve(expected_inverse, drive))
        checks.close(f"frequency/physical_EOM_residual/{index}", -1j*z*state, original_rhs(raw, kernel, frozen, masses, state, drive))
        checks.exposed(f"frequency/wrong_inertia_sign/{index}", frozen+cel+z*z*np.diag(masses), expected_inverse)
        checks.exposed(f"frequency/wrong_curvature_sign/{index}", frozen-cel-z*z*np.diag(masses), expected_inverse)

    # An explicitly prepared classical real Gaussian, independent of the mean.
    loading = rng.normal(size=(dimension, dimension))*.025
    loading[2:4] *= np.sqrt(masses)[:, None]
    coordinate_covariance = loading @ loading.T
    covariance = model.prepared_covariance(coordinate_covariance)
    checks.close("covariance/declared_preparation", covariance, (basis @ loading) @ (basis @ loading).conj().T)
    tau = .83
    propagated = model.propagate_covariance(covariance, tau)
    independent_loading = expm(tau*independent) @ basis @ loading
    checks.close("covariance/independent_loading_propagation", propagated, independent_loading @ independent_loading.conj().T)
    checks.close("covariance/coordinate_recovery", basis.conj().T @ covariance @ basis, coordinate_covariance)
    checks.truth("covariance/positive_semidefinite", np.linalg.eigvalsh(propagated).min() > -ALLOWANCE)
    final_coordinates = basis.conj().T @ propagated @ basis
    checks.close("covariance/ensemble_energy_conservation", np.trace(metric @ final_coordinates), np.trace(metric @ coordinate_covariance))
    # Deterministic symmetric sigma ensemble exactly realizes this covariance;
    # no Monte Carlo tolerance is substituted for the preregistered allowance.
    sigma = np.concatenate([np.sqrt(dimension)*loading, -np.sqrt(dimension)*loading], axis=1)
    sigma_states = expm(tau*independent) @ basis @ sigma
    checks.close("covariance/zero_mean_preparation", np.mean(sigma_states, axis=1), np.zeros(model.state_size))
    checks.close("covariance/explicit_ensemble", sigma_states @ sigma_states.conj().T/(2*dimension), propagated)
    checks.exposed("covariance/missing_adjoint_mutation", independent_loading @ independent_loading.T, propagated)
    details.update({"state_size": model.state_size, "real_tangent_dimension": dimension,
        "masses_ry": masses.tolist(), "time_samples": times.tolist(), "initial_energy_terms": energy_initial,
        "preparation": "real Gaussian tangent coordinates; seed 22521; covariance L L^T, displacement/electronic scale .025 and momentum rows scaled sqrt(M); zero ensemble mean independent of tested nonzero first moment",
        "covariance_is_coupled_equilibrium": False,
        "energy_metric_min_eigenvalue": float(np.linalg.eigvalsh(metric).min()),
        "ivp_rhs_evaluations": trajectory.nfev})


def null_sector_checks(checks):
    # Both-zero and equal nonzero populations: full EOM retained, no dividing by zero.
    g = np.array([[[0, 1, .2], [1, 0, .3], [.2, .3, 0]]], complex)
    for label, energies, occupations in (
        ("zero_population", [0, 1, 2], [1, 0, 0]),
        ("equal_population", [0, 1, 2], [.7, .3, .3]),
        ("energy_degeneracy", [0, 1, 1], [.7, .3, .3]),
        ("unequal_population_degenerate_energy", [0, 1, 1], [.7, .4, .2]),
    ):
        channel = FermionChannel(np.diag(energies), occupations, np.eye(3), g, g, 1, 1)
        model = InertialFermionDynamics([channel], [[.2]], [[1.4]], [3.2])
        basis = model.tangent_basis()
        state = basis @ np.linspace(.01, .06, basis.shape[1])
        checks.close("null/"+label+"/energy", model.energy_terms(model.propagate(state, .7))["total"], model.energy_terms(state)["total"])
        checks.truth("null/"+label+"/all_EOM_entries_retained", model.state_size == 11)
        outside = model.zero_state()
        outside[2] = .1  # population diagonal in the identity orbital frame
        try:
            model.energy_terms(outside)
        except ValueError:
            checks.truth("null/"+label+"/outside_orbit_rejected", True)
        else:
            checks.truth("null/"+label+"/outside_orbit_rejected", False)
        if occupations[1] == occupations[2]:
            equal_block = model.zero_state()
            equal_block[2+1*3+2] = equal_block[2+2*3+1] = .1
            checks.truth("null/"+label+"/no_fake_oscillator", basis.shape[1] == 6)
            try:
                model.energy_terms(equal_block)
            except ValueError:
                checks.truth("null/"+label+"/equal_population_coherence_rejected", True)
            else:
                checks.truth("null/"+label+"/equal_population_coherence_rejected", False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    cpu_started = time.process_time()
    checks, details = Checks(), {}
    report = {"preregistration_commit": PREREGISTRATION, "allowance": ALLOWANCE,
        "scope": "finite original-basis held-reference unitary-orbit tangent; physical +g; independently supplied frozen curvature; no material or equilibrium-noise claim",
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "sympy": sp.__version__},
        "budget": {"wall_cap_seconds": 115, "root_total_cap_seconds": 120, "BLAS_threads": 1, "CPU_allocation": 6},
        "missing_physical_inputs": ["stationary full-target material action mapping", "coupled-equilibrium covariance/noise and physical temperature identification", "justified environment/population-relaxation or loss law"],
        "energy_scope": "conserved action Hessian on admissible real unitary tangent; not guaranteed positive; excludes equal-population occupation/coherence perturbations; no work-integral definition",
        "source_bindings": [], "details": details}
    def deadline(signum, frame):
        raise TimeoutError("finite oracle wall cap 115 seconds exceeded")
    # Existing output is never overwritten. Root owns invocation and receipts.
    with args.output.open("x") as stream:
        previous = signal.signal(signal.SIGALRM, deadline)
        signal.setitimer(signal.ITIMER_REAL, 115)
        try:
            checks.section("bindings", lambda: report.update(source_bindings=source_bindings(checks)))
            checks.section("symbolic", lambda: symbolic_checks(checks))
            checks.section("numerical", lambda: numerical_checks(checks, details))
            checks.section("null_sectors", lambda: null_sector_checks(checks))
        except BaseException as exc:
            checks.truth("fatal", False, exception=repr(exc), traceback=traceback.format_exc())
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous)
            report["wall_seconds"] = time.monotonic()-started
            report["cpu_seconds"] = time.process_time()-cpu_started
            report["checks"] = checks.records
            report["failures"] = [row for row in checks.records if not row["passed"]]
            report["passed"] = bool(checks.records) and not report["failures"]
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write("\n")
    print(json.dumps({"output": str(args.output), "passed": report["passed"],
                      "checks": len(checks.records), "failures": len(report["failures"])}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
