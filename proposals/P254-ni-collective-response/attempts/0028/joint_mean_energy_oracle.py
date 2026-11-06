"""Independent exact canonical mean-energy oracle; no producer/consumer imports.

Grounding: streamed_joint_reduction.JointReduction.stage,
material_joint_field_scalar.coordinate_correction, and
physical_field_environment.UniformMaxwellField.energy. This checks only
baseline mean energy, not quantum fluctuations, heat, or excitation energy.
"""

from types import SimpleNamespace

import numpy as np
import sympy as sp


def _vector(name, size, *, positive=False):
    assumptions = {"positive": True} if positive else {"real": True, "nonzero": True}
    return np.array(sp.symbols(f"{name}0:{size}", **assumptions), dtype=object)


def _symmetric(name, size):
    result = np.empty((size, size), dtype=object)
    for i in range(size):
        for j in range(i, size):
            result[i, j] = result[j, i] = sp.Symbol(
                f"{name}{i}_{j}", real=True, nonzero=True)
    return result


def _dot(left, right):
    return sum((x * y for x, y in zip(left, right)), sp.S.Zero)


def _quadratic(matrix, vector):
    """Half a symmetric quadratic form, independently enumerating monomials."""
    diagonal = sum(matrix[i, i] * vector[i] ** 2 / 2 for i in range(len(vector)))
    cross = sum(matrix[i, j] * vector[i] * vector[j]
                for i in range(len(vector)) for j in range(i + 1, len(vector)))
    return diagonal + cross


def _bilinear(left, matrix, right):
    return sum(left[i] * matrix[i, j] * right[j]
               for i in range(len(left)) for j in range(len(right)))


def _exact(expression):
    return sp.cancel(sp.expand(sp.sympify(expression)))


def _identity(actual, expected):
    residual = _exact(actual - expected)
    return {"passed": bool(residual == 0), "residual": str(residual)}


def _identities(actual, expected):
    entries = [_identity(a, b) for a, b in zip(actual, expected)]
    return {"passed": all(item["passed"] for item in entries), "entries": entries}


def _fixture():
    """Two ions, each with all three Cartesian coordinates; no zero couplings."""
    n = sp.Symbol("N", positive=True, integer=True)
    measure, field_mass = sp.symbols("dV M_field", positive=True)
    u, p, a, pi = (_vector(name, size) for name, size in
                   (("u", 6), ("p", 6), ("a", 3), ("Pi", 3)))
    masses, charges = _vector("m", 2, positive=True), _vector("Z", 2, positive=True)
    reference = SimpleNamespace(
        nbar=_vector("nbar", 4).reshape(2, 2),
        ionic_force=_vector("F", 6), support_force=_vector("S", 6),
        curvature=_symmetric("C", 6), masses=np.repeat(masses, 3), ncells=n)
    coefficients = SimpleNamespace(
        delta_curvature=_symmetric("DeltaC", 6), contact=_symmetric("D", 3),
        mixed=np.array([[sp.Symbol(f"B{i}_{j}", real=True, nonzero=True)
                         for j in range(3)] for i in range(6)], dtype=object))
    stage = SimpleNamespace(
        density=_vector("rho", 4).reshape(2, 2), gradient=_vector("g", 6),
        field_gradient_electronic=_vector("ell", 3),
        potential=_vector("phi", 4).reshape(2, 2), u=u, p=p, shift=a, momentum=pi)
    return dict(reference=reference, coefficients=coefficients,
                field=SimpleNamespace(mass=field_mass), charges=charges, stage=stage,
                one_body=sp.Symbol("E0", real=True, nonzero=True), grid_measure=measure)


def _canonical(inputs):
    """Compose electronic, ionic, coupling and global field Hamiltonian blocks."""
    r, c, s = (inputs[key] for key in ("reference", "coefficients", "stage"))
    n, u, p, a = r.ncells, s.u, s.p, s.shift
    z = [inputs["charges"][i // 3] for i in range(6)]
    # Canonical momentum minus charge times Cartesian vector potential.
    mechanical = [p[i] - z[i] * a[i % 3] for i in range(6)]
    kinetic = sum(mechanical[i] ** 2 / (2 * r.masses[i]) for i in range(6))
    hartree = inputs["grid_measure"] * _dot(
        (s.density - r.nbar).flat, s.potential.flat) / 2
    electronic = inputs["one_body"] + hartree + _dot(u, s.gradient)
    ionic = -_dot(u, r.ionic_force + r.support_force)
    ionic += _quadratic(r.curvature + c.delta_curvature, u) + kinetic
    coupling = (_dot(a, s.field_gradient_electronic) + _quadratic(c.contact, a)
                + _bilinear(u, c.mixed, a))
    maxwell = _dot(s.momentum, s.momentum) / (2 * inputs["field"].mass)
    matter = electronic + ionic + coupling
    hamiltonian = n * matter + maxwell
    # Physical decomposition is checked separately from the composed Hamiltonian.
    terms = dict(
        bare_electronic=n * inputs["one_body"], hartree=n * hartree,
        electron_ionic=n * _dot(u, s.gradient),
        ionic_reference_linear=-n * _dot(u, r.ionic_force + r.support_force),
        ionic_C=n * _quadratic(r.curvature, u), ionic_kinetic_minimal=n * kinetic,
        delta_C=n * _quadratic(c.delta_curvature, u),
        electron_field=n * _dot(a, s.field_gradient_electronic),
        field_contact=n * _quadratic(c.contact, a),
        ion_field_mixed=n * _bilinear(u, c.mixed, a), maxwell=maxwell)
    return hamiltonian, matter, terms, mechanical, z


def _comparison(actual, expected, hamiltonian):
    checks = {name: _identity(actual[name], value) for name, value in expected.items()}
    checks["full_hamiltonian"] = _identity(sum(actual.values()), hamiltonian)
    return checks


def check_components(component_function):
    """Call the supplied implementation on exact arrays and return JSON-safe proof.

    Every identity and mutation is symbolic; no numerical sample banks a check.
    Labels select physical components, not an incidental spelling/length test.
    """
    inputs = _fixture()
    hamiltonian, matter, expected, mechanical, charges = _canonical(inputs)
    actual = {name: sp.sympify(value)
              for name, value in component_function(**inputs).items()}
    checks = _comparison(actual, expected, hamiltonian)
    total = sum(actual.values())
    r, c, s = (inputs[key] for key in ("reference", "coefficients", "stage"))
    n, u, a, pi = r.ncells, s.u, s.shift, s.momentum
    # Hamilton equations independently expose signs, canonical kinetic counting,
    # Cartesian charge replication, and the absence of N on the global field.
    checks["canonical_ionic_momentum_derivative"] = _identities(
        [sp.diff(total, p) for p in s.p],
        [n * mechanical[i] / r.masses[i] for i in range(6)])
    checks["ionic_potential_gradient"] = _identities(
        [sp.diff(total, q) for q in u],
        [n * (s.gradient[i] - r.ionic_force[i] - r.support_force[i]
              + sum((r.curvature[i, j] + c.delta_curvature[i, j]) * u[j]
                    for j in range(6))
              + sum(c.mixed[i, j] * a[j] for j in range(3))) for i in range(6)])
    checks["cartesian_field_gradient"] = _identities(
        [sp.diff(total, q) for q in a],
        [n * (s.field_gradient_electronic[j]
              + sum(c.contact[j, k] * a[k] for k in range(3))
              + sum(c.mixed[i, j] * u[i] for i in range(6))
              - sum(charges[i] * mechanical[i] / r.masses[i]
                    for i in range(j, 6, 3))) for j in range(3)])
    checks["global_maxwell_velocity"] = _identities(
        [sp.diff(total, q) for q in pi], [q / inputs["field"].mass for q in pi])
    checks["explicit_density_factor_at_fixed_saved_potential"] = _identities(
        [sp.diff(total, q) for q in s.density.flat],
        [n * inputs["grid_measure"] * q / 2 for q in s.potential.flat])
    checks["explicit_reference_factor_at_fixed_saved_potential"] = _identities(
        [sp.diff(total, q) for q in r.nbar.flat],
        [-n * inputs["grid_measure"] * q / 2 for q in s.potential.flat])
    checks["cell_extensivity"] = _identity(sp.diff(total, n), matter)
    checks["maxwell_has_no_extra_N_at_fixed_global_mass"] = _identity(
        sp.diff(actual["maxwell"], n), 0)
    baseline_passed = all(item["passed"] for item in checks.values())

    neutral = n * sum(s.p[i] ** 2 / (2 * r.masses[i]) for i in range(6))
    # Each physical defect has an exact mixed/second derivative witness.  Taking
    # derivatives of the mutation delta prevents unrelated baseline errors from
    # masquerading as sensitivity to the intended error.
    specifications = [
        ("missing_matter_N", {key: value / n for key, value in actual.items()
                              if key != "maxwell"}, [inputs["one_body"]], 1 - n),
        ("double_neutral_kinetic", {"ionic_kinetic_minimal":
          actual["ionic_kinetic_minimal"] + neutral}, [s.p[0], s.p[0]], n / r.masses[0]),
        ("mixed_sign", {"ion_field_mixed": -actual["ion_field_mixed"]},
         [u[0], a[0]], -2 * n * c.mixed[0, 0]),
        ("contact_sign", {"field_contact": -actual["field_contact"]},
         [a[0], a[0]], -2 * n * c.contact[0, 0]),
        ("support_sign", {"ionic_reference_linear": actual["ionic_reference_linear"]
          + 2 * n * _dot(u, r.support_force)}, [u[0], r.support_force[0]], 2 * n),
        ("extra_DeltaC", {"delta_C": actual["delta_C"]
          + n * _quadratic(c.delta_curvature, u)},
         [u[0], u[0], c.delta_curvature[0, 0]], n),
        ("spurious_N_on_global_maxwell", {"maxwell": n * actual["maxwell"]},
         [pi[0], pi[0]], (n - 1) / inputs["field"].mass),
    ]
    mutations = {}
    for name, replacements, variables, sensitivity in specifications:
        mutant = dict(actual)
        mutant.update(replacements)
        comparison = _comparison(mutant, expected, hamiltonian)
        delta = _exact(sum(mutant.values()) - total)
        witness = _exact(sp.diff(delta, *variables))
        witness_check = _identity(witness, sensitivity)
        failed = {key: value for key, value in comparison.items() if not value["passed"]}
        detected = (baseline_passed and not comparison["full_hamiltonian"]["passed"]
                    and witness_check["passed"] and witness != 0)
        mutations[name] = dict(
            detected=bool(detected), changed_components=list(replacements),
            failed_checks=failed, exact_defect=str(delta),
            sensitivity=dict(variables=[str(q) for q in variables],
                             derivative=str(witness), expected=str(sensitivity),
                             nonzero=bool(witness != 0), identity=witness_check))
    passed = baseline_passed and all(item["detected"] for item in mutations.values())
    return dict(status="passed" if passed else "failed", arithmetic="exact_sympy",
                checks=checks, mutations=mutations,
                derivative_scope="Density-factor partials hold the saved potential fixed, not the full K(delta n) variation. Partial N holds the supplied global field mass fixed; physical mass varies with total volume. Mean ionic velocity uses the existing 1/N Poisson bracket.")
