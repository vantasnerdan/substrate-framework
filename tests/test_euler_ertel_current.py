import sympy as sp

from substrate_framework.euler_ertel_current import (
    cao_azimuthal_phase_charge,
    cao_azimuthal_phase_density,
    cao_ertel_density,
    closed_vorticity_line_multiplier,
    cao_forced_magnetization_flux,
    ertel_charge_from_flux,
    ertel_current_parity,
    ertel_stretching_residual,
    forced_ertel_flux,
    forced_ertel_source,
    forced_lock_residual,
    geometric_azimuth_material_derivative,
    transported_lock_residual,
)


def test_ertel_stretching_cancels_for_arbitrary_velocity_gradient():
    entries = sp.symbols("a0:9")
    grad_u = sp.Matrix(3, 3, entries)
    omega = sp.Matrix(sp.symbols("w0:3"))
    grad_chi = sp.Matrix(sp.symbols("c0:3"))
    assert sp.simplify(ertel_stretching_residual(omega, grad_chi, grad_u)) == 0


def test_true_scalar_ertel_current_is_axial():
    assert ertel_current_parity("scalar") == {
        "density": "pseudoscalar",
        "spatial_current": "axial",
    }
    assert ertel_current_parity("pseudoscalar") == {
        "density": "scalar",
        "spatial_current": "polar",
    }


def test_ertel_charge_is_boundary_flux_and_vanishes_when_flux_does():
    flux = sp.Symbol("Phi_boundary", real=True)
    assert ertel_charge_from_flux(flux) == flux
    assert ertel_charge_from_flux(0) == 0


def test_forced_ertel_source_and_conserved_flux():
    curl_f = sp.Matrix(sp.symbols("cf0:3"))
    grad_chi = sp.Matrix(sp.symbols("gc0:3"))
    q = sp.Symbol("q")
    u = sp.Matrix(sp.symbols("u0:3"))
    force = sp.Matrix(sp.symbols("f0:3"))
    assert forced_ertel_source(curl_f, grad_chi) == curl_f.dot(grad_chi)
    assert forced_ertel_flux(q, u, force, grad_chi) == q * u - force.cross(grad_chi)


def test_lock_requires_advected_lambda_where_tag_is_nonzero():
    chi, dt_lambda = sp.symbols("chi dt_lambda", nonzero=True)
    assert transported_lock_residual(chi, dt_lambda) == -chi * dt_lambda
    assert transported_lock_residual(chi, 0) == 0


def test_forced_lock_has_exact_source_condition():
    chi, dt_lambda = sp.symbols("chi dt_lambda", nonzero=True)
    curl_f = sp.Matrix(sp.symbols("cf0:3"))
    grad_chi = sp.Matrix(sp.symbols("gc0:3"))
    assert forced_lock_residual(chi, dt_lambda, curl_f, grad_chi) == (
        chi * dt_lambda - curl_f.dot(grad_chi)
    )


def test_closed_line_monodromy_constant_and_zero_integral_cases():
    lam, period = sp.symbols("lam period", real=True)
    assert closed_vorticity_line_multiplier(lam * period) == sp.exp(lam * period)
    assert closed_vorticity_line_multiplier(0) == 1


def test_axisymmetric_cao_tag_has_zero_ertel_density():
    zeta = sp.Symbol("zeta", nonzero=True)
    assert cao_ertel_density(zeta, 0) == 0


def test_cao_azimuthal_circle_phase_has_exact_density_and_charge():
    zeta, kappa = sp.symbols("zeta kappa", real=True)
    assert cao_azimuthal_phase_density(zeta) == zeta
    assert cao_azimuthal_phase_charge(kappa) == 2 * sp.pi * kappa


def test_fixed_geometric_azimuth_is_material_only_without_swirl():
    u_theta, radius = sp.symbols("u_theta radius", nonzero=True)
    assert geometric_azimuth_material_derivative(u_theta, radius) == u_theta / radius
    assert geometric_azimuth_material_derivative(0, radius) == 0


def test_charged_cao_force_gives_magnetization_flux():
    g, rho, chi, chip = sp.symbols("g rho chi chip", nonzero=True)
    grad_phi = sp.Matrix(sp.symbols("p0:3"))
    grad_P = sp.Matrix(sp.symbols("P0:3"))
    got = cao_forced_magnetization_flux(g, rho, chi, chip, grad_phi, grad_P)
    want = g * chi * chip / rho * grad_phi.cross(grad_P)
    assert sp.simplify(got - want) == sp.zeros(3, 1)


def test_untyped_tag_parity_is_rejected():
    try:
        ertel_current_parity("density")
    except ValueError:
        pass
    else:
        raise AssertionError("untyped tag parity must be rejected")
