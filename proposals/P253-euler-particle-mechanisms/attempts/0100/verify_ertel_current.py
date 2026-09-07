#!/usr/bin/env python3
"""Exact checks for the P253/0100 Ertel-current supplier."""

import sympy as sp

from substrate_framework.euler_ertel_current import (
    cao_azimuthal_phase_charge,
    cao_azimuthal_phase_density,
    cao_forced_magnetization_flux,
    cao_ertel_density,
    closed_vorticity_line_multiplier,
    ertel_charge_from_flux,
    ertel_current_parity,
    ertel_stretching_residual,
    forced_ertel_flux,
    forced_ertel_source,
    forced_lock_residual,
    geometric_azimuth_material_derivative,
    transported_lock_residual,
)


def main() -> None:
    entries = sp.symbols("a0:9")
    grad_u = sp.Matrix(3, 3, entries)
    omega = sp.Matrix(sp.symbols("w0:3"))
    grad_chi = sp.Matrix(sp.symbols("c0:3"))
    assert sp.simplify(ertel_stretching_residual(omega, grad_chi, grad_u)) == 0
    print("PASS 1: Euler stretching and transported-gradient terms cancel exactly")

    x, y, z = sp.symbols("x y z", real=True)
    chi = sp.Function("chi")(x, y, z)
    w1 = sp.Function("w1")(x, y, z)
    w2 = sp.Function("w2")(x, y, z)
    w3 = sp.Function("w3")(x, y, z)
    div_chi_omega = sp.diff(chi * w1, x) + sp.diff(chi * w2, y) + sp.diff(chi * w3, z)
    q = w1 * sp.diff(chi, x) + w2 * sp.diff(chi, y) + w3 * sp.diff(chi, z)
    div_omega = sp.diff(w1, x) + sp.diff(w2, y) + sp.diff(w3, z)
    assert sp.simplify(div_chi_omega - q - chi * div_omega) == 0
    assert ertel_charge_from_flux(0) == 0
    print("PASS 2: q_E=div(chi omega) for div-free vorticity and zero flux gives zero total charge")

    f1 = sp.Function("f1")(x, y, z)
    f2 = sp.Function("f2")(x, y, z)
    f3 = sp.Function("f3")(x, y, z)
    force = sp.Matrix([f1, f2, f3])
    grad_chi = sp.Matrix([sp.diff(chi, x), sp.diff(chi, y), sp.diff(chi, z)])
    curl_f = sp.Matrix([
        sp.diff(f3, y) - sp.diff(f2, z),
        sp.diff(f1, z) - sp.diff(f3, x),
        sp.diff(f2, x) - sp.diff(f1, y),
    ])
    cross = force.cross(grad_chi)
    div_cross = sp.diff(cross[0], x) + sp.diff(cross[1], y) + sp.diff(cross[2], z)
    assert sp.simplify(div_cross - forced_ertel_source(curl_f, grad_chi)) == 0
    uvec = sp.Matrix(sp.symbols("u0:3"))
    qsym = sp.Symbol("q")
    assert forced_ertel_flux(qsym, uvec, force, grad_chi) == qsym * uvec - cross
    print("PASS 3: forced Ertel source is a divergence and fixes the conserved flux sign")

    assert ertel_current_parity("scalar") == {
        "density": "pseudoscalar",
        "spatial_current": "axial",
    }
    assert ertel_current_parity("pseudoscalar") == {
        "density": "scalar",
        "spatial_current": "polar",
    }
    print("PASS 4: true-scalar and pseudoscalar tag parities are distinct")

    tag, dt_lambda = sp.symbols("tag dt_lambda", nonzero=True)
    assert transported_lock_residual(tag, dt_lambda) == -tag * dt_lambda
    assert transported_lock_residual(tag, 0) == 0
    print("PASS 5: preservation of q_E=lambda chi requires advected lambda on nonzero tag support")

    forced_residual = forced_lock_residual(tag, dt_lambda, curl_f, grad_chi)
    assert forced_residual == tag * dt_lambda - curl_f.dot(grad_chi)
    assert sp.simplify(
        forced_residual.subs(dt_lambda, curl_f.dot(grad_chi) / tag)
    ) == 0
    print("PASS 6: forced lock obeys chi*D_t lambda=(curl f) dot grad chi")

    lam, period = sp.symbols("lam period", real=True)
    assert closed_vorticity_line_multiplier(lam * period) == sp.exp(lam * period)
    assert closed_vorticity_line_multiplier(0) == 1
    print("PASS 7: closed-line lock has the exact exponential monodromy")

    zeta = sp.Symbol("zeta", nonzero=True)
    assert cao_ertel_density(zeta, 0) == 0
    print("PASS 8: the axisymmetric Cao tag is orthogonal to toroidal vorticity")

    g, rho, tag_value, tag_prime = sp.symbols("g rho tag tagprime", nonzero=True)
    grad_phi = sp.Matrix(sp.symbols("gp0:3"))
    grad_P = sp.Matrix(sp.symbols("gP0:3"))
    force_cao = g * tag_value / rho * (-grad_phi - sp.Symbol("H") * grad_P)
    grad_tag = tag_prime * grad_P
    flux_direct = -force_cao.cross(grad_tag)
    flux_closed = cao_forced_magnetization_flux(
        g, rho, tag_value, tag_prime, grad_phi, grad_P
    )
    assert sp.simplify(flux_direct - flux_closed) == sp.zeros(3, 1)
    chi2 = sp.Function("chi2")(x, y, z)
    phi = sp.Function("phi")(x, y, z)
    grad_chi2 = sp.Matrix([sp.diff(chi2, x), sp.diff(chi2, y), sp.diff(chi2, z)])
    grad_phi_fun = sp.Matrix([sp.diff(phi, x), sp.diff(phi, y), sp.diff(phi, z)])
    magnetization = phi * grad_chi2
    curl_m = sp.Matrix([
        sp.diff(magnetization[2], y) - sp.diff(magnetization[1], z),
        sp.diff(magnetization[0], z) - sp.diff(magnetization[2], x),
        sp.diff(magnetization[1], x) - sp.diff(magnetization[0], y),
    ])
    assert sp.simplify(curl_m - grad_phi_fun.cross(grad_chi2)) == sp.zeros(3, 1)
    div_curl_m = sp.diff(curl_m[0], x) + sp.diff(curl_m[1], y) + sp.diff(curl_m[2], z)
    assert sp.simplify(div_curl_m) == 0
    print("PASS 9: charged Cao correction is an exact divergence-free magnetization superpotential")

    kappa = sp.Symbol("kappa", nonzero=True)
    assert cao_azimuthal_phase_density(zeta) == zeta
    assert cao_azimuthal_phase_charge(kappa) == 2 * sp.pi * kappa
    print("PASS 10: the azimuthal circle phase gives q_E=zeta and total charge 2*pi*kappa")

    u_theta, radius = sp.symbols("u_theta radius", nonzero=True)
    assert geometric_azimuth_material_derivative(u_theta, radius) == u_theta / radius
    assert geometric_azimuth_material_derivative(0, radius) == 0
    print("PASS 11: fixed geometric azimuth is material only on the no-swirl subspace")


if __name__ == "__main__":
    main()
