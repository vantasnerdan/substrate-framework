"""Covariant local kinetic diagonalization and long-wave remainder order."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0032-collective-map")
    inertia, stiffness = sp.symbols("I K", positive=True)
    micro, cage, collective, rigid = sp.symbols("Phi beta Psi rigid_angle", real=True)
    micro_rate, cage_rate, collective_rate = sp.symbols("Phi_dot beta_dot Psi_dot", real=True)
    a = 11 * inertia / 10
    b = -9 * inertia / 10
    c = 11 * inertia / 10
    forward = (a * micro + b * cage) / (a + b)
    inverse = sp.solve(sp.Eq(collective, forward), micro)[0]
    ledger.check("local physical angle map", sp.simplify(inverse - (2 * collective + 9 * cage) / 11) == 0)
    ledger.check("collective frame covariance", sp.simplify(forward.subs({micro: micro + rigid, cage: cage + rigid}) - forward - rigid) == 0)
    kinetic = (a * micro_rate**2 + 2 * b * micro_rate * cage_rate + c * cage_rate**2) / 2
    inverse_rate = inverse.subs({collective: collective_rate, cage: cage_rate})
    mapped_kinetic = sp.expand(kinetic.subs(micro_rate, inverse_rate))
    j_collective = sp.simplify(sp.diff(mapped_kinetic, collective_rate, 2))
    j_cage = sp.simplify(sp.diff(mapped_kinetic, cage_rate, 2))
    ledger.check("cross inertia eliminated exactly", sp.diff(mapped_kinetic, collective_rate, cage_rate) == 0)
    ledger.check("positive collective inertia", sp.simplify(j_collective - 2 * inertia / 55) == 0)
    ledger.check("positive retained cage inertia", sp.simplify(j_cage - 4 * inertia / 11) == 0)
    mapped_potential = sp.expand((stiffness * (micro - cage)**2 / 2).subs(micro, inverse))
    k_collective = sp.simplify(sp.diff(mapped_potential, collective, 2))
    ledger.check("physical relative-angle potential retained", sp.simplify(mapped_potential - k_collective * (collective - cage)**2 / 2) == 0)
    ledger.check("mapped stiffness", sp.simplify(k_collective - 4 * stiffness / 121) == 0)
    ledger.check("constrained cage frequency differs from free pair", sp.simplify(k_collective / j_collective - 10 * stiffness / (11 * inertia)) == 0)
    ledger.mutation_sensitive(
        "transformed rather than borrowed inertia",
        lambda value: sp.simplify(value - sp.diff(mapped_kinetic, collective_rate, 2)) == 0,
        j_collective,
        [inertia],
    )

    tension, radius_squared = sp.symbols("T r_squared", positive=True)
    micro_gradient, cage_gradient, collective_gradient = sp.symbols("Phi_s beta_s Psi_s", real=True)
    gradient_energy = 3 * tension * radius_squared * (micro_gradient**2 + cage_gradient**2) / 2
    mapped_gradient = sp.expand(gradient_energy.subs(micro_gradient, (2 * collective_gradient + 9 * cage_gradient) / 11))
    collective_twist = sp.diff(mapped_gradient, collective_gradient, 2)
    ledger.check("mapped physical helical twist", sp.simplify(collective_twist - 12 * tension * radius_squared / 121) == 0)
    ledger.check("gradient cross term retained in exact action", sp.diff(mapped_gradient, collective_gradient, cage_gradient) != 0)

    wave_number, frequency_squared, slope = sp.symbols("k z slope", real=True)
    bulk_mass, spin_mass, spring, spin_gradient, cage_mass = sp.symbols("rho j K0 C0 b0", positive=True)
    # A mass correction b0*k^2 changes the coupled transverse determinant by:
    determinant_change = cage_mass * wave_number**2 * frequency_squared * (spin_mass * frequency_squared - spring - spin_gradient * wave_number**2)
    acoustic_change = sp.expand(determinant_change.subs(frequency_squared, slope * wave_number**2))
    optical_change = sp.expand(determinant_change.subs(frequency_squared, spring / spin_mass + slope * wave_number**2))
    ledger.check("acoustic correction starts at k^4", acoustic_change.coeff(wave_number, 2) == 0 and acoustic_change.coeff(wave_number, 4) != 0)
    ledger.check("optical correction starts at k^4", optical_change.coeff(wave_number, 2) == 0 and optical_change.coeff(wave_number, 4) != 0)
    determinant_at_zero = bulk_mass * frequency_squared * (spin_mass * frequency_squared - spring)
    ledger.check("acoustic root simple in omega squared", sp.diff(determinant_at_zero, frequency_squared).subs(frequency_squared, 0) != 0)
    ledger.check("optical root simple in omega squared", sp.diff(determinant_at_zero, frequency_squared).subs(frequency_squared, spring / spin_mass) != 0)
    coupling, cross_gradient = sp.symbols("coupling cross_gradient", real=True)
    coupling_change = sp.expand((coupling * wave_number + cross_gradient * wave_number**3)**2 - coupling**2 * wave_number**2)
    ledger.check("cross-gradient correction starts at k^4", coupling_change.coeff(wave_number, 2) == 0 and coupling_change.coeff(wave_number, 4) != 0)
    print("Psi =", forward)
    print("Phi =", inverse)
    print("J_Psi =", j_collective)
    print("J_cage =", j_cage)
    print("K_Psi =", k_collective)
    print("C_Psi_log_leading =", collective_twist)
    print("constrained_cage_gap =", sp.simplify(k_collective / j_collective))
    print("Scope: exact local collective map; standard micropolar dispersion through k^2")
    print("requires an affine material cage, a bulk kinetic closure, and finite-core limits.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
