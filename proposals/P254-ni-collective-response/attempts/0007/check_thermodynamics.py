"""Exact Fock Kubo-Mori trace plus independent expm free-energy Hessians."""
import numpy as np
import sympy as sp

from substrate_framework.finite_fermion_state import fermi_dirac, fock_gibbs, second_quantize
from substrate_framework.finite_fermion_response import isolated_static, noise_lines
from substrate_framework.finite_fermion_thermodynamics import grandcanonical_curvature, mean_number_curvature
from fixtures import OPERATORS, ENERGIES, BETA, zero, serial
from free_energy_oracle import free_energy_functions, negative_hessian

def matrix_text(matrix):
    """Preserve readable exact sums of logarithms without log-combining huge integers."""
    return [[str(sp.expand_log(matrix[i, j], force=True)) for j in range(matrix.cols)]
            for i in range(matrix.rows)]



def fock_isothermal(rho, energies, lifted):
    """Full state divided differences from derivative of exp(-H), centered zero block."""
    centered = [a-sp.trace(rho*a)*sp.eye(rho.rows) for a in lifted]
    result = sp.zeros(len(lifted))
    for n, en in enumerate(energies):
        for m, em in enumerate(energies):
            gap = sp.simplify(em-en)
            coefficient = BETA*rho[n, n] if gap == 0 else (rho[n, n]-rho[m, m])/gap
            for a in range(len(lifted)):
                for b in range(len(lifted)):
                    result[a, b] += coefficient*centered[a][n, m]*centered[b][m, n]
    return result.applyfunc(sp.simplify)


def run(ledger):
    lifted = tuple(second_quantize(a) for a in OPERATORS)
    outputs = []
    for energies in ENERGIES:
        f = fermi_dirac(energies, BETA)
        rho, _, many = fock_gibbs(energies, BETA)
        static = isolated_static(energies, f, OPERATORS)
        gc = grandcanonical_curvature(energies, f, OPERATORS, BETA)
        direct = fock_isothermal(rho, many, lifted)
        ledger.check('isothermal GC curvature equals full Fock divided difference', zero(gc-direct))
        populations = noise_lines(energies, f, OPERATORS)[0]
        ledger.check('GC minus isolated is exactly zero-frequency occupation covariance',
                     zero(gc-static-BETA*populations))
        ledger.check('equating isolated zero response and isothermal GC rejected', not zero(gc-static))
        legendre = mean_number_curvature(energies, f, OPERATORS, BETA)
        canonical_rho, canonical_z, _ = fock_gibbs(energies, BETA, number=1)
        canonical = fock_isothermal(canonical_rho, many, lifted)
        ledger.check('finite canonical differs from mean-number Legendre curvature', not zero(canonical-legendre))
        ledger.check('mean-number projection changes GC curvature', not zero(legendre-gc))
        number = second_quantize(sp.eye(3))
        ledger.check('canonical number variance is exactly zero', zero(
            sp.trace(canonical_rho*number**2)-sp.trace(canonical_rho*number)**2))
        n0 = second_quantize(sp.diag(1, 0, 0))
        n1 = second_quantize(sp.diag(0, 1, 0))
        canonical_cov = sp.simplify(sp.trace(canonical_rho*n0*n1)
                                    -sp.trace(canonical_rho*n0)*sp.trace(canonical_rho*n1))
        ledger.check('canonical factorized independent occupations mutation rejected', canonical_cov != 0)
        functions = free_energy_functions([float(e) for e in energies])
        rows = []
        for label, function, exact in zip(('GC', 'mean_N_Legendre', 'canonical_N1'),
                                          functions[:3], (gc, legendre, canonical)):
            ledger.check(f'{label} curvature is exactly real symmetric',
                         zero(exact-exact.conjugate()) and zero(exact-exact.T))
            numerical = []
            errors = []
            for step in (0.002, 0.001, 0.0005):
                matrix = negative_hessian(function, step)
                error = float(np.max(np.abs(matrix-np.array(exact.evalf(), complex))))
                numerical.append({'step': step, 'curvature': matrix.tolist(), 'max_abs_error': error})
                errors.append(error)
            ledger.check(f'{label} independent expm free-energy derivative refinement',
                         errors[-1] < 0.0002 and errors[-1] < errors[0])
            rows.append({'ensemble': label, 'exact': matrix_text(exact), 'finite_differences': numerical})
        outputs.append({'energies': energies, 'mean_N_at_mu0': sum(f),
                        'isolated_chi_zero': matrix_text(static), 'zero_covariance': matrix_text(populations),
                        'grandcanonical_curvature': matrix_text(gc),
                        'mean_number_legendre_curvature': matrix_text(legendre),
                        'canonical_N1_curvature': matrix_text(canonical), 'canonical_partition': canonical_z,
                        'canonical_n0_n1_covariance': canonical_cov,
                        'Legendre_minus_canonical': matrix_text(legendre-canonical),
                        'independent_free_energy_derivatives': rows})
    return serial(outputs)
