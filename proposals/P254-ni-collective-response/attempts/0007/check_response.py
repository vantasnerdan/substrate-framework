"""Full Fock Lehmann traces independently expose every complex two-point entry."""
import numpy as np
import sympy as sp
from scipy.integrate import quad

from substrate_framework.finite_fermion_state import fermi_dirac, fock_gibbs, second_quantize
from substrate_framework.finite_fermion_response import retarded, noise_lines, symmetrized_lines
from fixtures import OPERATORS, ENERGIES, BETA, HBAR, Z, zero, serial


def direct_lines(rho, many_energies, lifted):
    """Insert full Fock eigenstate projectors in trace rho delta A(t) delta B."""
    centered = [op-sp.trace(rho*op)*sp.eye(rho.rows) for op in lifted]
    lines = {}
    for n, en in enumerate(many_energies):
        for m, em in enumerate(many_energies):
            if any(centered[a][n, m] != 0 for a in range(len(lifted))):
                gap = sp.expand_log(sp.simplify(em-en), force=True)
                weight = sp.Matrix(len(lifted), len(lifted), lambda a, b:
                                   rho[n, n]*centered[a][n, m]*centered[b][m, n])
                lines[gap] = lines.get(gap, sp.zeros(len(lifted)))+weight
    return lines


def direct_kernel(rho, energies, lifted, z):
    # Full-state commutator, not FD/Wick reduction.
    return sp.Matrix(len(lifted), len(lifted), lambda a, b: sum(
        (rho[n, n]-rho[m, m])*lifted[a][n, m]*lifted[b][m, n]
        /(energies[m]-energies[n]-HBAR*z)
        for n in range(rho.rows) for m in range(rho.rows)
        if energies[m] != energies[n] and lifted[a][n, m]*lifted[b][m, n] != 0))


def integral_check(rho, energies, lifted):
    """Numerically integrate +i/hbar actual Fock commutator, not spectral kernel."""
    e = np.array([float(x) for x in energies])
    r = np.array(rho, complex)
    a, b = [np.array(x, complex) for x in lifted]
    z = complex(Z)
    def value(t):
        evolution = np.exp(1j*e*t/float(HBAR))
        at = evolution[:, None]*a*evolution.conj()[None, :]
        return np.exp(1j*z*t)*1j/float(HBAR)*np.trace(r@(at@b-b@at))
    real, er = quad(lambda t: value(t).real, 0, np.inf, epsabs=1e-10, limit=300)
    imag, ei = quad(lambda t: value(t).imag, 0, np.inf, epsabs=1e-10, limit=300)
    return real+1j*imag, er+ei


def run(ledger):
    lifted = tuple(second_quantize(a) for a in OPERATORS)
    result = []
    for energies in ENERGIES:
        rho, _, many = fock_gibbs(energies, BETA)
        f = fermi_dirac(energies, BETA)
        kernel = retarded(energies, f, OPERATORS, Z, HBAR)
        direct = direct_kernel(rho, many, lifted, Z)
        ledger.check('complete complex retarded kernel equals independent Fock trace', zero(kernel-direct))
        ledger.check('retarded sign mutation rejected', not zero(-kernel-direct))
        wrong_conjugation = retarded(energies, f, (OPERATORS[0].T, OPERATORS[1]), Z, HBAR)
        ledger.check('complex conjugation mutation rejected', not zero(wrong_conjugation-direct))
        wrong_hbar = retarded(energies, f, OPERATORS, Z, 1)
        ledger.check('hbar normalization mutation rejected', not zero(wrong_hbar-direct))
        integral, error = integral_check(rho, many, lifted)
        residual = abs(integral-complex(sp.N(kernel[0, 1])))
        ledger.check('direct time commutator integral fixes perturbation sign and hbar', residual < 1e-8)
        lines = noise_lines(energies, f, OPERATORS)
        full = direct_lines(rho, many, lifted)
        ledger.check('all connected unsymmetrized line weights equal exact Fock trace',
                     set(lines) == set(full) and all(zero(lines[g]-full[g]) for g in lines))
        ledger.check('all line weights obey complex detailed balance', all(
            zero(lines[-g].T-sp.exp(-BETA*g)*w) for g, w in lines.items()))
        sym = symmetrized_lines(lines)
        ledger.check('symmetrized noise equals independent Fock anticommutator', all(
            zero(sym[g]-(full[g]+full[-g].T)/2) for g in lines))
        ledger.check('FDT including complex cross channels at every nonzero line', all(
            zero(sym[g]-(w-lines[-g].T)*sp.coth(BETA*g/2)/2)
            for g, w in lines.items() if g != 0))
        no_holes = {g: w.copy() for g, w in lines.items()}
        i, j = 0, 2
        g = energies[j]-energies[i]
        extra = sp.Matrix(2, 2, lambda a, b: f[i]*f[j]*OPERATORS[a][i, j]*OPERATORS[b][j, i])
        no_holes[g] += extra
        ledger.check('missing Pauli hole occupation factor rejected', not zero(no_holes[g]-full[g]))
        ledger.check('dropping zero-frequency population covariance rejected', not zero(full[0]))
        diagonal = sp.Matrix(2, 2, lambda a, b: sum(
            f[i]*(1-f[i])*OPERATORS[a][i, i]*OPERATORS[b][i, i] for i in range(3)))
        if energies[0] == energies[1]:
            ledger.check('dropping degenerate off-diagonal population covariance rejected',
                         not zero(diagonal-full[0]))
        result.append({'energies': energies, 'chi_at_z': kernel, 'noise_energy_lines': lines,
                       'symmetrized_energy_lines': sym, 'diagonal_only_zero_noise': diagonal,
                       'direct_integral_chi_AB': str(integral), 'integral_absolute_error': residual,
                       'quadrature_error_estimate': error, 'missing_hole_residual': extra})
    return serial(result)
