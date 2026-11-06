"""Expose CAR, grandcanonical Gibbs and non-Gaussian finite population behavior."""
import sympy as sp

from substrate_framework.finite_fermion_state import annihilators, fermi_dirac, fock_gibbs, second_quantize
from fixtures import A, B, BETA, ENERGIES, zero, serial


def run(ledger):
    c = annihilators(3)
    identity = sp.eye(8)
    ledger.check('fixture operators complex Hermitian and noncommuting',
                 A == A.H and B == B.H and not zero(A*B-B*A)
                 and any(sp.im(v) != 0 for v in A))
    ledger.check('explicit bit Fock CAR for all orbital pairs', all(
        zero(c[i]*c[j].H+c[j].H*c[i]-(identity if i == j else sp.zeros(8)))
        and zero(c[i]*c[j]+c[j]*c[i]) for i in range(3) for j in range(3)))
    bad = c[1].as_mutable()
    bad[1, 3] *= -1
    ledger.check('fermion sign mutation exposes CAR violation',
                 not zero(c[0]*bad+bad*c[0]))
    outputs = []
    for energies in ENERGIES:
        rho, partition, _ = fock_gibbs(energies, BETA)
        f = fermi_dirac(energies, BETA)
        ledger.check('FD partition equals full Fock enumeration', zero(
            partition-sp.prod(1+sp.exp(-e) for e in energies)))
        one_body = sp.Matrix(3, 3, lambda i, j: sp.trace(rho*c[i].H*c[j]))
        ledger.check('FD density equals exact Fock one-body trace', one_body == sp.diag(*f))
        ledger.check('physical Gibbs positive normalized state',
                     sp.trace(rho) == 1 and all(p.is_positive for p in rho.diagonal()))
        number = c[2].H*c[2]
        center = number-f[2]*identity
        fourth_cumulant = sp.simplify(sp.trace(rho*center**4)
                                      -3*sp.trace(rho*center**2)**2)
        ledger.check('population fourth cumulant differs from Gaussian zero', fourth_cumulant != 0)
        outputs.append({'energies': energies, 'partition': partition, 'f': f,
                        'rho_diagonal': list(rho.diagonal()), 'fourth_cumulant_n2': fourth_cumulant})
    alternate = (sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 10))
    odds = [sp.simplify(v/(1-v)) for v in alternate]
    residual = sp.simplify(odds[1]/odds[0]-sp.exp(-ENERGIES[0][1]))
    ledger.check('arbitrary smoothing occupations fail physical Gibbs ratio', residual != 0)
    # Even an effective temperature/chemical potential cannot fit these three odds.
    effective_beta_gap = sp.log(odds[0]/odds[1])/sp.log(2)-sp.log(odds[0]/odds[2])/sp.log(5)
    ledger.check('smoothing fixture is not FD at another common beta and mu',
                 abs(float(effective_beta_gap)) > 1e-6)
    return serial({'fixtures': outputs, 'non_FD_occupation_mutant': alternate,
                   'Gibbs_ratio_residual': residual, 'effective_beta_mismatch': effective_beta_gap})
