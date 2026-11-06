"""Independent dense matrix-exponential free energies, without FD response formulas."""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq

from substrate_framework.finite_fermion_state import second_quantize
from fixtures import OPERATORS


def free_energy_functions(energies):
    h = np.array(second_quantize(np.diag(energies)), complex)
    ops = [np.array(second_quantize(a), complex) for a in OPERATORS]
    n = np.diag([s.bit_count() for s in range(8)])
    indices = [s for s in range(8) if s.bit_count() == 1]
    def shifted(x):
        return h-sum((v*a for v, a in zip(x, ops)), np.zeros_like(h))
    def grand(x):
        return -np.log(np.trace(expm(-shifted(x))).real)
    def canonical(x):
        block = shifted(x)[np.ix_(indices, indices)]
        return -np.log(np.trace(expm(-block)).real)
    target = np.trace(n@expm(-h)).real/np.trace(expm(-h)).real
    def legendre(x):
        matrix = shifted(x)
        def residual(mu):
            weights = expm(-matrix+mu*n)
            return np.trace(n@weights).real/np.trace(weights).real-target
        mu = brentq(residual, -10, 10, xtol=1e-13)
        return -np.log(np.trace(expm(-matrix+mu*n)).real)+mu*target
    return grand, legendre, canonical, float(target)


def negative_hessian(function, step):
    origin = np.zeros(2)
    result = np.empty((2, 2))
    center = function(origin)
    for a in range(2):
        da = np.eye(2)[a]*step
        result[a, a] = -(function(da)-2*center+function(-da))/step**2
        for b in range(a):
            db = np.eye(2)[b]*step
            value = -(function(da+db)-function(da-db)-function(-da+db)
                      +function(-da-db))/(4*step**2)
            result[a, b] = result[b, a] = value
    return result
