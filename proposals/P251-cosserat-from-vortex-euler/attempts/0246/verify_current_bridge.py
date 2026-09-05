"""Exact current/virtual-work and independent-state identities for0246."""
import sympy as s

from substrate_framework.micropolar import (
    MicropolarCoefficients,
    isotropic_micropolar_energy,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger('P251-0246-periodic-current-action-bridge')
    t = s.Symbol('t', real=True)
    nu = s.Symbol('nu', positive=True)
    acoustic = s.Matrix([[1,t],[0,1]])
    optical = s.Matrix([[s.cos(nu*t),s.sin(nu*t)/nu],
                        [-nu*s.sin(nu*t),s.cos(nu*t)]])
    checks.check('Eq3: both leading physical branch state maps have unit Wronskian',
                 acoustic.det() == 1 and s.trigsimp(optical.det()-1) == 0)
    mu, alpha, cs, ca, ctr = s.symbols('mu alpha cs ca ctr', real=True)
    h = s.Matrix(3,3,lambda i,j:s.Symbol(f'h{i}{j}',real=True))
    g = s.Matrix(3,3,lambda i,j:s.Symbol(f'g{i}{j}',real=True))
    phi = s.Matrix(s.symbols('phi:3',real=True))
    coeff = MicropolarCoefficients(0,mu,alpha,ctr,cs,ca)
    energy = isotropic_micropolar_energy(h,phi,g,coeff)
    force = s.Matrix(3,3,lambda i,j:s.diff(energy,h[i,j]))
    couple = s.Matrix(3,3,lambda i,j:s.diff(energy,g[i,j]))
    ax = lambda m:s.Matrix([sum(s.LeviCivita(i,j,k)*m[j,k]
                               for j in range(3) for k in range(3)) for i in range(3)])
    checks.check('Eq8: differentiated local action has physical angular derivative ax(force)',
                 s.simplify(energy.diff(phi)-ax(force)) == s.zeros(3,1))
    expected_couple = 2*ctr*s.trace(g)*s.eye(3)+cs*(g+g.T)+ca*(g-g.T)
    checks.check('Eq8: couple stress is derived from every retained gradient term',
                 s.simplify(couple-expected_couple) == s.zeros(3))
    v = s.Matrix(s.symbols('v:3',real=True))
    rho = s.Symbol('rho',positive=True)
    transport = rho*v*v.T
    checks.check('Eq4: retaining actual convective transport does not change axial stress',
                 ax(transport) == s.zeros(3,1))
    x = s.symbols('x:3',real=True)
    F = s.Matrix(3,3,lambda i,j:s.Function(f'F{i}{j}')(*x))
    N = s.Matrix(3,3,lambda i,j:s.Function(f'N{i}{j}')(*x))
    virtual = s.Matrix([s.Function(f'v{i}')(*x) for i in range(3)])
    rotation = s.Matrix([s.Function(f'r{i}')(*x) for i in range(3)])
    gradient = lambda a:s.Matrix(3,3,lambda i,j:s.diff(a[i],x[j]))
    div = lambda m:s.Matrix([sum(s.diff(m[i,j],x[j]) for j in range(3)) for i in range(3)])
    contract = lambda a,b:sum(a[i,j]*b[i,j] for i in range(3) for j in range(3))
    work = contract(F,gradient(virtual))+contract(N,gradient(rotation))+ax(F).dot(rotation)
    bulk = -div(F).dot(virtual)+(-div(N)+ax(F)).dot(rotation)
    boundary_vector = F.T*virtual+N.T*rotation
    boundary_div = sum(s.diff(boundary_vector[j],x[j]) for j in range(3))
    checks.check('Eq8: complete periodic virtual-work difference is residual plus actual boundary',
                 s.simplify(work-bulk-boundary_div) == 0)
    checks.check('dropping antisymmetric-force torque changes the bulk virtual-work residual',
                 s.simplify(work-ax(F).dot(rotation)-bulk-boundary_div) != 0)
    checks.check('dropping boundary flux changes the pointwise identity on a cut domain',
                 s.simplify(work-bulk) != 0)
    return checks.finish()


if __name__ == '__main__':
    raise SystemExit(main())
