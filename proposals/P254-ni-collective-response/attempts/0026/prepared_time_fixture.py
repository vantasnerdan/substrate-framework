"""Tiny original-basis/Fock comparator. Never imports a production RHS.

One finite ion has all three canonical coordinates; two spin sectors have
three original orbitals each, two occupied-factor columns and a full complement.
"""
from types import SimpleNamespace
import numpy as np


def fixture():
    from substrate_framework.physical_factor_preparation import ActionFactors, FactorFunctional, FactorPreparation
    from substrate_framework.prepared_reference import PreparedReference, PreparedSector
    from substrate_framework.prepared_factor_dynamics import PreparedFactorDynamics
    rng = np.random.default_rng(2252603)
    def herm(scale):
        a = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        return scale*(a+a.conj().T)
    h = np.array([herm(.13), herm(.17)])
    g = np.array([[herm(.09) for _ in range(3)] for _ in range(2)])
    v = rng.normal(size=(2, 2, 3)) + 1j*rng.normal(size=(2, 2, 3))
    density = np.einsum('sri,srj->srij', v, v.conj())/11
    kernel = np.array([[.7,.13,-.11,.08],[.13,.9,.07,-.12],[-.11,.07,.6,.19],[.08,-.12,.19,.8]])
    transforms = (np.array([[1.2,.17+.09j],[.04j,.8]]), np.array([[.9,.12j],[.19,1.1]]))
    occupations = np.array([[.35,.35,0.],[.6,.2,0.]])
    f = SimpleNamespace(h=h, g=g, density=density, kernel=kernel, transforms=transforms,
                        occupations=occupations, weight=.25, ncells=4, dv=1.7, ndof=3)

    class Coupling:
        npw, nr, ndof = 3, 2, 3
        grid_measure_bohr3, cell_volume_bohr3 = 1.7, 3.4
        def __init__(self, spin):
            self.spin = spin
        def moments(self, left, right, coefficient=None):
            rho = left @ right.conj().T if coefficient is None else left @ coefficient @ right.conj().T
            return np.einsum('rij,ji->r', density[self.spin], rho), np.einsum('aij,ji->a', g[self.spin], rho)
        def potential_action(self, field, vectors):
            return f.dv*np.einsum('r,rij->ij', field, density[self.spin]) @ vectors
        def ionic_action(self, direction, vectors):
            return np.einsum('a,aij->ij', direction, g[self.spin]) @ vectors

    sectors = []
    for spin, transform in enumerate(transforms):
        psi = np.eye(3, dtype=complex)[:, :2] @ transform
        inv = np.linalg.inv(transform)
        d = inv @ np.diag(occupations[spin, :2]) @ inv.conj().T
        label = 'time-oracle-original-spin-' + str(spin)
        action = ActionFactors(psi, psi, 'I', 'I', '1', label, 'finite original identity')
        functional = FactorFunctional.from_actions(psi, {'I': action}, factor_identity=label+'-raw', pw_basis_identity=label)
        prep = FactorPreparation(functional, d, 'independent Bernoulli occupations and empty original complement')
        sectors.append(PreparedSector(prep, psi, Coupling(spin), lambda x, a=h[spin]: a@x,
                                      spin, f.weight, label, 'explicit finite Hermitian physical matrices'))
    f.reference = PreparedReference(
        nbar=np.array([[.31,.23],[.18,.36]]), ionic_force=np.array([.19,-.13,.05]),
        held_density=np.array([[.04,.02],[.03,.06]]),
        held_gradient=np.array([-.03,.05,.01]), curvature=np.array([[1.3,.21,.07],[.21,.85,-.09],[.07,-.09,1.1]]),
        masses=np.array([2.7,5.1,3.6]), support_force=np.array([.02,.01,-.01]),
        kernel=lambda n:(kernel@n.ravel()).reshape(2,2), ncells=4,
        expansion_authority='explicit finite joint-second-order scalar',
        ionic_authority='explicit synthetic linear ionic scalar; not Ni/H or equilibrium',
        coverage_authority='two mixed-feedback spins, complete original complement and held remainder',
        support_authority='explicit constant conservative synthetic support',
        coordinate_identity='one synthetic ion, all three canonical Ry coordinates',
        state_authority='independent quasifree Bernoulli spin marginals, not fixed particle number')
    f.model = PreparedFactorDynamics(f.reference, sectors)
    f.mean = f.model.initial_mean(np.array([.11,-.07,.04]), np.array([.08,-.13,.06]))
    f.b = np.array([herm(.21)+1j*herm(.09), herm(.19)+1j*herm(.08)])
    f.bi = np.array([-.1j,.25,.04+.03j,.17,-.08+.02j,.06j])
    f.displacement = np.array([.18,-.09,.04,.12,.07,-.03])
    widths = np.array([.7,1.1,.9])
    f.sigma = np.block([[np.diag(widths), .5j*np.eye(3)],[-.5j*np.eye(3), np.diag(1/(4*widths))]])
    return f


def fock_data(f):
    car = []
    for mode in range(3):
        a = np.zeros((8,8), complex)
        for bits in range(8):
            if bits & (1<<mode):
                a[bits^(1<<mode), bits] = (-1)**((bits & ((1<<mode)-1)).bit_count())
        car.append(a)
    def lift(a):
        return sum(a[i,j]*car[i].conj().T@car[j] for i in range(3) for j in range(3))
    weights = [np.array([np.prod([occ[j] if bits&(1<<j) else 1-occ[j] for j in range(3)])
                         for bits in range(8)]) for occ in f.occupations]
    w = np.kron(*weights)
    eb = np.sqrt(f.ncells)*f.weight*(np.kron(lift(f.b[0]),np.eye(8))+np.kron(np.eye(8),lift(f.b[1])))
    means = (np.array([(bits//8)&1 for bits in range(64)])-f.occupations[0,0])[:,None]*f.displacement
    return car, w, eb, means


def initial_sources(f):
    """Production source versus independently literal CAR/canonical/joint moments."""
    from substrate_framework.prepared_covariance_source import PreparedProbe, PreparedCovariance, FactorCrossPreparation
    from substrate_framework.physical_bilinear_preparation import IonicPreparation
    b = PreparedProbe(tuple(lambda x,a=a:a@x for a in f.b),
                      tuple(lambda x,a=a:a.conj().T@x for a in f.b), f.bi[:3], f.bi[3:], -.07j, 'complex B', 'synthetic')
    car, weights, eb, means = fock_data(f)
    def expect(a):
        return weights@np.diag(a)
    classical = means.T@(weights[:,None]*means)
    zero = np.zeros((2,3,3), complex)
    sources, direct = [], []
    for joint in (False, True):
        sigma = f.sigma + (classical if joint else 0)
        ions = IonicPreparation(sigma, 'full-line Gaussian; conditional displacement if joint')
        cross = FactorCrossPreparation.product(f.model, 'explicit electron/ion product')
        if joint:
            coeff = [np.zeros((6,2,2),complex) for _ in range(2)]
            inv = np.linalg.inv(f.transforms[0])
            raw = inv@np.diag([f.occupations[0,0]*(1-f.occupations[0,0]),0])@inv.conj().T
            coeff[0] = np.sqrt(f.ncells)*f.displacement[:,None,None]*raw
            cross = FactorCrossPreparation(tuple(coeff), 'literal occupation-conditioned displaced Gaussian mixture',
                                           'independent quasifree spin marginals preserved; cross in raw span')
        state_label = 'conditioned joint Gaussian mixture' if joint else 'explicit product Gaussian'
        prep = PreparedCovariance(f.model, electronic_authority='literal independent Bernoulli Fock weights; '+state_label,
                                  ionic=ions, cross=cross)
        for ordering in ('AB','BA'):
            sources.append(prep.source(b, ordering=ordering))
            electronic, crossrho = zero.copy(), zero.copy()
            for spin in range(2):
                for i in range(3):
                    for j in range(3):
                        op = car[j].conj().T@car[i]
                        x = np.sqrt(f.ncells)*(np.kron(op,np.eye(8)) if spin == 0 else np.kron(np.eye(8),op))
                        diagonal_product = np.einsum('ij,ji->i', x, eb) if ordering == 'AB' else np.einsum('ij,ji->i', eb, x)
                        electronic[spin,i,j] = weights@diagonal_product-expect(x)*expect(eb)
                        if joint:
                            crossrho[spin,i,j] = weights@((np.diag(x)-expect(x))*(means@f.bi))
            y = (sigma if ordering == 'AB' else sigma.T)@f.bi
            cy = means.T@(weights*(np.diag(eb)-expect(eb))) if joint else np.zeros(6,complex)
            direct.extend((pack(electronic,np.zeros(3),np.zeros(3)), pack(zero,y[:3],y[3:]), pack(crossrho,cy[:3],cy[3:])))
    return tuple(sources), np.array(direct), b


def pack(rho, u, p):
    return np.concatenate((np.asarray(rho).ravel(), u, p)).astype(complex)


def unpack(vector):
    return vector[:18].reshape(2,3,3), vector[18:21], vector[21:24]


def factor_vector(model, mean, tangent=None):
    rho = ([psi@s.d@psi.conj().T for s,psi in zip(model.sectors,mean.frames)] if tangent is None else
           [left@psi.conj().T+psi@right for psi,left,right in zip(mean.frames,tangent.ket,tangent.bra)])
    x = mean if tangent is None else tangent
    return pack(rho,x.u,x.p)


def fields(f, rho, u):
    density = f.weight*np.einsum('srij,sji->rs',f.density,rho)
    gradient = f.weight*np.einsum('saij,sji->a',f.g,rho)
    delta = density+f.reference.held_density-f.reference.nbar
    potential = (f.kernel@delta.ravel()).reshape(2,2)
    h = f.h+f.dv*np.einsum('rs,srij->sij',potential,f.density)+np.einsum('a,saij->sij',u,f.g)
    return h, gradient, delta


def original_rhs(f, vector, *, frozen_connection=False):
    """Complete original density matrices and full linear variations, no factors."""
    rows = vector.reshape(-1,24)
    rho,u,p = unpack(rows[0])
    h, gradient, _ = fields(f,rho,u)
    r = f.reference
    force = r.ionic_force+r.support_force-r.held_gradient-gradient-r.curvature@u
    result = [pack(-1j*(h@rho-rho@h),p/r.masses,force)]
    for row in rows[1:]:
        dr,du,dp = unpack(row)
        dn = f.weight*np.einsum('srij,sji->rs',f.density,dr)
        dg = f.weight*np.einsum('saij,sji->a',f.g,dr)
        potential = (f.kernel@dn.ravel()).reshape(2,2)
        dh = f.dv*np.einsum('rs,srij->sij',potential,f.density)+np.einsum('a,saij->sij',du,f.g)
        transport = 0 if frozen_connection else h@dr-dr@h
        result.append(pack(-1j*(transport+dh@rho-rho@dh),dp/r.masses,-r.curvature@du-dg))
    return np.asarray(result).ravel()


def original_energy(f, row):
    rho,u,p = unpack(row)
    _,gradient,delta = fields(f,rho,u)
    r = f.reference
    return (f.weight*np.einsum('sij,sji->',f.h,rho)+f.dv*delta.ravel()@f.kernel@delta.ravel()/2
            +u@(gradient+r.held_gradient-r.ionic_force-r.support_force)
            +u@r.curvature@u/2+np.sum(p*p/r.masses)/2)


def original_probe(f, row, name, *, connected=False):
    rho,u,p = unpack(row)
    if name == 'force':
        direction = np.array([1.,-.3,.2])
        value = -direction@(f.weight*np.einsum('saij,sji->a',f.g,rho)+f.reference.curvature@u)
        if not connected:
            value += direction@(f.reference.ionic_force+f.reference.support_force-f.reference.held_gradient)
        return value
    region = np.array([1.,.3])
    density = f.weight*np.einsum('srij,sji->rs',f.density,rho)
    if not connected:
        density = density+f.reference.held_density
    return f.dv*np.sum(region[:,None]*density)
