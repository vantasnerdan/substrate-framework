"""Independent SymPy scalar third and literal finite-N coordinate oracle.

The reference uses the actual BvK supercell scalar with all cell coordinates,
not a copied reciprocal-q formula. Its deliberately tiny lattice tests general
three-leg momenta, atom order, units and the canonical N^-1/2 normalization.
"""
from itertools import product, permutations
from math import pi
import json
import numpy as np
import sympy as sp
from cubic_ionic_donor import IonicThird, radial_third


def symbolic_radial():
    x = sp.symbols('x:3',real=True)
    r = sp.sqrt(sum(v*v for v in x))
    scalar = sp.erfc(r)/r
    expressions = [sp.diff(scalar,x[a],x[b],x[c]) for a,b,c in product(range(3),repeat=3)]
    return sp.lambdify([x],expressions,'scipy')


def literal_finite_n(cell,positions,charges,mesh,q,directions,radial,reciprocal_scale=1.):
    cells = np.array(list(product(*(range(n) for n in mesh))))
    n = len(cells); supercell = cell@np.diag(mesh)
    sites = np.array([cell@r+tau for r in cells for tau in positions])
    z = np.tile(charges,n)
    u = [np.array([directions[leg][atom]*np.exp(2j*pi*np.dot(q[leg]/mesh,r))/np.sqrt(n)
                   for r in cells for atom in range(len(charges))]) for leg in range(3)]
    # Uncentered literal image enumeration is independent of the production
    # centered primitive pair list. Generous bounds cover displaced atom pairs.
    inv = np.linalg.inv(supercell)
    bound = np.ceil(12*np.linalg.norm(inv,axis=1)+2).astype(int)
    images = np.array(list(product(*(range(-b,b+1) for b in bound))))@supercell.T
    real = 0j
    for i in range(len(sites)):
        for j in range(len(sites)):
            for image in images:
                d = sites[i]-sites[j]-image
                radius = np.linalg.norm(d)
                if 1e-10 < radius < 12.:
                    t = np.asarray(radial(d)).reshape(3,3,3)
                    real += z[i]*z[j]*np.einsum('abc,a,b,c->',t,*(v[i]-v[j] for v in u))
    reciprocal = 2*pi*np.linalg.inv(supercell).T
    bound = np.ceil(12*np.linalg.norm(supercell,axis=0)/(2*pi)).astype(int)
    ks = np.array(list(product(*(range(-b,b+1) for b in bound))))@reciprocal.T
    total = 0j
    for k in ks:
        k2 = k@k
        if not 0 < k2 <= 144.:
            continue
        weight = np.exp(-k2/4)/k2
        for i in range(len(sites)):
            for j in range(len(sites)):
                phase = np.exp(1j*k@(sites[i]-sites[j]))
                jet = np.prod([1j*k@(v[i]-v[j]) for v in u])
                total += z[i]*z[j]*weight*phase*jet
    return real, reciprocal_scale*4*pi/abs(np.linalg.det(supercell))*total


def main():
    radial = symbolic_radial(); errors = {}; mutations = {}
    d = np.array([.7,-1.1,.4])
    errors['scalar_radial'] = float(np.max(np.abs(radial_third(d)-np.asarray(radial(d)).reshape(3,3,3))))
    cell = np.diag([4.1,4.7,5.3])
    positions = np.array([[.1,.2,.3],[1.3,1.7,.8],[2.4,.6,2.1]])
    charges = np.array([2.,1.,3.]); mesh = np.array([3,1,1]); q = np.array([[1,0,0],[1,0,0],[1,0,0]])
    rng = np.random.default_rng(4312)
    directions = rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
    donor = IonicThird(cell,positions,charges,mesh)
    direct = literal_finite_n(cell,positions,charges,mesh,q,directions,radial)
    encoded = []
    for piece,reference in zip(('real','reciprocal'),direct):
        block = donor.block(q,piece=piece)
        value = np.einsum('abc,a,b,c->',block,*(v.ravel() for v in directions))
        encoded.append(value)
        errors[piece+'_finiteN'] = float(abs(value-reference))
    total = sum(encoded)
    for order in permutations(range(3)):
        block = donor.block(q[list(order)])
        value = np.einsum('abc,a,b,c->',block,*(directions[i].ravel() for i in order))
        errors['permutation_'+''.join(map(str,order))] = float(abs(value-total))
    # Distinct, genuinely non-Gamma q legs; opposite/Nyquist labels remain
    # exact integers, including their reciprocal equivalence on the mesh.
    q2 = np.array([[1,0,0],[0,0,0],[2,0,0]])
    reference = sum(literal_finite_n(cell,positions,charges,mesh,q2,directions,radial))
    value = np.einsum('abc,a,b,c->',donor.block(q2),*(v.ravel() for v in directions))
    errors['general_q'] = float(abs(value-reference))
    mutations['phase_sign'] = float(abs(np.einsum('abc,a,b,c->',donor.block(-q2),
                                                 *(v.ravel() for v in directions))-reference))
    mutations['missing_N'] = float(abs((np.sqrt(3)-1)*value))
    mutations['Ry_Hartree'] = float(abs(total/2-total))
    mutations['volume'] = float(abs(direct[1]*(abs(np.linalg.det(cell))-1)))
    wrong = IonicThird(cell,positions,charges[::-1],mesh)
    mutations['charge_order'] = float(abs(np.einsum('abc,a,b,c->',wrong.block(q2),
                                                   *(v.ravel() for v in directions))-reference))
    mutations['sign'] = float(abs(2*value))
    # Atomic folding is integer modulo, never nearest float/opposite conjugacy.
    errors['qfold'] = float(np.linalg.norm(donor.block(q2)-donor.block(q2+mesh*np.array([[1],[0],[-1]]))))
    if max(errors.values()) > 2e-8 or min(mutations.values()) < 1e-9:
        raise AssertionError(dict(errors=errors,mutations=mutations))
    print(json.dumps(dict(schema='p254-a43-cubic-ionic-oracle-v1',errors=errors,mutations=mutations),indent=2))


if __name__ == '__main__':
    main()
