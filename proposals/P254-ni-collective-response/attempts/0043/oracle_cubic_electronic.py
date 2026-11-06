"""SymPy placed-projector scalar oracle, independent of jet/product code.

Synthetic asymmetric complex dyads expose the position-phase/product layer;
this does not stand in for the source-bound Ni/H material CLI.
"""
from itertools import product, permutations
from types import SimpleNamespace as NS
import json
import numpy as np
import sympy as sp
from cubic_electronic_donor import moment_overlaps, product_third, local_third


def main():
    rng = np.random.default_rng(4323); errors = {}; mutations = {}
    p = np.array([[.2,-.3,.7],[-.5,.1,.4],[.3,.8,-.2],[.6,-.7,.1]])
    beta = rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
    left = rng.normal(size=(4,3))+1j*rng.normal(size=(4,3))
    right = rng.normal(size=(4,3))+1j*rng.normal(size=(4,3))
    d = np.array([[.4,.2+.1j],[-.3+.2j,.7]])
    shape = (1,2,3)
    bl = moment_overlaps(beta,left,p,shape); br = moment_overlaps(beta,right,p,shape)
    encoded = product_third(bl,br,d)[0]
    x = sp.symbols('x:3',real=True)
    # Direct scalar Tr[(L R†) B(x) D B(x)†], but D and L R† need not be
    # Hermitian. Only the known phase/B leg is conjugated, never the dyad.
    scalar = 0
    for a,b,h,k,c in product(range(4),range(4),range(2),range(2),range(3)):
        coefficient = right[a,c].conjugate()*beta[a,h]*d[h,k]*beta[b,k].conjugate()*left[b,c]
        scalar += complex(coefficient)*sp.exp(sp.I*sum(float(p[b,j]-p[a,j])*x[j] for j in range(3)))
    direct = np.array([complex(sp.diff(scalar,x[a],x[b],x[c]).subs(dict.fromkeys(x,0)))
                       for a,b,c in product(range(3),repeat=3)]).reshape(3,3,3)
    errors['all_eight_legs'] = float(np.max(np.abs(encoded-direct)))
    for order in permutations(range(3)):
        errors['permutation_'+''.join(map(str,order))] = float(np.linalg.norm(encoded-encoded.transpose(order)))
    # Explicit unequal sector weights and independent original complex states:
    # no implicit spin degeneracy, and weights are applied once OUTSIDE trace.
    second = product_third(moment_overlaps(beta,left*(.7+.2j),p,shape),br,d)[0]
    errors['explicit_weighting'] = float(np.linalg.norm(.2*encoded+.8*second-(.2+.8*(.7+.2j))*direct))
    mutations['spin2'] = float(np.linalg.norm(2*(encoded+second)-(.2*direct+.8*(.7+.2j)*direct)))
    only_end = -1j*np.einsum('ib,ib->',br[()].conj()[0],np.einsum('ij,jb->ib',d,bl[(0,0,0)][0]))
    only_end += 1j*np.einsum('ib,ib->',br[(0,0,0)].conj()[0],np.einsum('ij,jb->ib',d,bl[()][0]))
    mutations['missing_cross_legs'] = float(abs(only_end-direct[0,0,0]))
    mutations['phase_sign'] = float(np.linalg.norm(-encoded-direct))
    mutations['conjugated_density'] = float(np.linalg.norm(product_third(br,bl,d)[0]-direct))
    local_errors, local_mutations = local_oracle(rng)
    errors.update(local_errors); mutations.update(local_mutations)
    if max(errors.values()) > 3e-9 or min(mutations.values()) < 1e-8:
        raise AssertionError(dict(errors=errors,mutations=mutations))
    print(json.dumps(dict(schema='p254-a43-cubic-electronic-oracle-v1',errors=errors,mutations=mutations),indent=2))


def local_oracle(rng):
    grid = (3,3,3); volume = 7.3; tau = np.array([.13,-.21,.17])
    miller = np.array(list(product([-1,0,1],repeat=3))).T
    reciprocal = np.diag([.7,.9,1.1]); g = (reciprocal@miller).T
    vloc = rng.normal(size=len(g))
    density = rng.normal(size=27)+1j*rng.normal(size=27)
    ionic = NS(inputs=NS(density=NS(miller=miller),
               reference=NS(nat=1,atom_names=['X'],positions_bohr=tau[:,None])),g=g,vloc={'X':vloc})
    coupling = NS(ionic=ionic,factors=NS(fft=NS(nnr=27,physical=grid,omega=volume)))
    encoded = local_third(density,coupling)[0]
    # Literal rectangular-grid integration, no FFT in the independent route.
    sites = np.indices(grid).reshape(3,-1,order='F').T/np.array(grid)
    sampled = np.array([np.sum(density*np.exp(2j*np.pi*sites@m))/27 for m in miller.T])
    direct = np.zeros((3,3,3),complex)
    for a,b,c in product(range(3),repeat=3):
        direct[a,b,c] = volume*sum(sampled[j]*vloc[j]*np.exp(-1j*g[j]@tau)
                                  *(-1j*g[j,a])*(-1j*g[j,b])*(-1j*g[j,c]) for j in range(len(g)))
    return ({'local_grid_phase_units':float(np.max(np.abs(encoded-direct)))},
            {'local_volume':float(np.linalg.norm(direct/volume-direct)),
             'local_phase':float(np.linalg.norm(-encoded-direct))})


if __name__ == '__main__':
    main()
