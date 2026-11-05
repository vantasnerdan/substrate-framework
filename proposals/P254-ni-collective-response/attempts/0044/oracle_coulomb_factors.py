"""Independent full physical-grid lrpa and plane-wave normalization oracle."""
import numpy as np
from bare_coulomb_factors import coulomb_multiplier
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.spin_pbe_kernel import apply_spin_pbe_arrays


def run():
    shape = (8,8,8)
    miller = np.indices(shape).reshape(3,-1)
    miller = np.where(miller <= 3, miller, miller-8)
    omega, n, tpiba = 125., 64, 2*np.pi/5
    reciprocal = np.eye(3)*tpiba
    fft = PhysicalRhoFFT.create(shape,miller,reciprocal,omega).at_q(np.array([.17,.09,-.11]))
    weights, admitted = coulomb_multiplier(fft.wavevectors.T,tpiba)
    rng = np.random.default_rng(4401)
    density = rng.normal(size=(fft.nnr,2))+1j*rng.normal(size=(fft.nnr,2))
    actual = apply_spin_pbe_arrays(fft,density,grho=None,coefficients={},tpiba=tpiba,lrpa=True)
    expected = fft.from_g((weights*fft.to_g(density.sum(axis=1)[:,None])[:,0])[:,None])
    assert np.allclose(actual, np.repeat(expected,2,axis=1),rtol=2e-12,atol=2e-12)
    head, mask = coulomb_multiplier(np.array([[0.,0.,0.],[tpiba*1e-5,0,0],[tpiba*2e-4,0,0]]),tpiba)
    assert np.array_equal(mask,[False,False,True]) and np.array_equal(head[:2],[0,0])
    # Plane waves normalized over the full N-cell supercell give 1/(N Omega).
    d = 1/np.sqrt(omega)
    vertex = weights[1]*d*d/n
    direct_integral = weights[1]/(n*omega)
    assert abs(vertex-direct_integral) < 1e-14
    mutations = dict(extra_volume=vertex/omega, missing_N=vertex*n,
                     spin2=2*vertex, repeated_N=vertex/n)
    assert all(abs(v-direct_integral)>abs(direct_integral)*.1 for v in mutations.values())
    assert not mask[0] and (8*np.pi/(tpiba*1e-5)**2)>0
    return dict(schema='p254-a44-factor-oracle-v1',physical_lrpa_relative=float(np.linalg.norm(actual-expected)/np.linalg.norm(actual)),
                normalization_vertex=float(vertex),mutations_exposed=list(mutations)+['exact-zero-only-head'])
