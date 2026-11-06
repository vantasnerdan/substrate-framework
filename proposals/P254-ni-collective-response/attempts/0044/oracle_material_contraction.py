"""Expose factorized spin-block gather/scatter against a global dense small oracle."""
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import numpy as np
from bare_coulomb_factors import BareCoulombFactors
from second_born_functional import terms_from_vertices, evaluate
from material_second_born import (evaluate_material, jvp_material, transpose_material,
                                  evaluate_material_time_grid)
from reduced_projection_cache import ProjectionCache


def run():
    rng = np.random.default_rng(4404)
    random = lambda shape:rng.normal(size=shape)+1j*rng.normal(size=shape)
    d = random((2,3,2,2))
    d[:,0] = 0
    d[:,2] = d[:,1].conj().swapaxes(-1,-2)
    c = SimpleNamespace(n=1,ns=2,nq=1,m=2,density=(d,),slot_spin=np.array([1,2]),slot_k=np.zeros(2,int),
        k_indices=np.zeros((1,3),int),q_indices=np.zeros((1,3),int),
        target=np.array([[0,1]]),minus=np.array([0]),reciprocal=np.eye(3),tpiba=1.,
        feature_g=(np.array([[0,0,0],[1,0,0],[-1,0,0]]),),
        feature_partner=(np.array([0,2,1,3,5,4]),))
    factors = BareCoulombFactors(c)
    forward = (random((2,2,2)),random((2,2,2)))
    backward = (random((2,2,2)),random((2,2,2)))
    direction = (random((2,2,2)),random((2,2,2)))
    back_direction = (random((2,2,2)),random((2,2,2)))
    def global_blocks(value):
        out = np.zeros((4,4),complex)
        out[:2,:2],out[2:,2:] = value
        return out
    o = np.zeros((3,4,4),complex)
    for spin in range(2):
        o[:,2*spin:2*spin+2,2*spin:2*spin+2] = d[spin].conj()
    x = o*factors.weights[0][:,None,None]
    y = o[[0,2,1]]
    terms = terms_from_vertices(((x,y),)*2,((x,y),)*2)
    legs = dict(greater=global_blocks(forward[0]),lesser=global_blocks(forward[1]),
                reverse_greater=global_blocks(backward[0]),reverse_lesser=global_blocks(backward[1]))
    dense = evaluate(terms,legs)
    defects = []
    for external in range(2):
        actual,count = evaluate_material(factors,external,forward,backward)
        assert count['momentum_spin_quartets'] == 3
        for order in dense:
            for channel in dense[order]:
                expected = dense[order][channel][2*external:2*external+2,2*external:2*external+2]
                defects.append(float(np.linalg.norm(actual[order][channel]-expected)))
                assert np.allclose(actual[order][channel],expected,rtol=2e-11,atol=2e-11)
        tangent = jvp_material(factors,external,forward,backward,direction,back_direction)
        dual = {order:{'total':random((2,2))} for order in dense}
        gf,gb = transpose_material(factors,external,forward,backward,dual)
        lhs = sum(np.sum(dual[order]['total']*tangent[order]['total']) for order in dense)
        rhs = sum(np.sum(x*y) for x,y in zip((*gf,*gb),(*direction,*back_direction),strict=True))
        assert abs(lhs-rhs)<2e-9*max(1,abs(lhs))
    # The real q-chunk consumer and immutable resume path on a small bank.
    electrons = {(0,0):forward,(0,1):forward,(1,0):backward,(1,1):backward}
    with TemporaryDirectory(prefix='p254-a44-q-chunk-oracle-') as directory:
        cache = ProjectionCache(directory+'/outputs',dict(oracle='entrywise-bound-q-chunk-v1'))
        grid,counts = evaluate_material_time_grid(factors,0,electrons,[0.,.4],cache)
        for ti in range(2):
            for si in range(2):
                one,_ = evaluate_material(factors,0,electrons[ti,si],electrons[si,ti],retarded=ti>=si)
                for order,channels in one.items():
                    for channel,value in channels.items():
                        assert np.allclose(grid[order+'_'+channel][ti,si],value)
        resumed = ProjectionCache(directory+'/outputs',dict(oracle='entrywise-bound-q-chunk-v1'))
        repeated,reuse = evaluate_material_time_grid(factors,0,electrons,[0.,.4],resumed)
        assert counts['new_quartets'] == 3 and reuse['new_quartets'] == 0 and reuse['reused_quartets'] == 3
        assert counts['new_q_chunks'] == 1 and reuse['reused_q_chunks'] == 1
        assert all(np.array_equal(grid[key],repeated[key]) for key in grid)
    # All original64q directed maps, including Nyquist; independent integer arithmetic.
    from reduced_projection_algebra import q_feature_geometry
    mesh = np.indices((4,4,4)).reshape(3,-1).T
    counts = set()
    for q in range(64):
        gs,p = q_feature_geometry(q,27)
        qm = int(np.ravel_multi_index(tuple(-mesh[q]%4),(4,4,4)))
        gm,pm = q_feature_geometry(qm,27)
        fq = mesh[q]/4-np.floor(mesh[q]/4+.5)
        fm = mesh[qm]/4-np.floor(mesh[qm]/4+.5)
        assert np.array_equal(pm[p],np.arange(len(p)))
        assert np.allclose(gs+fq,-(gm[p[:len(gs)]]+fm))
        counts.add(len(gs))
    assert counts == {27,36,42,46}
    original = SimpleNamespace(n=64,ns=128,k_indices=mesh,slot_k=np.tile(np.arange(64),2))
    original_factors = object.__new__(BareCoulombFactors)
    original_factors.c = original
    for external in (0,37,64,101):
        actual = set(original_factors.quartets(external))
        independent = set()
        for a in range(128):
            for b in range(128):
                spin_c = external//64+a//64-b//64
                if spin_c not in (0,1):
                    continue
                kc = (mesh[external%64]+mesh[a%64]-mesh[b%64])%4
                independent.add((a,b,64*spin_c+int(np.ravel_multi_index(tuple(kc),(4,4,4)))))
        assert len(actual) == 12288 and actual == independent
        q_counts = np.zeros(64,dtype=int)
        for _,b,_ in actual:
            q_counts[original_factors.transfer(external,b)] += 1
        assert np.array_equal(q_counts,np.full(64,192))
    nyquist_g,nyquist_partner = q_feature_geometry(2,27)
    nyquist_q = np.array([0,0,-.5])
    assert not np.allclose(nyquist_g+nyquist_q,-(nyquist_g+nyquist_q))
    assert not np.array_equal(nyquist_partner,np.arange(len(nyquist_partner)))
    return dict(schema='p254-a44-material-bridge-oracle-v1',dense_block_defect=max(defects),
                closure_feature_counts=sorted(counts),scope='global dense explicit-spin bridge and full64q integer feature closure')
