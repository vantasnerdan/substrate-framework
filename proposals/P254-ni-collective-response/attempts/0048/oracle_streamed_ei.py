"""New all-external output/recovery boundary, using actual unchanged E-I physics."""
import json
from pathlib import Path
import sys
import tempfile
import shutil
from contextlib import contextmanager
from types import SimpleNamespace
import numpy as np
from run_all_external_ei import compute_streamed, make_store, own_banner
from material_electron_ion_inputs import FrozenReference
from run_material_electron_ion_kernel import compute
from reduced_projection_cache import ProjectionCache, atomic_json
from substrate_framework.streamed_material_state import file_binding

@contextmanager
def owning_directory():
    directory = Path(tempfile.mkdtemp(prefix='p254-a48-owning-'))
    try:
        yield directory
    except BaseException:
        print(json.dumps(dict(failed_physical_fixture_retained=str(directory))), flush=True)
        raise
    else:
        shutil.rmtree(directory)



def fixture():
    rng = np.random.default_rng(4801)
    ns, m, d, nq = 4, 2, 3, 2
    h = rng.normal(size=(ns,m,m)) + 1j*rng.normal(size=(ns,m,m))
    h = .11*(h+h.conj().swapaxes(-1,-2))
    rho = np.asarray([[[.25+.02*s,.10+.04j],[.10-.04j,.72-.03*s]] for s in range(ns)])
    target = np.asarray([[0,1,2,3],[1,0,3,2]])
    banks = []
    for q in range(nq):
        g = .03*(rng.normal(size=(ns,d,m,m))+1j*rng.normal(size=(ns,d,m,m)))
        if q == 0:
            g = g+g.conj().swapaxes(-1,-2)
        else:
            for a,b in ((0,1),(2,3)):
                g[b] = g[a].conj().swapaxes(-1,-2)
        banks.append(SimpleNamespace(ionic=g, ionic_reverse=g.conj(),
            density=np.zeros(1), potential=np.zeros(1), kernel=np.zeros(1)))
    j = np.block([[np.zeros((d,d)),np.eye(d)],[-np.eye(d),np.zeros((d,d))]])
    generators, w0 = [], []
    for q in range(nq):
        frequencies = np.asarray([.7,1.1,1.4])+.2*q
        generators.append(j@np.diag(np.r_[frequencies**2,np.ones(d)]))
        w0.append(np.diag(np.r_[1/(2*frequencies),frequencies/2])+.5j*j)
    c = SimpleNamespace(ns=ns,nq=nq,m=m,d=d,n=2,h0=h,rho0=rho,target=target,
        slot_spin=np.asarray([1,1,2,2]),slot_k=np.asarray([0,1,0,1]),
        minus=np.asarray([0,1]),q_indices=np.asarray([[0],[1]]),
        mechanical=SimpleNamespace(j=np.asarray([j,j]),generators=np.asarray(generators),
                                   w0=np.asarray(w0)),q=lambda q:banks[q])
    return c, FrozenReference(c,[0.,1.,10.,30.])


class Interrupted(RuntimeError):
    pass


def main():
    print(json.dumps(own_banner()), flush=True)
    destination = Path(sys.argv[1])
    c, reference = fixture()
    nt = len(reference.times)
    spec = dict(fixture='actual-complex-two-spin-nonzero-Nyquist-q-EI',nq=2,
                sigma_shape=[nt,nt,c.ns,c.m,c.m],pi_shape=[nt,nt,c.d,c.d],
                original_rho='fixed-complex-coherent',
                source=file_binding(__file__))
    with owning_directory() as directory:
        dense, count, reused, _ = compute(c,reference,list(range(c.ns)),[0,1],
                                         ProjectionCache(directory/'old',spec))
        def interrupt(_):
            raise Interrupted('actual physical q committed before simulated interruption')
        store = make_store(directory/'streamed',spec)
        try:
            compute_streamed(c,reference,store,interrupt)
        except Interrupted:
            pass
        _, first_done, _ = make_store(directory/'streamed',spec).restore()
        assert first_done == 1
        try:
            compute_streamed(c,reference,make_store(directory/'streamed',spec),interrupt)
        except Interrupted:
            pass
        before = make_store(directory/'streamed',spec)
        _, second_done, newest = before.restore()
        assert second_done == 2
        # A real corrupt newest numeric archive must use the authenticated previous generation.
        Path(newest['archive']['path']).write_bytes(b'corrupted newest owned temporary archive')
        recovered = make_store(directory/'streamed',spec)
        _, fallback_done, _ = recovered.restore()
        assert fallback_done == 1 and recovered.recoveries
        actual, resumed = compute_streamed(c,reference,recovered,lambda _:None)
        sigma_errors = {key:float(np.linalg.norm(actual['sigma'][key]-dense['sigma_'+key]))
                        for key in actual['sigma']}
        pi_errors = {key:float(np.linalg.norm(np.stack([
            recovered.pi_values(q)[key] for q in range(2)],axis=2)-dense['pi_'+key]))
            for key in actual['sigma']}
        assert not any(sigma_errors.values()) and not any(pi_errors.values())
        completed, completed_counts = compute_streamed(c,reference,
            make_store(directory/'streamed',spec),lambda _:None)
        assert completed_counts['new_edge_pairs'] == 0 and completed_counts['reused_q'] == 2
        assert resumed['new_q'] == 1 and resumed['reused_q'] == 1
        source_mismatch = False
        try:
            make_store(directory/'streamed',dict(spec,fixture='wrong-state'))
        except ValueError:
            source_mismatch = True
        assert source_mismatch
        # Make an internally hash-consistent but physically incomplete result; coverage must refuse it.
        mutated = make_store(directory/'streamed',spec)
        mutated.pi.manifest['items'].pop(mutated.key(1))
        mutated.pi.publish()
        result_path = mutated.directory/'result.json'
        result = json.loads(result_path.read_bytes())
        result['pi_manifest'] = file_binding(mutated.pi.path)
        atomic_json(result_path,result)
        missing_q = False
        try:
            mutated.consume()
        except ValueError:
            missing_q = True
        assert missing_q
        report = dict(schema='p254-a48-owning-output-recovery-oracle-v1',status='passed',
            physical_fixture=dict(explicit_spin_k=4,q=2,bands=2,modes=3,ordered_time_pairs=16,
                                  original_coherence_retained=True),
            unchanged_dense_edge_pairs=count, dense_reused_pairs=reused,
            sigma_errors=sigma_errors,pi_errors=pi_errors,
            actual_partial_done=first_done,actual_second_committed_before_failure=second_done,
            corrupt_newest_previous_generation=fallback_done,resumed_counts=resumed,
            completed_counts=completed_counts,source_mismatch_refused=source_mismatch,
            internally_consistent_missing_q_refused=missing_q,
            source=file_binding(__file__),
            carried='Original A41 Fock/SymPy functional/JVP/transpose proof, not rerun')
    with destination.open('x') as stream:
        json.dump(report,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(result=file_binding(destination),status='passed',
                          sigma_errors=sigma_errors,pi_errors=pi_errors)),flush=True)


if __name__ == '__main__':
    main()
