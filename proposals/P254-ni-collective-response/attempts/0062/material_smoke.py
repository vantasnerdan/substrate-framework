"""Single bounded real-material q1/pair smoke; no full-history physical run."""
from pathlib import Path
import os
import sys
import time
import json
import resource
import subprocess
import traceback
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from history_functional import produce
import numpy as np
from history_banks import OwnedKernels, binding, authenticate, save, ORDERS
from history_sources import open_reference
from joint_material_inputs import load_inputs
from dressed_functionals import common_pair, SOURCE
from material_electron_ion_inputs import edge_pair
from electron_ion_functional import evaluate_components
from cubic_contour import CubicContour, ordered_matrix
from streamed_ee import read_blocks

OUT = Path('/var/tmp/p254-a62-history-smoke-001')
DATASET = Path('/var/tmp/p254-a58-common-feedback-step-002/actual-material-step.npz')
BUNDLE = Path('/var/tmp/p254-a52-requested-reference-001/producer/bundle.json')
PAID = HERE.parent / '0060/pilot-001/dressed3.npyblocks'
UNIT = os.environ['A62_NATIVE_UNIT']
RECEIPTS = Path(os.environ['A62_RECEIPTS_PATH'])


def main():
    start = time.monotonic()
    query = subprocess.run(['systemctl','--user','show',UNIT,
        '--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,InvocationID,ActiveState,SubState'], capture_output=True, text=True)
    props = dict(line.split('=',1) for line in query.stdout.splitlines() if '=' in line)
    cg = Path('/sys/fs/cgroup') / props['ControlGroup'].lstrip('/')
    limits = {k:(cg/k).read_text().strip() for k in ('memory.max','memory.swap.max')}
    try:
        limits['cpuset.cpus.effective'] = (cg/'cpuset.cpus.effective').read_text().strip()
    except FileNotFoundError as error:
        limits['cpuset.cpus.effective'] = dict(unavailable=str(error), affinity_evidence='os.sched_getaffinity')
    banner = dict(argv=sys.argv, pid=os.getpid(), affinity=sorted(os.sched_getaffinity(0)), properties=props,
        kernel=limits, query_exit=query.returncode, query_stdout=query.stdout, query_stderr=query.stderr,
        threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        HERDR={k:v for k,v in os.environ.items() if k.startswith('HERDR_')},
        implementation_bindings=[binding(p) for p in sorted(HERE.glob('*.py'))])
    save(RECEIPTS/'native-start.json',banner)
    print('A62_NATIVE_START ' + json.dumps(banner), flush=True)
    code = 0
    try:
        if os.environ.get('HERDR_ENV') != '1' or banner['affinity'] != [9] or limits['memory.max'] != '3221225472' or limits['memory.swap.max'] != '0':
            raise RuntimeError('actual allocation differs')
        record = produce(OUT/'kernels', dataset=DATASET, owned_q=(1,), owned_pairs=((2,0),),
                         partial=True, pair_batch=1, external_batch=2, static_slots=16)
        paid = json.loads(Path(str(PAID)+'.json').read_bytes())
        if paid['manifest']['source_bindings'][1] != binding(DATASET):
            raise ValueError('paid EE source binding mismatch')
        expected_sigma = {o:np.zeros((128,12,12),complex) for o in ORDERS}
        counts = np.zeros(128,int)
        for block in read_blocks(PAID):
            if block.q != 1:
                raise ValueError('paid q scope differs')
            for j,pair in enumerate(block.pairs):
                if pair == (2,0):
                    counts[list(block.external)] += 1
                    for o in ORDERS:
                        expected_sigma[o][list(block.external)] += block.values[o+'_total'][:,j]
        if not np.all(counts == 1):
            raise ValueError('paid pair does not cover all128 slots exactly once')
        c, baseline, old, bindings = load_inputs()
        old.close()
        ref = open_reference(c,bindings['manifest'],dataset=DATASET)
        try:
            data=c.q(1)
            for source in range(c.ns):
                target=int(c.target[1,source])
                pair=edge_pair(data,source,target,ref.electron(2,0),ref.electron(0,2),ref.boson(1,2,0))
                values=evaluate_components(pair,n=c.n,t=ref.times[2],s=ref.times[0],components=('sigma',))
                for o in ORDERS:
                    expected_sigma[o][target] += values['sigma_'+o]
            full=common_pair(c,ref,2,0)
            ionic=common_pair(c,ref,2,0,True)
            expected_pi={o:full['pi_'+o].copy() for o in ORDERS}
            for o in ORDERS:
                expected_sigma[o] += full['sigma_'+o]-ionic['sigma_'+o]
            cubic=CubicContour.from_source(SOURCE)
            partner,left,right=cubic.pair(0,1)
            first,second=ref.boson(1,2,0),ref.boson(partner,2,0)
            cg=ordered_matrix(left,right,first[0],second[0])
            cl=ordered_matrix(left,right,first[1],second[1])
            for o,value in (('greater',cg),('lesser',cl),('retarded',cg-cl)):
                expected_pi[o][:111,:111] += value
        finally:
            ref.close()
        errors={}
        with OwnedKernels(OUT/'kernels/catalog.json',require_complete=False) as saved:
            for kind,expected in (('sigma',expected_sigma),('pi',expected_pi)):
                bank=getattr(saved,kind)
                for o in ORDERS:
                    error=float(np.max(np.abs(bank[o][2,0]-expected[o])))
                    errors[kind+'_'+o]=error
                    if not np.allclose(bank[o][2,0],expected[o],rtol=3e-12,atol=3e-12):
                        raise AssertionError((kind,o,error))
                    for a in range(3):
                        for b in range(3):
                            if (a,b)!=(2,0) and np.any(bank[o][a,b]):
                                raise AssertionError('unowned pair populated')
            maps=[x for bank in saved.banks.values() for x in bank.values()]
        if not all(x._mmap.closed for x in maps):
            raise AssertionError('reader mappings retained')
        rejected=False
        try:
            with OwnedKernels(OUT/'kernels/catalog.json'):
                pass
        except ValueError as e:
            rejected='partial' in str(e)
        if not rejected:
            raise AssertionError('partial admitted for full solve')
        # Genuine U31 reader exercised at original nodes, no propagation or
        # full31 functional calculation; selected nodes are never fabricated.
        bare=open_reference(c,bindings['manifest'],bundle=BUNDLE,times=[0.,17.,30.])
        try:
            bg,bl=bare.electron(2,1)
            bare_evidence=dict(times=bare.times.tolist(),electron_shape=list(bg.shape),
                greater_norm=float(np.linalg.norm(bg)),lesser_norm=float(np.linalg.norm(bl)),
                u_integration_calls=bare.u_integration_calls,bindings=bare.bindings)
        finally:
            bare.close()
        result=dict(status='passed',scope='actual q1, ordered pair(2,0), genuine times[0,.5,1]; all128 slots/all12 bands/inner64/all114 Gamma; explicitly partial',
            catalog=binding(OUT/'kernels/catalog.json'), paid_ee_catalog=binding(str(PAID)+'.json'),
            paid_numeric=paid['numeric'], paid_index=paid['index'], input_dataset=binding(DATASET),
            comparison='EE from source-bound A60 admitted saved q1 blocks; EI/joint/cubic via existing owning contractions independently assembled for owned pair',
            maxabs_errors=errors,partial_full_solve_rejected=rejected,unowned_pairs_zero=True,
            closed_reader_mappings=True,bare_u31_adapter=bare_evidence,
            producer_timings=record['timings'],producer_elapsed_seconds=record['elapsed_seconds'],
            producer_peak_rss_bytes=record['peak_rss_bytes'],full31_produced=False)
        save(HERE/'actual-material-smoke-001.json',result)
    except BaseException:
        traceback.print_exc()
        code=1
    save(RECEIPTS/'native-exit.json',dict(exitcode=code,elapsed_seconds=time.monotonic()-start,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        cgroup_peak_bytes=(Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')/'memory.peak').read_text().strip()))
    print('A62_NATIVE_EXIT_'+str(code),flush=True)
    return code


if __name__=='__main__':
    sys.exit(main())
