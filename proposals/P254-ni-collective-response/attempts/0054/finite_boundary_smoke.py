"""Changed algorithm only: exact full quadratic comparison and paid actual input."""
from pathlib import Path
import json
import sys
import numpy as np
from scipy.integrate import quad_vec
from supplied_keldysh import solve_supplied, load_solution, A49
from finite_gaussian_reference import fermion_case, boson_case
from command_receipts import binding, run_recorded
from consume_solution import consume

PAID = Path('/var/tmp/p254-a53-fixed-kernel-material-001/fermions/slot-000/input.npz')
PAID_SHA = '776184aa7bdce42247ef1df86c2ceba9e025f7daf149e8eb59d6583f1df175af'


def error(a,b):
    return float(np.max(np.abs(a-b)))


def solve_case(factory, stat, nt, order, binary, root, label, omission=0):
    times = np.linspace(0,.6,nt+1)
    case = factory(times)
    path = solve_supplied(times,*case['source'],*case['kernel'],statistics=stat,order=order,
                          binary=binary,workdir=root/label,omission=omission)
    with np.load(path,allow_pickle=False) as saved:
        result = {name:saved[name] for name in ('ret','less','f_ret','fcc_ret','f_less','fcc_less')}
    report = {'solution':binding(path), 'ret_error':error(result['ret'],case['exact'][0]),
              'less_error':error(result['less'],case['exact'][1]),
              'source_R_diagonal_copy_max':max(error(result['ret'][n,n],case['source'][0][n,n]) for n in range(nt+1)),
              'raw_diagnostics':json.loads((root/label/'output.bin.diagnostics.json').read_text())}
    try:
        report['consumer'] = consume(path)
        report['strict_reader'] = 'accepted'
    except ValueError as exc:
        report['strict_reader'] = 'rejected: '+str(exc)
    return result,case,report


def main(binary,root):
    root = Path(root)
    reports = {}
    report_path = root/'changed-boundary-report.json'
    def persist():
        report_path.write_text(json.dumps(reports,indent=2)+'\n')
    for label,factory,stat in [('fermion',fermion_case,-1),('boson_position',boson_case,1)]:
        for nt in (10,20,40):
            key = f'{label}-{nt}'
            result,case,report = solve_case(factory,stat,nt,3,binary,root,key)
            reports[key] = report
            persist()
        times = np.linspace(0,.6,41)
        def integrand(t):
            probe = factory(np.array([.3,t,.6]))
            q,k = probe['source'][0],probe['kernel'][0]
            return np.stack([-q[2,1]@k[1,0],-k[2,1]@q[1,0]])
        exact,estimate = quad_vec(integrand,.3,.6,epsabs=1e-12,epsrel=1e-12)
        report.update(F_quadrature_error=error(result['f_ret'][-1,20],exact[0]),
                      Fcc_quadrature_error=error(result['fcc_ret'][-1,20],exact[1]),
                      swapped_Fcc_error=error(result['f_ret'][-1,20],exact[1]),
                      quadrature_estimate=float(estimate),
                      preparation_discard_error=error(result['less'],case['discarded_preparation'][1]))
        def lesser_first(t):
            probe = factory(np.array([t,.3,.6]))
            qr,ql = probe['source']
            kr,kl = probe['kernel']
            return np.stack([-qr[1,0]@kl[0,2],-kr[1,0]@ql[0,2]])
        def lesser_second(t):
            probe = factory(np.array([t,.3,.6]))
            qr,ql = probe['source']
            kr,kl = probe['kernel']
            return np.stack([-ql[1,0]@kr[2,0].conj().T,-kl[1,0]@qr[2,0].conj().T])
        first,_ = quad_vec(lesser_first,0.,.3,epsabs=1e-12,epsrel=1e-12)
        second,_ = quad_vec(lesser_second,0.,.6,epsabs=1e-12,epsrel=1e-12)
        lesser_exact = first+second
        report.update(F_lesser_quadrature_error=error(result['f_less'][20,-1],lesser_exact[0]),
                      Fcc_lesser_quadrature_error=error(result['fcc_less'][20,-1],lesser_exact[1]),
                      saved_F_Fcc_companion_exchange_max=error(result['f_less'],-result['fcc_less'].swapaxes(0,1).swapaxes(2,3).conj()))
        for omission,name in [(1,'omit-initial'),(2,'omit-noise'),(3,'wrong-companion')]:
            mutated,_,mutation = solve_case(factory,stat,20,3,binary,root,f'{label}-{name}',omission)
            with np.load(root/f'{label}-20'/'solution.npz',allow_pickle=False) as base:
                mutation['same_grid_lesser_difference'] = error(mutated['less'],base['less'])
            reports[f'{label}-{name}'] = mutation
            persist()
        # Strong coarse-step sensitivity is distinct from the fine invariant test.
        _,_,reports[f'{label}-coarse-order1'] = solve_case(factory,stat,1,1,binary,root,f'{label}-coarse-order1')
        persist()
        run_recorded([sys.executable,str(Path(__file__).with_name('consume_solution.py')),
                      str(root/f'{label}-40'/'solution.npz'),str(root/f'{label}-consumer.json')],
                     root,f'{label}-saved-consumer',60)
    paid_binding = binding(PAID)
    if paid_binding['sha256'] != PAID_SHA or paid_binding['bytes'] != 38914:
        raise ValueError('paid actual input binding changed')
    with np.load(PAID,allow_pickle=False) as original:
        data = {key:original[key] for key in original.files}
    actual = solve_supplied(data['times'],*[data[k] for k in ('source_ret','source_less','kernel_ret','kernel_less')],
                            statistics=int(data['statistics']),order=int(data['order']),beta=float(data['beta']),
                            binary=binary,workdir=root/'actual-slot000')
    # Save raw numbers even if the strict reader rejects; never repair the solution.
    reports['actual-slot000'] = {'original_input':paid_binding,'solution':binding(actual),
        'raw_diagnostics':json.loads((root/'actual-slot000'/'output.bin.diagnostics.json').read_text())}
    persist()
    reports['actual-slot000']['consumer'] = consume(actual)
    reports['actual-slot000']['strict_reader'] = 'accepted'
    persist()
    run_recorded([sys.executable,str(Path(__file__).with_name('consume_solution.py')),str(actual),
                  str(root/'actual-slot000-consumer.json')],root,'actual-saved-consumer',60)
    for label in ('fermion','boson_position'):
        for nt in (10,20):
            coarse,fine = reports[f'{label}-{nt}'],reports[f'{label}-{nt*2}']
            coarse['R_refinement_ratio'] = coarse['ret_error']/max(fine['ret_error'],1e-30)
            coarse['less_refinement_ratio'] = coarse['less_error']/max(fine['less_error'],1e-30)
        persist()
        fine = reports[f'{label}-40']
        assert fine['strict_reader']=='accepted' and fine['ret_error']<2e-6 and fine['less_error']<2e-6, fine
        assert fine['source_R_diagonal_copy_max']==0
        assert fine['F_quadrature_error']<2e-6 and fine['Fcc_quadrature_error']<2e-6
        assert reports[f'{label}-20']['less_refinement_ratio']>2
        assert reports[f'{label}-omit-initial']['less_error']>.05
        assert reports[f'{label}-omit-noise']['less_error']>1e-5
        assert fine['F_lesser_quadrature_error']<2e-6 and fine['Fcc_lesser_quadrature_error']<2e-6
        assert fine['saved_F_Fcc_companion_exchange_max']<1e-12
        assert reports[f'{label}-wrong-companion']['same_grid_lesser_difference']>1e-7
    print('A54_CHANGED_BOUNDARY_AND_PAID_INPUT_COMPLETE',flush=True)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])
