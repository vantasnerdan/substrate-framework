"""Focused new-boundary exercise and actual material consumer, no physics replay."""
from pathlib import Path
import json
import sys
import numpy as np
from production_transport import (FinalBank, ResultView, solve_supplied, load_solution,
    validate, HEADER, binding, atomic_json, run_recorded)

HERE = Path(__file__).resolve().parent
OLD = Path('/var/tmp/p254-a58-common-feedback-step-002')


def free_control(count, dimension=3):
    times = np.arange(count,dtype='<f8')*.125
    h = np.diag(np.arange(1,dimension+1)*.07).astype(complex)
    for a in range(dimension):
        for b in range(a):
            h[a,b] = .01*(a+b+1)+.013j*(a-b)
            h[b,a] = h[a,b].conjugate()
    w,v = np.linalg.eigh(h)
    u = np.asarray([(v*np.exp(-1j*w*t))@v.conj().T for t in times])
    rho = np.diag(np.linspace(.2,.7,dimension)).astype(complex)
    if dimension > 1:
        rho[0,1],rho[1,0] = .021+.017j,.021-.017j
    ret = np.zeros((count,count,dimension,dimension),dtype='<c16')
    less = np.empty_like(ret)
    for n in range(count):
        for m in range(count):
            if n <= m:
                less[n,m] = 1j*u[n]@rho@u[m].conj().T
            else:
                less[n,m] = -less[m,n].conj().T
            if n >= m:
                # Noncanonical diagonal tests rawRdiag transport, no replacement.
                ret[n,m] = (-1j-.03)*(u[n]@u[m].conj().T)
    zero = np.zeros_like(ret)
    return times,[ret,less,zero,zero]


def expect_rejection(label, fn, findings):
    try:
        fn()
    except (ValueError,RuntimeError,FileNotFoundError) as exc:
        findings.append(dict(label=label,rejected=True,error=str(exc)))
    else:
        raise AssertionError('mutation accepted: '+label)


def close_saved(saved):
    for key in ('ret','less'):
        saved[key]._mmap.close()


def boundary(binary, root):
    findings = []
    times, arrays = free_control(5)
    # Nonnative endian, Fortran/strided input test actual byte/axis conversion.
    arrays = [np.asarray(a,dtype='>c16',order='F') for a in arrays]
    base = root/'boundary-bank'
    bank = FinalBank(base,times,3,2,statistics=-1,order=2)
    solve_supplied(times,*arrays,statistics=-1,order=2,binary=binary,bank=bank,slot=0,source_bindings=[binding(__file__)])
    assert np.array_equal(bank.ret[0],arrays[0]) and np.array_equal(bank.less[0],arrays[1])
    assert not list((base/'.staging').glob('*/input.bin')) and not list((base/'.staging').glob('*/result.bin'))
    findings.append(dict(label='complex-F-order-big-endian-rawRdiag',max_ret=0.,max_less=0.,packed_released_after_commit=True))
    expect_rejection('partial-bank-not-loadable',lambda:load_solution(base),findings)
    bank.close()
    with FinalBank(base,resume=True) as bank:
        assert set(bank.completed)=={0}
        expect_rejection('committed-slot-no-overwrite',lambda:solve_supplied(times,*arrays,statistics=-1,order=2,binary=binary,bank=bank,slot=0,source_bindings=[binding(__file__)]),findings)
        solve_supplied(times,*arrays,statistics=-1,order=2,binary=binary,bank=bank,slot=1,source_bindings=[binding(__file__)])
        bank.finalize()
    saved = load_solution(base)
    assert isinstance(saved['ret'],np.memmap) and not saved['ret'].flags.writeable
    assert np.array_equal(saved['ret'][1],arrays[0]) and np.array_equal(saved['less'][1],arrays[1])
    findings.append(dict(label='committed-slot-resume-and-readonly-final-mapping',passed=True,shape=list(saved['ret'].shape)))
    close_saved(saved)
    for label,ts,aa in (
        ('complex64-rejected',times,[a.astype('c8') for a in arrays]),
        ('float32-time-rejected',times.astype('f4'),arrays),
        ('axis-shape-rejected',times,[a[:,:,:2] for a in arrays]),
        ('nonuniform-time-rejected',np.array([0.,.125,.25,.4,.5]),arrays)):
        expect_rejection(label,lambda ts=ts,aa=aa:validate(ts,aa,-1,2,1.),findings)
    aa = [a.copy() for a in arrays]; aa[1][4,1,0,1] += .1j
    expect_rejection('independent-lesser-leg-rejected',lambda:validate(times,aa,-1,2,1.),findings)
    aa = [a.copy() for a in arrays]; aa[0][0,1,1,0] = .01j
    expect_rejection('acausal-retarded-rejected',lambda:validate(times,aa,-1,2,1.),findings)
    # Persist finite asymmetric layout payload; map triangles directly, mutate only new bytes.
    good = root/'boundary-minimal.bin'
    with good.open('xb') as stream:
        stream.write(HEADER.pack(b'P254VIE1',4,3,2,-1,1.));stream.write(times.tobytes())
        for index,a in enumerate(arrays[:2]):
            for n in range(5):
                stream.write(np.asarray(a[:n+1,n] if index%2 else a[n,:n+1],dtype='<c16',order='C').tobytes())
    with ResultView(good) as view:
        assert np.shares_memory(view.ret,view._mapped) and np.shares_memory(view.less,view._mapped)
        r,l = np.empty_like(arrays[0],dtype='<c16'),np.empty_like(arrays[1],dtype='<c16')
        view.copy_into(r,l)
        assert np.array_equal(r,arrays[0]) and np.array_equal(l,arrays[1])
        expect_rejection('destination-dtype-rejected',lambda:view.copy_into(r.astype('c8'),l),findings)
    assert view._mapped is None
    findings.append(dict(label='packed-zero-copy-lifetime',passed=True))
    raw = good.read_bytes()
    for label,data in (('truncated',raw[:-1]),('trailing',raw+b'x'),('badmagic',b'BADMAGIC'+raw[8:])):
        p = root/('boundary-'+label+'.bin');p.write_bytes(data)
        expect_rejection(label,lambda p=p:ResultView(p),findings)
    # This is an input extent guard test, not a physical solve failure.
    try:
        run_recorded([str(binary),str(good),str(root/'invalid-native-output.bin')],root,'expected-native-extent-rejection',10)
    except RuntimeError:
        receipt=json.loads((root/'expected-native-extent-rejection.json').read_text())
        assert receipt['exit']==2
        findings.append(dict(label='native-input-extent-rejected',exit=receipt['exit']))
    else:
        raise AssertionError('native accepted two-component input')
    atomic_json(root/'boundary-findings.json',dict(findings=findings,scope='new IO/layout/lifetime only; free finite fixture, not material evidence'))


def actual(binary,root):
    findings=[]
    for label,old,d,stat in (('fermion-slot000',OLD/'fermions'/'slot-000',12,-1),('full114Gamma3node',OLD/'intrinsic',114,1)):
        original = json.loads((old/'manifest.json').read_text())
        for key in ('input','transport','solution'):
            assert binding(original[key]['path']) == original[key]
        with np.load(old/'input.npz',allow_pickle=False) as data:
            times=data['times']; arrays=[data[k] for k in ('source_ret','source_less','kernel_ret','kernel_less')]
            order,beta=int(data['order']),float(data['beta'])
        target=root/label
        with FinalBank(target,times,d,1,statistics=stat,order=order,beta=beta) as bank:
            receipt=solve_supplied(times,*arrays,statistics=stat,order=order,beta=beta,binary=binary,bank=bank,slot=0,source_bindings=[binding(old/'manifest.json'),original['input'],original['transport'],original['solution']],timeout=180)
            bank.finalize()
        saved=load_solution(target)
        with np.load(old/'solution.npz',allow_pickle=False) as prior:
            dr=saved['ret'][0]-prior['ret'];dl=saved['less'][0]-prior['less']
            np.savez(root/(label+'-complete-findings.npz'),ret_difference=dr,less_difference=dl,raw_ret_diagonal=np.asarray([saved['ret'][0,n,n] for n in range(len(times))]),raw_lesser_diagonal=np.asarray([saved['less'][0,n,n] for n in range(len(times))]))
            assert np.array_equal(saved['ret'][0],prior['ret']) and np.array_equal(saved['less'][0],prior['less'])
            entry=dict(label=label,shape=list(saved['ret'].shape),ret_max=float(np.max(np.abs(dr))),less_max=float(np.max(np.abs(dl))),bitwise_equal=True,old_input=original['input'],old_output=original['output'],old_solution=original['solution'],new_packed_bytes=receipt['provenance']['output']['bytes'],old_packed_bytes=original['output']['bytes'],new_final_bytes=sum((target/(n+'.npy')).stat().st_size for n in ('ret','less')),old_final_bytes=original['solution']['bytes'],findings=binding(root/(label+'-complete-findings.npz')),native=json.loads((target/'.staging'/'slot-000'/'native.json').read_text()),profile=binding(target/'.staging'/'slot-000'/'native.time'))
        close_saved(saved)
        # Reopen and validate durable completion: no packed stage is necessary.
        with FinalBank(target,resume=True) as recovered:
            assert set(recovered.completed)=={0}
        findings.append(entry)
    atomic_json(root/'actual-material-findings.json',dict(findings=findings,scope='exact original paid A58 supplied source/kernels; complete R/L and raw diagonal/state unchanged; not a new selfconsistent/material31 solve'))


def nt31(binary,root):
    times,arrays=free_control(31,2)
    target=root/'Nt31-free-finite-control'
    with FinalBank(target,times,2,1,statistics=-1,order=2) as bank:
        solve_supplied(times,*arrays,statistics=-1,order=2,binary=binary,bank=bank,slot=0,source_bindings=[binding(__file__)])
        bank.finalize()
    saved=load_solution(target)
    dr,dl=saved['ret'][0]-arrays[0],saved['less'][0]-arrays[1]
    np.savez(root/'Nt31-complete-findings.npz',ret_difference=dr,less_difference=dl,times=times)
    assert np.array_equal(saved['ret'][0],arrays[0]) and np.array_equal(saved['less'][0],arrays[1])
    atomic_json(root/'Nt31-findings.json',dict(sourcekind='explicitly-free-finite-Hermitian-two-level-control; zero selfenergy by definition',shape=list(saved['ret'].shape),ret_max=float(np.max(np.abs(dr))),less_max=float(np.max(np.abs(dl))),missing_material31_collision_history=False,scope='arbitrary uniform time count transport only, not original missing material31 collision history',findings=binding(root/'Nt31-complete-findings.npz')))
    close_saved(saved)


if __name__=='__main__':
    binary,root=Path(sys.argv[1]),Path(sys.argv[2])
    boundary(binary,root)
    actual(binary,root)
    nt31(binary,root)
    print('A61_ACTUAL_CONSUMER_COMPLETE',flush=True)
