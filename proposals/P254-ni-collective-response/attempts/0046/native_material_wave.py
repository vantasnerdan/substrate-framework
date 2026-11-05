"""Root-authorized WaveII; full original material with measured q-boundary decisions.

No scientific approximation or A45 transition repair. A first-q forecast which
cannot fit the remaining grant pauses only AFTER durable numeric commitment.
"""
import cProfile
import contextlib
import datetime
import io
import json
import os
from pathlib import Path
import pstats
import resource
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
GENERATION = os.environ.get('A46_MATERIAL_GENERATION','001')
OUT = HERE/('material-wave-'+GENERATION)
UNIT = 'p254-a46-moving-material-'+GENERATION+'.service'
PRIOR_SECONDS = float(os.environ.get('A46_MATERIAL_PRIOR_SECONDS','0'))


def save(name,value):
    p = OUT/name
    p.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False))


class PauseAtBound(Exception):
    pass


def main():
    start = time.monotonic()
    OUT.mkdir(exist_ok=True)
    manager_argv = ['systemctl','--user','show',UNIT,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,CPUAffinity,ActiveState,SubState,InvocationID']
    m = subprocess.run(manager_argv,capture_output=True,text=True)
    props = dict(x.split('=',1) for x in m.stdout.splitlines() if '=' in x)
    cg = Path('/sys/fs/cgroup')/props.get('ControlGroup','').lstrip('/')
    kernel = {}
    for name in ('memory.max','memory.swap.max','cpuset.cpus.effective'):
        try:
            kernel[name] = dict(available=True,value=(cg/name).read_text().strip())
        except OSError as e:
            kernel[name] = dict(available=False,reason=str(e))
    save('actual-startup.json',dict(unit=UNIT,manager_argv=manager_argv,
        manager_exit=m.returncode,manager_stdout=m.stdout,manager_stderr=m.stderr,
        properties=props,kernel=kernel,affinity=sorted(os.sched_getaffinity(0)),
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    code,status = 0,'complete'
    profile = cProfile.Profile()
    ordered_calls = []
    def hard_pause(signum,frame):
        raise PauseAtBound('580s owned finite limit; current uncommittedq discarded, durable prefixes retained')
    signal.signal(signal.SIGALRM,hard_pause)
    signal.setitimer(signal.ITIMER_REAL,580)
    with (OUT/'science.rawstdout.txt').open('w',buffering=1) as out, (OUT/'science.rawstderr.txt').open('w',buffering=1) as err:
        with contextlib.redirect_stdout(out),contextlib.redirect_stderr(err):
            try:
                import numpy as np
                import run_all_slot_coulomb as run
                import all_slot_coulomb_contractions as contraction
                from all_slot_coulomb_checkpoint import binding,authenticate
                from material_electron_ion_inputs import SavedTransitions
                if GENERATION != '001':
                    # Only the historical decimal-integer binding schema changed.
                    # Exercise its actual original saved-result/input branch.
                    receipt = json.loads((HERE/'material-wave-001'/'actual-input-admission.json').read_bytes())
                    historical = receipt['coefficient_bindings']['background_receipt']
                    authenticate(historical)
                    integer = dict(historical,bytes=int(historical['bytes']))
                    authenticate(integer)
                    rejected = False
                    try:
                        authenticate(dict(historical,bytes=str(int(historical['bytes'])+1)))
                    except ValueError:
                        rejected = True
                    if not rejected:
                        raise AssertionError('historical byte-count mutation escaped')
                    digest_rejected = False
                    try:
                        authenticate(dict(historical,sha256='0'*64))
                    except ValueError:
                        digest_rejected = True
                    if not digest_rejected:
                        raise AssertionError('historical digest mutation escaped')
                    save('new-lossless-binding-schema-oracle.json',dict(status='passed',
                        original=historical,integer=integer,wrong_byte_count_rejected=rejected,
                        wrong_digest_rejected=digest_rejected,
                        existing_physical_oracles='carried, not replayed'))
                c = run.load_coulomb_inputs(density_cap_bytes=900*1024**2)
                times = np.array([0.,1.,10.,30.])
                bundle = '/var/tmp/p254-a45-moving-reference-001/producer/bundle.json'
                reference = SavedTransitions(c,times,bundle,c.bindings['manifest'])
                save('actual-input-admission.json',dict(coefficient_bindings=c.bindings,
                    actual_reference_bindings=reference.bindings,
                    approximation=json.loads(Path(bundle).read_bytes())['approximation'],
                    A45_unrepaired_diagnostics=dict(U_unitarity_frobenius=2.3984622285716737e-6,
                        rho_reconstruction_frobenius=1.4640124062085698e-6,
                        rho_reconstruction_maxabs=6.749818615148229e-8),
                    accuracy='diagnostics only, NOT rigorous global accuracy certification',
                    feature_counts=[len(g) for g in c.feature_g]))
                original_ordered = contraction.ordered_q
                def timed_ordered(factors,q,*args,**kwargs):
                    ts = time.monotonic()
                    result = original_ordered(factors,q,*args,**kwargs)
                    event = dict(q=q,ordered_call=len(ordered_calls)+1,seconds=time.monotonic()-ts,
                        cumulative_native_seconds=time.monotonic()-start)
                    ordered_calls.append(event)
                    save('ordered-call-costs.json',ordered_calls)
                    print(json.dumps(event),flush=True)
                    return result
                contraction.ordered_q = timed_ordered
                original_evaluate = run.evaluate_time_grid
                eval_start = time.monotonic()
                last_q_boundary = [eval_start]
                def measured_evaluate(factors,electrons,times,checkpoint,**kwargs):
                    callback = kwargs.pop('after_commit',None)
                    def committed(q):
                        elapsed = time.monotonic()-start
                        q_seconds = sum(x['seconds'] for x in ordered_calls if x['q']==q-1)
                        if not ordered_calls:
                            q_seconds = time.monotonic()-last_q_boundary[0]
                        last_q_boundary[0] = time.monotonic()
                        remaining_seconds = max(0.,3600-PRIOR_SECONDS-elapsed)
                        feature_counts = [len(g) for g in c.feature_g]
                        predicted_remaining = q_seconds*sum(feature_counts[q:])/feature_counts[q-1]
                        decision = dict(completed_q=q,all128externals=True,all16pairs=True,
                            complete_current_q_seconds=q_seconds,cumulative_native_seconds=elapsed,
                            remaining_grant_seconds=remaining_seconds,
                            predicted_remaining_full_result_seconds=predicted_remaining,
                            forecast_model='dominant exchange q-feature width * fixed full inner-h work; direct included, no promise',
                            action='pause-after-durable-prefix' if predicted_remaining>remaining_seconds else 'continue')
                        save('allocation-decision-q-%02d.json'%q,decision)
                        print(json.dumps(decision),flush=True)
                        if callback is not None:
                            callback(q)
                        if predicted_remaining>remaining_seconds:
                            raise PauseAtBound('first completedq predicts full result outside remaining3600s; checkpoint-safe pause')
                        if elapsed>520:
                            raise PauseAtBound('initial600s Native cap approaching; pause after durableq')
                    return original_evaluate(factors,electrons,times,checkpoint,
                        after_commit=committed,**kwargs)
                run.evaluate_time_grid = measured_evaluate
                profile.enable()
                try:
                    result = run.run_reference(c,reference,times,OUT/'result',external_batch=4,
                        extra_sources=[binding(__file__)])
                    save('actual-complete-consumption.json',result)
                finally:
                    profile.disable()
                    reference.close()
            except PauseAtBound as e:
                code,status = 2,'checkpoint-safe-paused'
                print(str(e),file=sys.stderr,flush=True)
            except BaseException:
                code,status = 1,'failed'
                traceback.print_exc()
            finally:
                signal.setitimer(signal.ITIMER_REAL,0)
                profile.disable()
                profile.dump_stats(str(OUT/'actual-contractions.profile'))
                text = io.StringIO()
                pstats.Stats(profile,stream=text).sort_stats('cumulative').print_stats(35)
                (OUT/'actual-contractions-profile.txt').write_text(text.getvalue())
    terminal = dict(unit=UNIT,status=status,exit=code,seconds=time.monotonic()-start,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        ordered_calls_completed=len(ordered_calls),end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    save('actual-native-exit.json',terminal)
    print(json.dumps(terminal),flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
