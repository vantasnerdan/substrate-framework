"""Resume paid actual q0 prefix with original direct q1..63 contractions."""
from refinement_common import *
import time
from refine_cubic import (_checkpoint,_consume,KEYS,CubicContour,load_material,SOURCE,FINALIZATION_RESERVE)
from cubic_contour import ordered_matrix

def run(output,wall_cap):
    start=time.monotonic();output=Path(output);output.mkdir(parents=True,exist_ok=False)
    initial=HERE/'material/cubic'
    old_scope=json.loads((initial/'scope.json').read_bytes())
    source_binding=binding(initial/'scope.json');archive=old_scope['archive']
    if old_scope['status']!='partial' or old_scope['completed_internal_q']!=[0]:raise ValueError('exact paid q0 prefix required')
    with np.load(authenticate(archive),allow_pickle=False) as z:
        arrays={key:z[key] for key in KEYS}
        if not np.array_equal(z['completed_internal_q'],[0]) or not np.array_equal(z['new_time_pairs'],NEW_PAIRS):raise ValueError('paid q0 mask/coverage differs')
        counts=z['pair_internal_q_counts']
        if any(int(counts[a,b])!=1 for a,b in NEW_PAIRS):raise ValueError('genuine allfive paid q0 contributions required')
    c,manifest,coefficient_binding=load_material();provider=CubicContour.from_source(SOURCE);ref=reference(c,coefficient_binding)
    inputs=dict(old_scope['inputs'],paid_initial_scope=source_binding,paid_initial_archive=archive,continuation_grant=binding(HERE/'direct-cubic-continuation-grant.json'),continuation_implementation=binding(__file__),initial_source_unchanged=binding(HERE/'refine_cubic.py'))
    timings=list(old_scope['timings']);completed=[0];decisions=[]
    try:
        for q in range(1,64):
            unit=time.monotonic();l,left,right=provider.pair(0,q);vertex_seconds=time.monotonic()-unit
            for a,b in NEW_PAIRS:
                dg1,dl1=ref.boson(q,a,b);dg2,dl2=(dg1,dl1) if l==q else ref.boson(l,a,b)
                for key,first,second in (('greater',dg1,dg2),('lesser',dl1,dl2)):
                    arrays[key][a,b]+=ordered_matrix(left,right,first,second)
                arrays['retarded'][a,b]=(arrays['greater'][a,b]-arrays['lesser'][a,b]) if TIMES[a]>=TIMES[b] else 0.
            completed.append(q);seconds=time.monotonic()-unit
            timings.append(dict(internal_q=q,routed_partner_q=int(l),seconds=seconds,vertex_seconds=vertex_seconds,new_time_pair_count=5,old_time_pair_count=0))
            scope=_checkpoint(output,arrays,completed,timings,inputs)
            elapsed=time.monotonic()-start;peak=max(t['seconds'] for t in timings[1:]);forecast=elapsed+(64-len(completed))*peak+FINALIZATION_RESERVE
            print(json.dumps(dict(direct_completed_q=q,reused_q0=True,seconds=seconds,elapsed_seconds=elapsed,forecast_seconds=forecast)),flush=True)
            if q==1:
                decision=dict(first_continued_q=timings[-1],elapsed_seconds=elapsed,forecast_seconds=forecast,wall_cap_seconds=wall_cap,action='continue original direct algorithm within amended grant',paid_q0_recomputed=False)
                save(output/'first-continued-q-decision.json',decision);decisions.append(decision)
            if len(completed)<64 and elapsed+max(1.5*peak,10.)+FINALIZATION_RESERVE>wall_cap:
                decisions.append(dict(action='preserve genuine direct partial before next q',elapsed_seconds=elapsed,completed_internal_q=completed));break
    finally:ref.close()
    consumer=_consume(output,old_scope['inputs']['primary_old_archive'])
    result=dict(schema='p254-a55-direct-cubic-continuation-v1',status=scope['status'],archive=scope['archive'],scope=binding(output/'scope.json'),actual_saved_consumer=binding(output/'saved-consumer.json'),consumer=consumer,inputs=inputs,completed_internal_q=completed,required_internal_q=list(range(64)),reused_internal_q=[0],new_internal_q=completed[1:],initial_paid_q0_archive=archive,initial_paid_q0_scope=source_binding,paid_q0_recomputed=False,old_pairs_computed_per_q=0,new_time_pairs=NEW_PAIRS,times_inverse_ry=TIMES.tolist(),timings=timings,decisions=decisions,elapsed_seconds=time.monotonic()-start,limits=['Original direct all64q full111Gamma only; no field cubic invention','Five new timepairs; oldfour final paid values literal','Not coupled evolution or T30/globalaccuracy/physical calibration/rates/fullE'])
    save(output/'result.json',result);return result
