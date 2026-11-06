"""Read existing failed-generation recipes; no H0 reconstruction or response replay."""
import hashlib, io, json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
report_path=HERE/'pilot-001/result.json'
raw=report_path.read_bytes()
report=json.loads(raw)
totals={branch:np.zeros((7,8,8),complex) for branch in ('plus','minus')}
expected=np.zeros(2)
rows=[]
for name,receipt in report['committed_recipes'].items():
    source=receipt['archive']; data=Path(source['path']).read_bytes()
    if hashlib.sha256(data).hexdigest()!=source['sha256'] or str(len(data))!=source['bytes']:
        raise ValueError('changed original archive')
    with np.load(io.BytesIO(data),allow_pickle=False) as saved:
        branch=receipt['metadata']['branch']; identity=receipt['metadata']['identity']
        spin=int(identity.split(':')[0].removeprefix('spin'))-1
        f=saved['source_occupations']; ft=saved['target_occupations']; w=float(saved['weight'])
        es=saved['source_energies']; et=saved['target_energies']
        operators=saved['operators']; number=operators[6+spin]
        correlation=saved['ordered_correlation']; times=saved['times_tau_ry']
        totals[branch]+=correlation
        variance=float(w*np.sum(f*(1-f)))
        if branch=='plus': expected[spin]+=variance
        norm=np.sum(np.abs(number)**2,axis=0)
        mappedf=np.sum(np.abs(number)**2*ft[:,None],axis=0)
        rows.append({'name':name,'identity':identity,'branch':branch,'archive':source,
          'raw_variance':variance,'saved_number_correlation_t0':float(correlation[3,6+spin,6+spin].real),
          'density_column_norm_maxabs_error':float(np.max(np.abs(norm-1))),
          'target_mapped_probability_maxabs_error':float(np.max(np.abs(mappedf-f))),
          'number_energy_action_maxabs':float(np.max(np.abs(number*(et[:,None]-es[None,:])))),
          'target_dimension':str(et.size),'source_dimension':str(es.size)})
reverse=[int(np.flatnonzero(times==-tau).item()) for tau in times]
comm=totals['plus']-totals['minus'][reverse].swapaxes(1,2)
result={'source_result':{'path':str(report_path),'bytes':str(len(raw)),'sha256':hashlib.sha256(raw).hexdigest()},
 'expected_normalized_spin_variance':expected.tolist(),
 'observed_plus_number_t0':totals['plus'][3,6:,6:].real.tolist(),
 'number_covariance_maxabs_error':float(np.max(np.abs(totals['plus'][:,6:,6:]-np.diag(expected)))),
 'number_commutator_maxabs':float(max(np.max(np.abs(comm[:,6:,:])),np.max(np.abs(comm[:,:,6:])))),
 'rows':rows,'scope':'Existing immutable recipe/output inspection only; no new frame/Native calculation or relaxed allowance.'}
(HERE/'gamma-failure-diagnostic-001.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({key:value for key,value in result.items() if key!='rows'},indent=2))
print(json.dumps(sorted(rows,key=lambda x:x['density_column_norm_maxabs_error'],reverse=True)[:2],indent=2))
