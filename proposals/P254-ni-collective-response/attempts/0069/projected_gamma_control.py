"""Original represented coupled Gaussian Ni1/x Gamma control, one probe only.

The A29 true bilinear adjoint and original product CAR/W0 contraction own the
calculation.  Covariance[a,b] is <delta Q(t_a) delta Q(t_b)> for Q=u_Gamma,0/sqrt(N).
No mean is added; this is not a matched evolved-mean control or intrinsic closure.
Complex JSON scalars are {real,imag}; NPZ retains the unmodified complex128 data.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/'0050')]
from joint_material_inputs import load_inputs
from material_electron_ion_inputs import authenticate, file_binding
from reduced_electronic_correlations import ActiveBandCorrelations
from reduced_coherence_readout import normalized_ordered


def jsonable(value):
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=float(value.real), imag=float(value.imag))
    if isinstance(value, np.generic):
        return value.item()
    return value


def save_json(path, record):
    with Path(path).open('x') as stream:
        json.dump(jsonable(record), stream, indent=2, allow_nan=False)
        stream.write('\n')


def source_compatibility(c, baseline, bindings, reference_path):
    """Compare actual source archives/layout/baseline, not descriptive labels."""
    reference_binding = file_binding(reference_path)
    frozen = json.loads((HERE/'preregistration.json').read_bytes())
    declared = frozen['verified_raw_input_bindings'][0]
    if reference_binding != declared:
        raise ValueError('reference bytes differ from frozen A52 binding')
    reference = json.loads(authenticate(reference_binding).read_bytes())
    if reference['schema'] != 'p254-a52-requested-time-reference-v1':
        raise ValueError('actual original A52 reference required')
    if reference['coefficient_manifest'] != bindings['manifest']:
        raise ValueError('A52 and coupled control coefficient manifest mismatch')
    if baseline.binding not in reference['physical_history_bindings']:
        raise ValueError('A52 reference does not bind this actual coupled baseline')
    for entry in reference['physical_history_bindings']:
        authenticate(entry)
    if not all(t in reference['times_inverse_ry'] for t in (0., 30.)):
        raise ValueError('genuine endpoint reference absent')
    spec = baseline.receipt['coefficient_specification']
    for key in ('k_indices', 'q_indices', 'slot_spin', 'slot_k', 'target', 'minus'):
        if not np.array_equal(getattr(c, key), np.asarray(spec['indexing'][key])):
            raise ValueError('actual original baseline layout mismatch: '+key)
    snapshot_binding = baseline.receipt['source']['original_coefficient_manifest_snapshot']
    snapshot = json.loads(authenticate(snapshot_binding).read_bytes())
    compared = []
    for key, item in bindings['consumed_items'].items():
        if item['archive'] != snapshot['items'][key]['archive']:
            raise ValueError('original baseline/source archive mismatch: '+key)
        compared.append(key)
    if (c.ns,c.n,c.d,c.f,c.m,c.nspin) != (128,64,111,3,12,2):
        raise ValueError('original represented scope differs')
    if int(c.minus[0]) != 0 or not np.array_equal(c.target[0],np.arange(c.ns)):
        raise ValueError('Gamma self-pairing source map differs')
    block = c.q(0)
    if not np.any(block.kernel) or not np.any(block.density) or not np.any(block.potential):
        raise ValueError('required original density/K response is absent')
    return dict(compatible=True, reference=reference_binding, baseline=baseline.binding,
        original_snapshot=snapshot_binding, exact_archive_matches=compared,
        layout_matches=['k_indices','q_indices','slot_spin','slot_k','target','minus'],
        density_kernel_potential_shapes=[list(v.shape) for v in (block.density,block.kernel,block.potential)],
        same_model='A29 coupled tangent about the exact A42 moving baseline used by A52; A50 represented electron density/K channel retained, original mechanical/field generators and initial W0/rho0 unchanged',
        distinction='A52 bare S supplies a reference only; this control uses coupled electron feedback, not bare S or undressed A41 polarization')


def produce(output, reference_path='/var/tmp/p254-a52-requested-reference-001/producer/bundle.json'):
    output = Path(output)
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    c, baseline, _old_reference, bindings = load_inputs()
    compatibility = source_compatibility(c,baseline,bindings,reference_path)
    print('A69_SOURCE_COMPATIBLE '+str(time.monotonic()-started),flush=True)
    model = ActiveBandCorrelations(c,baseline,rtol=2e-8,atol=2e-10,method='DOP853')
    probe = model.ionic_probes(0,(0,))/c.sqrt_n
    initial = model.backward_linear_probes(0,0.,probe)
    initial_variance = model.initial_ordered(0,initial,initial).item()
    expected_initial = c.mechanical.w0[0][0,0]/c.n
    if not np.allclose(initial_variance,expected_initial,rtol=4*np.finfo(float).eps,atol=0.):
        raise ValueError('canonical initial projected variance differs from bound W0/N')
    save_json(output/'source-compatibility.json',compatibility)
    pulled = model.backward_linear_probes(0,30.,probe)
    duals = (initial,pulled)
    covariance = np.asarray([[model.initial_ordered(0,left,right).item()
                              for right in duals] for left in duals],dtype=np.complex128)
    # Independent ordered contractions: no conjugation or Hermitian repair.
    raw = output/'projected-covariance.npz'
    with raw.open('xb') as stream:
        np.savez(stream,covariance=covariance,times=np.array([0.,30.]))
    dual_path = output/'projected-pullbacks.npy'
    with dual_path.open('xb') as stream:
        np.save(stream,np.stack(duals),allow_pickle=False)
    sources = sorted({Path(module.__file__).resolve() for module in list(sys.modules.values())
        if getattr(module,'__file__',None) and str(HERE.parent) in str(module.__file__)
        and str(module.__file__).endswith('.py')})
    result = dict(schema='p254-projected-gamma-control-v1',status='complete-at-stated-scope',
        kind='original-represented-coupled-Gaussian',reference=compatibility['reference'],
        baseline=baseline.binding,probe=dict(name='canonicalGammaNi1/x',index=0,
        normalization='1sqrtN',coefficient=float(1/c.sqrt_n),q=0),n=int(c.n),times=[0.,30.],
        dimensions=dict(slots=c.ns,bands=c.m,spins=c.nspin,ionic=c.d,field=c.f,
                        joint_positions=c.d+c.f,coupled_state=model.state_size(0),propagated_probe_columns=1),
        raw_archive=file_binding(raw),raw_complex_covariance=covariance,
        pullback_archive=file_binding(dual_path),source_compatibility=file_binding(output/'source-compatibility.json'),
        normalized_ordered=normalized_ordered(covariance[1,0],covariance[1,1],covariance[0,0]),
        initial_preparation=dict(expected_W0_over_N=expected_initial,actual_variance=initial_variance,
            matched=True,mean_added=False,preparation='literal original product W0 and full128 coherent CAR rho0'),
        settings=dict(method='DOP853',rtol=model.rtol,atol=model.atol,precision='complex128',
            temporal_error='unresolved beyond original tolerance-controlled adaptive solve; no refinement claim'),
        carried_proof=dict(binding=file_binding(HERE.parent/'0050/actual-original-joint-gamma-final-binding-004.json'),
            role='Existing A50 source-specific Schur/CAR/true transpose evidence carried unchanged; no oracle replay'),
        input_bindings=bindings,source_bindings=[file_binding(path) for path in sources],
        preregistration=file_binding(HERE/'preregistration.json'),
        findings=[dict(status='established-at-stated-scope',claim='genuine endpoint connected projected covariance from full original represented coupled tangent with density/K and mechanical/field feedback'),
            dict(status='carried',claim='existing source-specific Schur/CAR proof unchanged'),
            dict(status='unresolved',claim='timestep/refinement error, intrinsic full dynamics, matched evolved mean, realspace q sum/ranking, physical rate/heat are not determined')],
        blockers=[],native_elapsed_seconds=time.monotonic()-started)
    save_json(output/'projected-control.json',result)
    return output/'projected-control.json'


def consume(path):
    """Actual saved-path normalization; raw archive is the only covariance input."""
    record = json.loads(Path(path).read_bytes())
    if (record['schema'],record['status'],record['kind']) != (
        'p254-projected-gamma-control-v1','complete-at-stated-scope','original-represented-coupled-Gaussian'):
        raise ValueError('complete original projected control required')
    if record['probe']['index'] != 0 or record['probe']['normalization'] != '1sqrtN' or record['times'] != [0.,30.]:
        raise ValueError('required canonical endpoint probe differs')
    with np.load(authenticate(record['raw_archive']),allow_pickle=False) as archive:
        covariance = archive['covariance']
        if covariance.shape != (2,2) or covariance.dtype != np.complex128 or not np.array_equal(archive['times'],[0.,30.]):
            raise ValueError('exact raw complex endpoint covariance required')
        readout = normalized_ordered(covariance[1,0],covariance[1,1],covariance[0,0])
    if jsonable(readout) != record['normalized_ordered']:
        raise ValueError('saved-path normalized readout differs')
    return dict(schema='p254-projected-control-saved-path-exercise-v1',control=file_binding(path),
        probe=record['probe'],times=record['times'],n=record['n'],normalized_ordered=readout,
        raw_complex_covariance=covariance,eligible=True,
        scope='actual persisted covariance read and normalized through existing canonical A29 readout')
