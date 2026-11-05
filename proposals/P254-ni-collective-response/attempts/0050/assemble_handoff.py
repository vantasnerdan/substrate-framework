"""Assemble full actual A50 evidence, costs, findings and declared limits."""
import json
from pathlib import Path
import sys
from joint_gamma_adapter import HERE
from joint_material_inputs import file_binding
from moving_reference_inputs import atomic_json


def main(generation):
    base=Path('/var/tmp/p254-a50-joint-gamma-'+generation)
    receipt=Path(str(base)+'-receipts')
    def load(path):return json.loads(Path(path).read_bytes())
    runs=[]
    for path in sorted(base.parent.glob('p254-a50-joint-gamma-*-receipts/exit.json')):
        directory=path.parent
        runs.append(dict(outcome=load(path),manager_startup=load(directory/'manager-startup.json'),
            exact_argv=load(directory/'argv.json'),raw_stdout=file_binding(directory/'stdout'),raw_stderr=file_binding(directory/'stderr')))
    result_path=base/'material/result.json'
    complete=result_path.exists() and (base/'consumer.json').exists()
    files=[file_binding(p) for directory in sorted(base.parent.glob('p254-a50-joint-gamma-*'))
           if directory.is_dir() for p in sorted(directory.rglob('*')) if p.is_file()]
    sources=[file_binding(p) for p in sorted(HERE.glob('*.py'))]
    sources += [file_binding(p) for p in sorted(HERE.glob('*.original-*'))]
    sources += [file_binding(p) for p in sorted(HERE.glob('*.json')) if not p.name.startswith('actual-original-joint-gamma-handoff')]
    inherited=[('0029','reduced_electronic_correlations.py'),('0029','reduced_correlation_transpose.py'),
        ('0029','reduced_electronic_dynamics.py'),('0029','reduced_coefficients.py'),
        ('0029','reduced_projection_storage.py'),('0041','electron_ion_functional.py'),
        ('0041','electron_ion_global_functional.py'),('0041','material_electron_ion_inputs.py'),
        ('0042','durable_active_baseline.py')]
    sources += [file_binding(HERE.parent/a/b) for a,b in inherited]
    root=HERE.parents[3]
    sources += [file_binding(root/'src/substrate_framework'/name) for name in
        ('material_clamped_array.py','material_joint_field_scalar.py','physical_field_environment.py')]
    result=load(result_path) if result_path.exists() else None
    consumer=load(base/'consumer.json') if (base/'consumer.json').exists() else None
    evidence=load(base/'oracle.json') if (base/'oracle.json').exists() else None
    input_path=base/'material/actual-inputs.json'
    handoff=dict(schema='p254-a50-complete-original-machine-handoff-v1',status='complete' if complete else 'incomplete-measured-boundary',
        actual_native_runs=runs,cumulative_native_seconds=sum(r['outcome']['elapsed_seconds'] for r in runs),
        actual_sources=sources,original_inputs=load(input_path) if input_path.exists() else None,
        scientific_result=result,final_numeric_consumer=consumer,new_boundary_oracle=evidence,
        runtime_files=files,runtime_output_bytes=sum(item['bytes'] for item in files),
        failure_dispositions=[load(HERE/name) for name in
            ('actual-first-oracle-failure-001.json','actual-original-oom-002.json')],
        diagnostic_label_disposition='Executed result historical key forward_reverse_gamma_dagger_defect computes ||stored_reverse-forward.T||, NOT a dagger. Numeric value unchanged; future source key corrected to stored_reverse_minus_forward_transpose_norm. Executed source retained byte-for-byte as run_joint_gamma.py.original-before-diagnostic-label-004.',
        numeric_readout_units='Norms are raw canonical-array Frobenius diagnostics, not physical rates or a homogeneous-unit energy scalar. Pi entries retain derivative_i*derivative_j units for original ionic/field coordinates; identical-component new-minus-represented only.',
        runtime_measurement_caveat='Raw systemd terminal Memory peak is retained verbatim; no claim it independently bounds each child peak. Original002 cgroup OOM is explicit, not hidden by final003 success.',
        findings=[
            dict(finding='A41 original ionic-only material assembly omitted Gamma field and mixed coordinates',disposition='joint canonical full114 adapter added; ionic block retained as numerical comparison',evidence='joint_gamma_adapter.py and common.npz consumer'),
            dict(finding='unscaled ionic + unscaled field with retained1/N has incorrect field normalization',disposition='actual derivatives ionic/sqrtN and field_modes; normalized edge N=1',evidence='new scalar SymPy/global-edge oracle'),
            dict(finding='A29 density/K feedback prevents identifying represented channel as undressed A41 bubble',disposition='electron-only tangent with K retained and retarded response independently extracted',evidence='represented_electron_channel.py and represented.npz'),
            dict(finding='response and orderednoise must be removed/replaced together for the same sector',disposition='directed cutover contract and Pi-only identical-unit difference API; no Sigma subtraction',evidence='CUTOVER/polarization_difference and explicit coupled Fock sensitivity oracle'),
            dict(finding='frozen mechanical Hessian/contact is not represented bubble',disposition='full original mechanical generators/j/W0 and all mixed/minimalcoupling/Maxwell blocks unchanged',evidence='original input bindings + joint full SW0ST selection'),
            dict(finding='ab-initio PBE/XC residual matching prescription absent',disposition='unresolved; no invented subtraction or physicalcalibration claim',evidence='source contract + result limits')],
        earned_scope='all111+3 Gamma/all128 CAR/all12 bands/original4times16pairs joint common component and K-dressed represented channel' if complete else 'implementation plus only actually completed persisted boundaries; never promotes a slice to full Gamma',
        quantitative_uncertainty='Finite numerical oracle tolerances and original moving mean rtol2e-8/atol2e-10 retained; represented DOP853 rtol2e-8/atol2e-10. No rigorous enclosure, continuum convergence or material calibration earned.',
        untouched='No A28-48 edits; no protected process inspection; no source/preparation replay, clipping/refill, eigenvalue/symplectic repair or field mass/frequency scaling',
        blockers=[] if complete else ['Actual full material completion absent; inspect exact raw exit and last persisted boundary for measured work, no partial result called full'],
        unresolved_physical=['PBE/XC double-count/correlation residual matching','selfconsistent joint physicalclosure','fullE/physical rates/loss','passive/complement dynamics/noise'])
    atomic_json(HERE/('actual-original-joint-gamma-handoff-'+generation+'.json'),handoff)
    print(json.dumps(dict(status=handoff['status'],native_seconds=handoff['cumulative_native_seconds'],consumer=consumer),allow_nan=False))


if __name__=='__main__':main(sys.argv[1])
