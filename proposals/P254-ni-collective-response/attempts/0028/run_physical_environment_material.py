#!/usr/bin/env python3
"""Authenticate completed initial force noise and admitted-cell q0 field units.

This consumer does NOT calculate new force actions, a material trajectory,
finite-q photons, a field initial state, or unequal-time material data.
It preserves the actual unweighted spin1/sourceik2 checkpoint and tests its
ordered covariance through the same interface used by future two-time data.
Root chooses the finite periodic replication count explicitly; a unit-cell
normalization is always also reported. No cavity or external pump is assumed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
for directory in (ROOT/'src',HERE.parent/'0025',HERE.parent/'0024'):
    sys.path.insert(0,str(directory))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--periodic-cells',type=int,required=True,
                        help='declared finite periodic replication count, not a fitted mode volume')
    parser.add_argument('--output-max-bytes',type=int,required=True)
    args=parser.parse_args()
    target=args.output.resolve(); target.mkdir(parents=True,exist_ok=False)
    began=time.monotonic()
    report=dict(schema='P254-physical-environment-material-v1',status='failed',
                scope='actual initial equal-time checkpoint and Maxwell normalization ONLY')
    try:
        import numpy as np
        from material_run_contract import binding,bound_json,recover_block,runtime_context
        from material_reference_output import load_arrays
        from substrate_framework.material_bloch_output import preload_parsers
        from substrate_framework.qe_saved_reference import load_saved_reference
        from substrate_framework import physical_field_environment as field_module
        from substrate_framework import nonstationary_force_kernel as kernel_module
        if args.periodic_cells<1 or args.output_max_bytes<65536:
            raise ValueError('positive explicit periodic cell count and >=65536 byte output capacity required')
        registration_path=HERE/'physical-environment-registration.json'
        registration_binding=binding(registration_path)
        registration=bound_json(registration_binding)
        for source in registration['source_bindings']:
            current=binding(source['path'])
            if current['sha256']!=source['sha256'] or int(current['bytes'])!=source['size_bytes']:
                raise ValueError('registered source bytes changed: '+source['path'])
        checkpoint_path=HERE/'material-force-noise-001/force-noise.json'
        checkpoint_binding=binding(checkpoint_path)
        receipt=bound_json(checkpoint_binding)
        archive=Path(receipt['archive']['path'])
        recovered=recover_block(archive.parent,archive.stem,receipt['recipe_sha256'],receipt['metadata'])
        if recovered!=receipt:
            raise ValueError('completed initial force-noise archive/array provenance differs')
        arrays=load_arrays(receipt)
        metadata=receipt['metadata']
        if (metadata['source_checkpoint_metadata']['spin'],
                metadata['source_checkpoint_metadata']['source_ik'])!=(1,2):
            raise ValueError('registered actual noise source identity differs')
        prep_binding=metadata['recipe']['preparation_terminal']
        prep=bound_json(prep_binding)
        if prep.get('status')!='complete':
            raise ValueError('completed admitted preparation required')
        bloch_binding=prep['recipe']['inputs']['bloch_terminal']
        bloch=bound_json(bloch_binding)
        if bloch.get('status')!='complete':
            raise ValueError('completed admitted Bloch source required')
        xml_binding=bloch['recipe']['inputs']['xml']
        preload_parsers()
        reference=load_saved_reference(xml_binding)
        authority='admitted-XML:'+xml_binding['sha256']
        cell_field=field_module.UniformMaxwellField.from_cell(reference.cell_bohr,
            polarizations=np.eye(3),cell_authority=authority)
        finite_field=field_module.UniformMaxwellField(
            cell_field.volume_bohr3*args.periodic_cells,np.eye(3),
            authority+'; explicit periodic copies='+str(args.periodic_cells))
        kernel=kernel_module.NonstationaryForceKernel.from_initial_checkpoint(arrays,metadata,
            checkpoint_authority=checkpoint_binding['sha256']+':'+receipt['archive']['sha256'])
        def packed(value):
            value=np.asarray(value)
            return dict(real=value.real.tolist(),imag=value.imag.tolist())
        pairs=[]
        for left in kernel.initial_names:
            for right in kernel.initial_names:
                pair=kernel.initial_pair(left,right)
                pairs.append(dict(left=left,right=right,ordered_AB=packed(pair.ordered_ab),
                    ordered_BA=packed(pair.ordered_ba),symmetrized=packed(pair.symmetrized),
                    commutator=packed(pair.commutator),
                    retarded_minus_source_equal_time_theta_half=packed(
                        pair.retarded(perturbation_sign=-1,equal_time_theta=0.5))))
        source_paths=[Path(__file__),Path(field_module.__file__),Path(kernel_module.__file__),registration_path]
        implementation_sources=[binding(p) for p in source_paths]
        # Expose practical arithmetic/storage counts without performing any new
        # full-PW physics. Actual source action files were temporary and are NOT
        # reclassified as available propagated material factors.
        rows,rank=metadata['full_pw_dimension'],metadata['factor_count']
        report.update(status='complete',runtime_context=runtime_context(),
            registration=registration_binding,checkpoint=checkpoint_binding,
            authenticated_archive=receipt['archive'],preparation_terminal=prep_binding,
            bloch_terminal=bloch_binding,admitted_xml=xml_binding,
            implementation_sources=implementation_sources,
            unit_cell_field=cell_field.normalization(),finite_periodic_field=finite_field.normalization(),
            periodic_cells=args.periodic_cells,
            actual_data=dict(time_inverse_ry=0.,operator_order=list(kernel.initial_names),
                operator_units={name:metadata['operators'][name]['units'] for name in kernel.initial_names},
                means=packed(kernel.initial_means),ordered_covariance=packed(kernel.initial_covariance),
                ordered_pairs=pairs,force_actions_recomputed=False,
                covariance_source='stored completed force-noise.npz; validated, not replaced',
                full_physical_pw_rows=rows,all_raw_columns=rank,
                source_slot_measure=metadata['source_slot_measure'],source_slot_measure_applied=False,
                full_material_total=False),
            future_two_time_data=dict(available=False,actual_times=[],
                required=['same prepared D and initial raw frame identity',
                    'U(t,0)^dagger F(t) U(t,0) Psi and independent adjoint actions, all PW rows/columns',
                    'shared authenticated quadratic propagation history and absolute times',
                    'explicit complete independent sector coverage and observable normalization'],
                response_scope='bare quasifree quantum kernel; screened moving-mean response requires PreparedBlochResponse',
                no_inference=['FDT','temperature','stationarity','irreversible loss','finite-q photon vertices']),
            missing_for_joint_actual_trajectory=[
                'declared initial canonical field coordinate/momentum or quantum field state',
                'current/contact from same unanchored physical UPF scalar along actual trajectory',
                'field-dependent ionic/electronic derivatives of that SAME scalar',
                'full-source/slot propagation and requested time/readout grid'],
            work_storage=dict(actual_input_archive_bytes=int(receipt['archive']['bytes']),
                new_full_PW_operator_applications=0,
                initial_contractions='m^2 small rank matrices; O(m^2*r^3), no PW rows loaded',
                future_one_action_complex128_bytes=16*rows*rank,
                future_independent_action_adjoint_pair_bytes=32*rows*rank,
                future_one_kernel_pair_work='O(P*r^2+r^3), row-blocked; no P-square matrix',
                future_all_times_storage='at most 32*n_times*n_probes*P*r bytes; replay/stream to bound storage'),
            energy_scope='finite conservative periodic q0 environment; no fitted damping or fabricated two-time data')
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc())
    report['elapsed_seconds']=time.monotonic()-began
    text=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(text.encode())>args.output_max_bytes:
        raise RuntimeError('physical-environment report exceeds granted output capacity')
    with (target/'summary.json').open('x') as stream:
        stream.write(text)
    print(json.dumps(dict(status=report['status'],summary=str(target/'summary.json'))),flush=True)
    return 0 if report['status']=='complete' else 1


if __name__=='__main__':
    raise SystemExit(main())
