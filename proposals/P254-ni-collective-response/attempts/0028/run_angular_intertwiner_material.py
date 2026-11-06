#!/usr/bin/env python3
"""Original Ni/H UPF angular Q for explicit saved unitary cubic operation IDs.

Writes ONLY opINDEX_SPECIES matrices to intertwiners.npz, directly consumable by
run_operator_covariance_material.py --intertwiners/--intertwiners-sha256.
result.json binds the original inputs, channel ledgers, selected operations,
source bytes, independent angular quadrature and both production center checks.
No projector fit, radial channel permutation, D repair or atom/phase adjustment.
"""
import argparse
import json
import os
from pathlib import Path
from time import perf_counter

from run_angular_intertwiner_oracle import (
    ROOT, np, registration, binding, authenticate, sources, publish,
    sphere_quadrature, independent_harmonics, residual, threadpool_limits)
from substrate_framework.material_angular_intertwiners import (
    signed_permutation, original_channel_intertwiner, intertwiner_diagnostics)
from substrate_framework.qe_nc_upf import read_nc_upf
from substrate_framework.qe_nc_radial import build_radial_tables
from substrate_framework.qe_nc_projectors import centered_projectors
from substrate_framework.qe_material_field_derivatives import UPFRadial, centered_derivatives
from substrate_framework.qe_saved_reference import load_saved_reference


def array_binding(value):
    from hashlib import sha256
    value = np.ascontiguousarray(value)
    return dict(shape=list(value.shape), dtype=value.dtype.str,
                sha256=sha256(memoryview(value).cast('B')).hexdigest(),
                encoding='C-order raw array bytes; shape/dtype separately bound')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--operation', required=True, action='append', type=int,
                        help='one-based saved XML operation ID; repeat for exact selected set')
    parser.add_argument('--oracle', required=True, type=Path)
    parser.add_argument('--oracle-sha256', required=True)
    parser.add_argument('--cpu', required=True, type=int, choices=range(6, 11))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    start = perf_counter()
    record = registration()
    if len(set(args.operation)) != len(args.operation):
        parser.error('duplicate original operation ID')
    if args.output.resolve().is_relative_to(ROOT):
        parser.error('material output directory must be outside Git/repository')
    oracle_binding = binding(args.oracle)
    if oracle_binding['sha256'] != args.oracle_sha256:
        raise ValueError('oracle digest mismatch')
    oracle = json.loads(args.oracle.read_text())
    if (oracle.get('schema') != 'P254-angular-intertwiner-oracle-v1'
            or oracle.get('status') != 'passed' or not oracle.get('checks')
            or not all(oracle['checks'].values()) or oracle.get('source_bindings') != sources(record)):
        raise ValueError('passed independent oracle bound to these exact source bytes required')
    input_bindings = record['original_inputs']
    for entry in input_bindings.values():
        authenticate(entry)
    with threadpool_limits(limits=1):
        reference = load_saved_reference(input_bindings['xml'])
        saved = reference.saved_symmetries
        if saved is None or any(not 1 <= i <= saved.nsym for i in args.operation):
            raise ValueError('operation ID absent from original saved crystal group')
        if {s.name for s in reference.species} != {'Ni', 'H'}:
            raise ValueError('original loaded Ni/H species required')
        operations = [saved.operations[i-1] for i in args.operation]
        op_records = []
        reciprocal = reference.reciprocal_bohr_inv
        for operation in operations:
            if operation.effective_time_reversal is not False:
                raise ValueError('antiunitary/unresolved operation cannot supply positive-time covariance Q')
            r = signed_permutation(operation.rotation)
            cartesian = reciprocal @ r @ np.linalg.inv(reciprocal)
            # No rounding arbitrary Cartesian rotations into the supported group.
            # Declare the saved integer S itself; record the floating metric defect.
            # Unsupported non-aligned cells are rejected instead of repaired.
            defect = float(np.max(np.abs(cartesian-r)))
            if defect > 2e-13:
                raise ValueError('saved operation is not the same Cartesian signed permutation in this cell')
            op_records.append(dict(index=operation.index, rotation=r.tolist(),
                cartesian_rotation=cartesian.tolist(), cartesian_declaration_defect_max=defect,
                determinant=int(round(np.linalg.det(r))), info=operation.info,
                fractional_translation=None if operation.fractional_translation is None else operation.fractional_translation.tolist(),
                equivalent_atoms=None if operation.equivalent_atoms is None else operation.equivalent_atoms.tolist(),
                raw_info_attributes=dict(operation.raw_info_attributes),
                raw_rotation_attributes=dict(operation.raw_rotation_attributes),
                raw_time_reversal=operation.raw_time_reversal,
                effective_time_reversal=operation.effective_time_reversal,
                effective_policy=operation.effective_policy))
        upfs = {}
        for species in reference.species:
            entry = input_bindings[species.name]
            if Path(entry['path']).name != species.pseudo_file:
                raise ValueError('saved original species pseudo filename differs')
            upf = read_nc_upf(entry['path'], expected_sha256=entry['sha256'], expected_bytes=entry['bytes'])
            if upf.element != species.name:
                raise ValueError('original UPF element differs from saved species')
            upfs[species.name] = upf
        args.output.mkdir(parents=True, exist_ok=False)
        publish(args.output/'started.json', dict(schema='P254-angular-intertwiner-material-started-v1',
            source_bindings=sources(record), oracle=oracle_binding, original_inputs=input_bindings,
            selected_operation_ids=args.operation, cpu_affinity=sorted(os.sched_getaffinity(0))))
        direction, weight = sphere_quadrature()
        y = independent_harmonics(direction)
        radii = np.array([.17, .63, 1.41])
        p = np.concatenate([radius*direction for radius in radii])
        omega = float(abs(np.linalg.det(reference.cell_bohr)))
        matrices, reports, ledgers = {}, [], {}
        for species, upf in upfs.items():
            tables = build_radial_tables(upf, omega, float(radii.max()))
            radial = UPFRadial(tables)
            direct = centered_derivatives(radial, p, np.zeros(3))[0]
            legacy = centered_projectors(tables, p)
            moduli = np.linalg.norm(p, axis=1)
            independent = (independent_harmonics(p)[:, upf.channel_lm]
                * moduli[:, None]**upf.channel_l
                * radial.evaluate(moduli*moduli)[0][:, upf.channel_radial]
                * (-1j)**upf.channel_l)
            ledgers[species] = dict(angular_l=upf.angular_l.tolist(),
                channel_radial_zero_based=upf.channel_radial.tolist(), channel_l=upf.channel_l.tolist(),
                channel_lm_zero_based=upf.channel_lm.tolist(), nh=upf.nh,
                raw_dij=array_binding(upf.raw_dij), expanded_d=array_binding(upf.expanded_d),
                independent_center_relative_residual=residual(direct, independent),
                radial_identity_policy='each original PP_BETA index retained; no mixing repeated l blocks')
            for operation, op_record in zip(operations, op_records):
                r = np.asarray(op_record['rotation'])
                cartesian = np.asarray(op_record['cartesian_rotation'])
                q = original_channel_intertwiner(upf, r)
                key = 'op'+str(operation.index)+'_'+species
                matrices[key] = q
                yr = independent_harmonics(direction @ cartesian.T)
                integrated = np.zeros_like(q)
                first = 0
                for ell in upf.angular_l:
                    last = first+2*int(ell)+1
                    sl = slice(ell*ell, (ell+1)**2)
                    integrated[first:last, first:last] = y[:, sl].T @ (weight[:, None]*yr[:, sl])
                    first = last
                target_p = p @ cartesian.T
                target_direct = centered_derivatives(radial, target_p, np.zeros(3))[0]
                target_legacy = centered_projectors(tables, target_p)
                origin = centered_derivatives(radial, np.zeros((1, 3)), np.zeros(3))[0]
                angular_defect = yr[:, upf.channel_lm]-y[:, upf.channel_lm] @ q
                reports.append(dict(key=key, operation=operation.index, species=species,
                    matrix=array_binding(q),
                    angular_quadrature_relative_residual=residual(q, integrated),
                    angular_max_residual=float(np.max(np.abs(angular_defect))),
                    direct_UPF_center_relative_residual=residual(target_direct, direct @ q),
                    qe_interpolated_center_relative_residual=residual(target_legacy, legacy @ q),
                    origin_center_relative_residual=residual(origin, origin @ q),
                    direct_UPF_center_norm=float(np.linalg.norm(direct)),
                    qe_interpolated_center_norm=float(np.linalg.norm(legacy)),
                    **intertwiner_diagnostics(q, upf.expanded_d)))
        archive = args.output/'intertwiners.npz'
        with archive.open('xb') as stream:
            np.savez(stream, **matrices)
        expected_keys = {'op'+str(i)+'_'+name for i in args.operation for name in upfs}
        with np.load(archive, allow_pickle=False) as restored:
            roundtrip = (set(restored.files) == expected_keys and
                         all(np.array_equal(restored[key], matrices[key]) for key in expected_keys))
        if not roundtrip:
            raise ValueError('saved NPZ original-channel/key round trip failed')
        archive_binding = binding(archive)
    result = dict(schema='P254-angular-intertwiner-material-v1', status='measured',
        source_bindings=sources(record), oracle=oracle_binding, original_inputs=input_bindings,
        archive=archive_binding, archive_keys=sorted(expected_keys), saved_archive_roundtrip=roundtrip,
        operation_ids=args.operation, operations=op_records, channel_ledgers=ledgers, diagnostics=reports,
        sphere_quadrature=dict(z_order=8, azimuth_points=32, azimuth_offset=.317),
        center_samples=dict(radii_bohr_inv=radii.tolist(), directions=len(direction),
                            volume_bohr3=omega, origin_checked=True),
        cpu_affinity=sorted(os.sched_getaffinity(0)), elapsed_seconds=perf_counter()-start,
        convention='pt=R ps (column Cartesian); Bt=Bs Q, Q in ORIGINAL radial-then-angular channels',
        consumer_arguments=['--intertwiners', str(archive.resolve()), '--intertwiners-sha256', archive_binding['sha256']],
        scope='declared angular intertwiner only; all geometry/phase/projector and active-D defects remain in existing complete covariance bound; no future nonlinear or antiunitary positive-time reuse')
    publish(args.output/'result.json', result)
    print(json.dumps(dict(archive=archive_binding, result=str((args.output/'result.json').resolve()))))


if __name__ == '__main__':
    main()
