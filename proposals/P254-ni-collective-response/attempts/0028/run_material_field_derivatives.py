#!/usr/bin/env python3
"""Full-PW selected-column UPF field derivatives and controlled C2 comparison.

No DFT, population, field geometry or damping is chosen. Directions and steps
are explicit numerical probes, not physical fields. Accepted background and
wavefunction bytes come from the existing optical registration. The old
nonconvergent interpolant stencil is not rerun. Exact radial quadrature
coefficients are compared to C2 reciprocal representations at dq,dq/2,dq/4;
only that genuinely differentiable family is finite-differenced.
"""
from pathlib import Path
import argparse
import hashlib
import io
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT/'src', HERE.parent/'0025'):
    sys.path.insert(0,str(directory))
import numpy as np
from substrate_framework.qe_material_field_derivatives import MaterialFieldDerivatives
from substrate_framework.qe_nc_radial import DQ
from material_run_contract import runtime_context


def digest(path):
    path=Path(path); raw=path.read_bytes()
    return dict(path=str(path),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())


def interpolation_diagnostics(hamiltonian, p):
    """Local polynomial slopes/jumps, not a repeat of the rejected optical run."""
    q=np.linalg.norm(p,axis=1)
    result=dict(origin_rows=int(np.count_nonzero(q==0)),
                minimum_distance_to_radial_knot_bohr_inv=float(np.min(np.abs(q/DQ-np.rint(q/DQ)))*DQ),
                species={})
    for species in hamiltonian.factors.species:
        table=hamiltonian.background.inputs.tables[species.name]
        b=table.beta_values
        right=(-11*b[:-3]+18*b[1:-2]-9*b[2:-1]+2*b[3:])/(6*DQ)
        left=(-2*b[:-3]-3*b[1:-2]+6*b[2:-1]-b[3:])/(6*DQ)
        jump=right[1:]-left[:-1]
        row_index=np.floor(q/DQ).astype(int)
        near=np.abs(q/DQ-np.rint(q/DQ))<1e-10
        result['species'][species.name]=dict(
            radial_angular_l=table.upf.angular_l.tolist(),
            forward_cubic_origin_radial_slope=right[0].tolist(),
            max_interior_first_derivative_jump_by_radial_channel=np.max(np.abs(jump),axis=0).tolist(),
            row_interval_indices_range=[int(row_index.min()),int(row_index.max())],
            rows_numerically_on_knots=int(np.count_nonzero(near)),
            explanation='A nonzero l=0 radial slope at q=0 produces a |a| cusp; l>0 also requires regular q^l parity. Nonzero interior jumps preclude a global C1 interpolant. These diagnostics establish obstructions, not exclusive attribution of the prior finite-step errors.')
    return result


def consume(hamiltonian, vectors, *, direction, mixed_direction, step_bohr_inv, spacing_bohr_inv):
    vectors=np.asarray(vectors)
    e=np.asarray(direction,dtype=float); f=np.asarray(mixed_direction,dtype=float)
    h=float(step_bohr_inv); spacing=float(spacing_bohr_inv)
    if vectors.ndim!=2 or not 1<=vectors.shape[1]<=8:
        raise ValueError('one through eight full-row columns required')
    if e.shape!=(3,) or f.shape!=(3,) or not np.isfinite([*e,*f,h,spacing]).all() or h<=0 or spacing<=0:
        raise ValueError('finite directions and positive numerical step/spacing required')
    exact=MaterialFieldDerivatives(hamiltonian)
    p=exact.momenta
    # Exact quadrature and analytic solid harmonics; never finite differences.
    linear,contact=exact.derivative_actions(vectors,e)
    _,mixed=exact.derivative_actions(vectors,e,f)
    h0=hamiltonian.apply(vectors)
    difference=exact.nonlocal_delta_action(vectors,np.zeros(3),anchored=False)
    arrays=dict(input_columns=vectors,direction=e,mixed_direction=f,
                accepted_H0_action_ry=h0,physical_UPF_zero_minus_A25_ry=difference,
                exact_linear_ry_bohr=linear,exact_contact_ry_bohr2=contact,
                exact_mixed_contact_ry_bohr2=mixed,charge_current_au=-linear/2)
    qmax=float(np.max(np.linalg.norm(p,axis=1))+h*(np.linalg.norm(e)+np.linalg.norm(f))+spacing)
    levels=[]
    previous=None
    for level in range(3):
        dq=spacing/2**level
        rep=MaterialFieldDerivatives(hamiltonian,spacing_bohr_inv=dq,qmax_bohr_inv=qmax)
        l,c=rep.derivative_actions(vectors,e)
        _,cm=rep.derivative_actions(vectors,e,f)
        delta=rep.nonlocal_delta_action(vectors,np.zeros(3),anchored=False)
        record=dict(spacing_bohr_inv=dq,linear_error_vs_UPF=float(np.linalg.norm(l-linear)),
                    contact_error_vs_UPF=float(np.linalg.norm(c-contact)),
                    mixed_error_vs_UPF=float(np.linalg.norm(cm-mixed)),
                    zero_field_error_vs_UPF=float(np.linalg.norm(delta-difference)),
                    zero_field_difference_vs_A25=float(np.linalg.norm(delta)))
        if previous is not None:
            record['successive_contact_difference']=float(np.linalg.norm(c-previous))
        previous=c
        levels.append(record)
        arrays[f'representation_{level}_linear']=l
        arrays[f'representation_{level}_contact']=c
        arrays[f'representation_{level}_mixed_contact']=cm
        arrays[f'representation_{level}_zero_minus_A25']=delta
    # This step series belongs ONLY to the new finest C2 family, not A25 cubic.
    finite=[]
    for level in range(3):
        step=h/2**level
        plus=rep.field_increment_action(vectors,step*e)
        minus=rep.field_increment_action(vectors,-step*e)
        direct_l=(plus-minus)/(2*step)
        direct_c=(plus+minus)/step**2
        corners=[rep.field_increment_action(vectors,step*(a*e+b*f)) for a,b in ((1,1),(1,-1),(-1,1),(-1,-1))]
        direct_m=(corners[0]-corners[1]-corners[2]+corners[3])/(4*step**2)
        finite.append(dict(step_bohr_inv=step,linear_error_vs_analytic=float(np.linalg.norm(direct_l-l)),
                           contact_error_vs_analytic=float(np.linalg.norm(direct_c-c)),
                           mixed_error_vs_analytic=float(np.linalg.norm(direct_m-cm)),
                           second_order_Taylor_remainder=float(np.linalg.norm(plus-step*l-.5*step**2*c))))
        arrays[f'C2_finite_{level}_linear']=direct_l
        arrays[f'C2_finite_{level}_contact']=direct_c
        arrays[f'C2_finite_{level}_mixed']=direct_m
    # One finite shift directly compares the quadrature and finest family.
    exact_increment=exact.field_increment_action(vectors,h*e)
    finest_increment=rep.field_increment_action(vectors,h*e)
    arrays['exact_UPF_probe_increment_ry']=exact_increment
    arrays['finest_C2_probe_increment_ry']=finest_increment
    summary=dict(schema='P254-material-field-derivatives-v1',npw=len(p),ncolumn=vectors.shape[1],
                 ispin=hamiltonian.ispin,direction=e.tolist(),mixed_direction=f.tolist(),
                 exact_zero_field_difference_vs_A25_norm=float(np.linalg.norm(difference)),
                 exact_zero_field_difference_relative_to_H0=float(np.linalg.norm(difference)/max(np.linalg.norm(h0),np.finfo(float).tiny)),
                 anchored_zero_increment_norm=float(np.linalg.norm(exact.field_increment_action(vectors,np.zeros(3)))),
                 finite_probe_C2_minus_UPF_norm=float(np.linalg.norm(finest_increment-exact_increment)),
                 representations=levels,C2_finite_difference_refinement=finite,
                 original_interpolation_diagnostics=interpolation_diagnostics(hamiltonian,p),
                 family='H_UPF(a) uses exact Bessel functions on original radial Simpson mesh. H_anchor(a)=H_A25(0)+H_UPF(a)-H_UPF(0). C2 approximations replace only F_l(s) by explicit quintic Hermite interpolation.',
                 current='electron charge current -L/2 in atomic units; no population/k weight/volume applied',
                 reference_approximation='Zero-field A25 H0 is preserved only by a field-independent operator offset. This is not a claim of exact Peierls identity with the saved cubic Hamiltonian.',
                 convergence='No automatic contact convergence claim. Inspect independent oracle, representation errors and new-family step errors; original adverse series is neither rerun nor overridden.',
                 limitations=['Full physical rows but only selected unweighted columns.',
                              'Original radial mesh quadrature retained; no continuum radial error certificate.',
                              'q=0 uniform vector potential only; no physical direction, temperature, linewidth or population is inferred.',
                              'Finite H0 representation discrepancy must be accepted explicitly before using anchored optical coefficients.'])
    return arrays,summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',type=Path,default=HERE/'field-derivative-registration.json')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--first-column',type=int,default=0)
    parser.add_argument('--columns',type=int,default=2)
    parser.add_argument('--direction',nargs=3,type=float,required=True)
    parser.add_argument('--mixed-direction',nargs=3,type=float,required=True)
    parser.add_argument('--step-bohr-inv',type=float,required=True)
    parser.add_argument('--spacing-bohr-inv',type=float,required=True)
    args=parser.parse_args()
    if args.first_column<0 or not 1<=args.columns<=8:
        parser.error('nonnegative first column and 1..8 columns required')
    from substrate_framework.qe_saved_native import read_bound_bytes,load_wavefunctions
    from substrate_framework.qe_material_reference import load_material_inputs
    from substrate_framework.qe_material_hamiltonian import build_pw_factors,MaterialHamiltonian
    from material_reference_output import restore_background
    registration=json.loads(args.registration.read_text())
    for source in registration['original_bytes']:
        read_bound_bytes(source)
    binding=registration['accepted_material']
    metadata_bytes,_=read_bound_bytes(binding['background_metadata'])
    metadata=json.loads(metadata_bytes)
    archive_bytes,_=read_bound_bytes(metadata['archive'])
    with np.load(io.BytesIO(archive_bytes),allow_pickle=False) as archive:
        background_arrays={key:archive[key] for key in archive.files}
    del archive_bytes
    sources=binding['inputs']
    inputs=load_material_inputs(sources['xml'],sources['density'],sources['upfs'])
    background=restore_background(inputs,background_arrays)
    del background_arrays
    wave=load_wavefunctions(binding['wavefunction'])
    stop=args.first_column+args.columns
    if stop>wave.nbnd:
        parser.error('selected columns exceed original saved band count')
    args.output.mkdir(parents=True,exist_ok=False)
    context=runtime_context()
    (args.output/'runtime-context.json').write_text(json.dumps(context,indent=2)+'\n')
    hamiltonian=MaterialHamiltonian(background,build_pw_factors(inputs,wave),wave.ispin,block_size=args.columns)
    arrays,summary=consume(hamiltonian,wave.coefficients[:,args.first_column:stop],
                           direction=args.direction,mixed_direction=args.mixed_direction,
                           step_bohr_inv=args.step_bohr_inv,spacing_bohr_inv=args.spacing_bohr_inv)
    summary.update(original_wavefunction=dict(wave.source),selected_columns_zero_based=list(range(args.first_column,stop)),
                   original_sources=registration['original_bytes'],accepted_material=binding,runtime_context=context,
                   implementation_sources=[digest(p) for p in (Path(__file__),HERE/'run_field_derivative_oracle.py',
                                                               args.registration,ROOT/'src/substrate_framework/qe_material_field_derivatives.py')])
    with (args.output/'field-derivatives.npz').open('xb') as stream:
        np.savez(stream,**arrays)
    with (args.output/'field-derivatives.json').open('x') as stream:
        json.dump(summary,stream,indent=2); stream.write('\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
