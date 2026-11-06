#!/usr/bin/env python3
"""Independent SciPy/SymPy equations, multi-case weights and saved-output smoke.

Optional --material-slot consumes a real completed target receipt; absence is
reported explicitly and is not material evidence. No source-text/wiring tests.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace
from matched_reference_io import (binding, save_json, save_arrays, runtime,
                                  registration, receipt_arrays)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--array-root', type=Path, required=True)
    parser.add_argument('--material-slot', type=Path)
    args = parser.parse_args()
    runtime(args.array_root)
    import numpy as np
    import sympy as sy
    from scipy.linalg import eigh
    from substrate_framework.material_reference_state import prepare_source
    from substrate_framework.material_reference_transport import raw_state_transport
    from substrate_framework.material_physical_field_moments import overlap_moments
    from substrate_framework.material_mixed_field_ionic import mixed_overlap_coefficient
    from substrate_framework.material_matched_reference import (aggregate_complete, assemble_frozen,
                                                               slot_identity, KEYS)
    from substrate_framework.material_force_lattice import force_lattice_to_q_grid
    reg_binding, reg = registration()
    rng = np.random.default_rng(254903)
    errors, sensitivity, saved = {}, {}, {}
    worst_anti, worst_wrong, worst_weight = 0., 0., 0.
    for case in range(4):
        # Deliberately nonorthogonal complex frame; independent SciPy spectral map.
        psi = rng.normal(size=(5, 3))+1j*rng.normal(size=(5, 3))
        pop = np.array([.91, .43, .08])/(1+case/5)
        state = prepare_source(psi, pop, source_identity='synthetic-'+str(case),
                               basis_identity='independent', authority='oracle literals')
        values, vectors = eigh(psi.conj().T @ psi)
        inverse = (vectors/np.sqrt(values)) @ vectors.conj().T
        d = inverse @ np.diag(pop) @ inverse
        errors['density-map-'+str(case)] = float(np.max(abs(d-state.density_coefficients)))
        slots, records, expected = [], [], None
        for j, weight in enumerate((.17, .29, .54)):
            anti = j != 0
            spin = 1 if j < 2 else 2
            slot = SimpleNamespace(ispin=spin, mesh_index=(j, 0, 0), source_ik=case+1,
                source_binding=dict(path='synthetic', bytes='0', sha256=str(case)), weight=weight,
                operation=SimpleNamespace(index=1), orbital_conjugation=anti,
                reciprocal_shift=np.zeros(3, int))
            slots.append(slot)
            # U is a genuine complex unitary, not a real fixture hiding conjugation.
            q, _ = np.linalg.qr(rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)))
            frame = q @ (psi.conj() if anti else psi)
            target = raw_state_transport(state, antiunitary=anti)
            rho = frame @ target.density_coefficients @ frame.conj().T
            reference_rho = q @ ((psi @ d @ psi.conj().T).conj() if anti else
                                  (psi @ d @ psi.conj().T)) @ q.conj().T
            worst_anti = max(worst_anti, float(np.max(abs(rho-reference_rho))))
            if anti:
                wrong = frame @ state.density_coefficients @ frame.conj().T
                worst_wrong = max(worst_wrong, float(np.max(abs(wrong-reference_rho))))
            # Independent physical target observable, no assumed moment covariance.
            h = rng.normal(size=(5,5))+1j*rng.normal(size=(5,5))
            h = h+h.conj().T
            value = np.trace(rho @ h)
            arrays = dict(density_bohr3=np.diag(rho), positive_gradient_ry_bohr=np.full(3,value),
                local_onsite_ry_bohr2=np.full((1,3,3), value/3),
                nonlocal_onsite_ry_bohr2=np.full((1,3,3), value/7),
                field_gradient=np.array([value,2*value,-value]),
                current_hartree=np.array([-value/2,-value,value/2]),
                contact=np.eye(3)*2*value, mixed_ionic_field=np.full((1,3,3),value/5),
                number=np.asarray(np.trace(rho)))
            records.append((slot_identity(slot), arrays))
            if expected is None:
                expected = {k: np.zeros_like(v) for k,v in arrays.items() if k!='density_bohr3'}
                expected['density_bohr3'] = np.zeros((5,2),complex)
            for k,v in arrays.items():
                if k=='density_bohr3':
                    expected[k][:,spin-1] += weight*np.diag(reference_rho)
                else:
                    # Dense trace independently evaluated rather than aggregate output.
                    ratio = np.trace(reference_rho @ h)/value if k!='number' else 1
                    expected[k] += weight*(v*ratio if k!='number' else np.trace(reference_rho))
        records.sort(key=lambda row:(row[0]['spin'],tuple(row[0]['mesh_index'])))
        grid = SimpleNamespace(slots=slots)
        actual = aggregate_complete(grid, records)
        errors['multi-case-weight-'+str(case)] = max(float(np.max(abs(actual[k]-expected[k]))) for k in KEYS)
        worst_weight = max(worst_weight, float(abs(sum(r[1]['number'] for r in records)-actual['number'])))
        for label, corrupt in (('missing',records[:-1]), ('duplicate',records+[records[0]]),
                               ('reordered',list(reversed(records)))):
            try:
                aggregate_complete(grid, corrupt)
            except ValueError:
                pass
            else:
                raise AssertionError(label+' coverage corruption accepted')
        ionic = np.arange(3,dtype=float)+.7
        iq = np.zeros((2,1,1,1,3,1,3),complex)
        iq[0,0,0,0,:,0,:] = np.eye(3)*.4
        iq[1,0,0,0,:,0,:] = np.eye(3)*1.3
        joined = assemble_frozen(actual, ionic, iq)
        errors['force-sign-'+str(case)] = float(np.max(abs(joined['force_at_origin_ry_bohr']-
                                                        (ionic-expected['positive_gradient_ry_bohr']))))
        errors['q-roundtrip-'+str(case)] = float(np.max(abs(force_lattice_to_q_grid(
            joined['common_force_lattice_ry_bohr2'])-joined['common_q_grid_ry_bohr2'])))
        saved['case_'+str(case)+'_force'] = joined['force_at_origin_ry_bohr']
    errors['antiunitary-dense'] = worst_anti
    sensitivity.update(missing_D_conjugation=worst_wrong, missing_original_weights=worst_weight)

    # Exact symbolic projector scalar with independent complex channel D and rho.
    # Its derivatives check all four mixed legs, force/current signs and factors.
    u, a = sy.symbols('u a', real=True)
    b = sy.Matrix([[1+sy.I,2],[sy.Rational(1,3),1-sy.I]])
    ba = sy.Matrix([[2,sy.I],[1-sy.I,sy.Rational(2,5)]])
    baa = sy.Matrix([[sy.I,1],[2,1+sy.I]])
    p = sy.diag(sy.Rational(2,3),-sy.Rational(3,4))
    channel = sy.Matrix([[2,sy.I/3],[-sy.I/3,3]])
    rho = sy.Matrix([[sy.Rational(3,4),sy.I/5],[-sy.I/5,sy.Rational(2,5)]])
    phase = sy.diag(sy.exp(-sy.I*p[0,0]*u),sy.exp(-sy.I*p[1,1]*u))
    bu = phase*(b+a*ba+a*a*baa/2)
    scalar = sy.trace(rho*bu*channel*bu.conjugate().T)
    def numeric(x):
        return np.asarray(x.evalf(),complex)
    def at_origin(x):
        return complex(x.subs({u:0,a:0}).evalf())
    beta, first, second = map(numeric,(b,ba,baa))
    cd, rd = numeric(channel), numeric(rho)
    # Psi=I gives overlaps B^dagger, with no production jet generator reused.
    zero = np.zeros_like(beta)
    field_jets = [beta,first,zero,zero,second,zero,zero,zero,zero,zero]
    _, gradient, contact = overlap_moments(np.asarray([x.conj().T for x in field_jets]),cd,rd)
    expected_g = at_origin(sy.diff(scalar,a))
    expected_c = at_origin(sy.diff(scalar,a,2))
    errors['symbolic-field-gradient'] = abs(gradient[0]-expected_g)
    errors['symbolic-contact'] = abs(contact[0,0]-expected_c)
    pj = numeric(p)
    uj = -1j*pj@beta
    uaj = -1j*pj@first
    mixed_jets = [beta,first,zero,zero,uj,zero,zero,uaj,zero,zero,zero,zero,zero,zero,zero,zero]
    mixed = mixed_overlap_coefficient(np.asarray([x.conj().T for x in mixed_jets]),cd,rd)
    expected_m = at_origin(sy.diff(scalar,u,a))
    errors['symbolic-mixed'] = abs(mixed[0,0]-expected_m)
    ionic_second = -pj@pj@beta
    ionic_jets = [beta,uj,zero,zero,ionic_second,zero,zero,zero,zero,zero]
    _, ionic_gradient, ionic_hessian = overlap_moments(
        np.asarray([x.conj().T for x in ionic_jets]),cd,rd)
    errors['symbolic-ionic-gradient'] = abs(ionic_gradient[0]-at_origin(sy.diff(scalar,u)))
    errors['symbolic-absolute-onsite'] = abs(ionic_hessian[0,0]-at_origin(sy.diff(scalar,u,2)))
    sensitivity.update(contact_half_factor=abs(expected_c/2), mixed_omission=abs(expected_m),
                       wrong_current_sign=abs(expected_g), wrong_force_sign=2*abs(at_origin(sy.diff(scalar,u))))
    saved.update(symbolic_field_gradient=gradient, symbolic_contact=contact, symbolic_mixed=mixed)
    if max(errors.values()) > 2e-10 or min(sensitivity.values()) < 1e-5:
        raise AssertionError(dict(errors=errors,sensitivity=sensitivity))
    args.output.mkdir(parents=True,exist_ok=False)
    archive = save_arrays(args.array_root/'matched-oracle-arrays.npz',saved,1 << 20)
    from substrate_framework.material_control_reference import bound_arrays
    reread = bound_arrays(archive)
    if set(reread)!=set(saved) or any(not np.array_equal(reread[k],v) for k,v in saved.items()):
        raise AssertionError('actual saved-output roundtrip failed')
    smoke = dict(status='not-exercised', reason='no --material-slot supplied; synthetic smoke is not material evidence')
    if args.material_slot:
        receipt, arrays = receipt_arrays(binding(args.material_slot))
        if (receipt['context']['registration'] != reg_binding or
                receipt['metadata']['geometry_policy'] != 'direct-original-target' or
                any(k not in arrays or not np.isfinite(arrays[k]).all() for k in KEYS)):
            raise AssertionError('actual material output contract failed')
        residual = float(np.max(abs(arrays['current_hartree']+arrays['field_gradient']/2)))
        if residual>2e-10 or np.max(abs(arrays['number_pair']-arrays['number']))>1e-7:
            raise AssertionError('actual saved material charge/current identities failed')
        smoke = dict(status='passed-saved-output-identities',receipt=binding(args.material_slot),current_residual=residual)
    save_json(args.output/'result.json',dict(status='passed',registration=reg_binding,
        errors=errors,sensitivity=sensitivity,archive=archive,saved_output_smoke='passed',
        material_output_smoke=smoke, scope='independent equations and generic multi-case reduction; not full material execution'))


if __name__ == '__main__':
    main()
