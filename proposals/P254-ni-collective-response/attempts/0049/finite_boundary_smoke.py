"""Actual new C++ supplied API vs independent SciPy full phase-space evolution."""
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad_vec
from finite_gaussian_reference import fermion_case, boson_case
from supplied_vie2 import solve_supplied, load_solution
from command_receipts import run_recorded, binding


def error(actual, expected):
    return float(np.max(np.abs(actual-expected)))


def run_case(factory, statistics, nt, binary, root, label, beta=1.):
    times = np.linspace(0,.6,nt+1)
    case = factory(times)
    path = solve_supplied(times,*case['source'],*case['kernel'], statistics=statistics,
                          order=3,binary=binary,workdir=root/label,beta=beta)
    result = load_solution(path)
    report = {'solution':binding(path),
              'ret_error':error(result['ret'],case['exact'][0]),
              'less_error':error(result['less'],case['exact'][1]),
              'source_ret_diagonal_copy_error':max(error(result['ret'][n,n],case['source'][0][n,n]) for n in range(nt+1)),
              'discarded_preparation_error':error(result['less'],case['discarded_preparation'][1]),
              'wrong_lesser_sign_error':error(result['less'],-case['exact'][1]),
              'f_fcc_difference':error(result['f_ret'],result['fcc_ret'])}
    # Independently adaptive quadrature of the actual continuum retarded products.
    # It checks each noncommuting order separately, not just a norm discrepancy.
    final, earlier = times[-1], times[nt//2]
    def integrand(s):
        probe = factory(np.array([earlier,s,final]))
        q, kernel = probe['source'][0], probe['kernel'][0]
        return np.stack([-q[2,1]@kernel[1,0],-kernel[2,1]@q[1,0]])
    convolution, quadrature_error = quad_vec(integrand,earlier,final,epsabs=1e-12,epsrel=1e-12)
    report['f_ret_quadrature_error'] = error(result['f_ret'][-1,nt//2],convolution[0])
    report['fcc_ret_quadrature_error'] = error(result['fcc_ret'][-1,nt//2],convolution[1])
    report['wrong_fcc_order_oracle_error'] = error(result['f_ret'][-1,nt//2],convolution[1])
    report['adaptive_quadrature_estimate'] = float(quadrature_error)
    if statistics == 1:
        report['wrong_mass_reference_error'] = error(result['less'],case['wrong_mass'][1])
        report['wrong_oscillator_scaling_error'] = error(result['less'],case['exact'][1]*np.sqrt(case['masses'][:2,None]*case['masses'][None,:2]))
        report['uncertainty_min_eigenvalue'] = case['uncertainty_min_eigenvalue']
        report['coordinate_covariance_stiffness_commutator'] = case['coordinate_covariance_stiffness_commutator']
        report['initial_position_covariance'] = (1j*result['less'][0,0]).real.tolist()
    else:
        report['rho_eigenvalues'] = case['rho_eigenvalues'].tolist()
        report['initial_coherence_abs'] = float(abs(case['rho'][0,1]))
    return path,result,case,report


def main(binary, root):
    root = Path(root)
    reports = {}
    saved = {}
    for label,factory,statistics in [('fermion',fermion_case,-1),('boson_position',boson_case,1)]:
        for nt in (20,40):
            key = f'{label}-{nt}'
            saved[key] = run_case(factory,statistics,nt,binary,root,key)
            reports[key] = saved[key][3]
            (root/'finite-boundary-report.json').write_text(json.dumps(reports,indent=2)+'\n')
        coarse,fine = reports[f'{label}-20'],reports[f'{label}-40']
        coarse['ret_to_fine_error_ratio'] = coarse['ret_error']/max(fine['ret_error'],1e-30)
        coarse['less_to_fine_error_ratio'] = coarse['less_error']/max(fine['less_error'],1e-30)
    _,result,case,_ = saved['boson_position-20']
    times = np.linspace(0,.6,21)
    wrong_sign = solve_supplied(times,case['source'][0],-case['source'][1],case['kernel'][0],-case['kernel'][1],
                               statistics=1,order=3,binary=binary,workdir=root/'boson-wrong-lesser-sign')
    reports['actual_wrong_lesser_sign_error'] = error(load_solution(wrong_sign)['less'],case['exact'][1])
    # Zero mixed storage makes the contour statistics sign inert numerically.
    alternate = solve_supplied(times,*case['source'],*case['kernel'],statistics=-1,order=3,
                               binary=binary,workdir=root/'boson-inert-statistics-beta',beta=3.7)
    reports['zero_mixed_wrong_statistics_beta_difference'] = error(load_solution(alternate)['less'],result['less'])
    bad_ret = case['source'][0].copy()
    for n in range(len(times)):
        bad_ret[n,n] = -1j*np.eye(2)
    bad_diagonal = solve_supplied(times,bad_ret,case['source'][1],*case['kernel'],statistics=1,order=3,
                                 binary=binary,workdir=root/'boson-wrong-canonical-diagonal')
    # Deliberately incompatible canonical replacement can also violate the
    # upstream solution's contour-Hermiticity; the public reader must reject it.
    try:
        load_solution(bad_diagonal)
    except ValueError as exc:
        reports['wrong_source_diagonal_consumer_rejection'] = str(exc)
    with np.load(bad_diagonal,allow_pickle=False) as raw:
        reports['actual_wrong_source_diagonal_error'] = error(raw['ret'],case['exact'][0])
    for label, altered_times, altered_less in [
            ('nonuniform',np.array([0.,1.,10.,30.]),case['source'][1][:4,:4]),
            ('independent-lesser',times,case['source'][1].copy())]:
        if label == 'independent-lesser':
            altered_less[1,0,0,1] += .1
        count = len(altered_times)
        try:
            solve_supplied(altered_times,case['source'][0][:count,:count],altered_less,
                           case['kernel'][0][:count,:count],case['kernel'][1][:count,:count],
                           statistics=1,order=3,binary=binary,workdir=root/f'reject-{label}')
        except ValueError as exc:
            reports[f'{label}_rejection'] = str(exc)
        else:
            raise AssertionError(f'{label} was silently accepted')
    for label in ('fermion','boson_position'):
        solution = saved[f'{label}-40'][0]
        run_recorded([sys.executable,str(Path(__file__).with_name('consume_solution.py')),str(solution),
                      str(root/f'{label}-consumer.json')],root,f'{label}-saved-consumer')
    # Persist observed numbers BEFORE acceptance, including any failing run.
    (root/'finite-boundary-report.json').write_text(json.dumps(reports,indent=2)+'\n')
    for label in ('fermion','boson_position'):
        fine = reports[f'{label}-40']
        assert fine['ret_error'] < 2e-6 and fine['less_error'] < 2e-6, (label,fine)
        assert fine['source_ret_diagonal_copy_error'] == 0
        assert fine['f_ret_quadrature_error'] < 2e-6 and fine['fcc_ret_quadrature_error'] < 2e-6
        assert fine['wrong_fcc_order_oracle_error'] > 1e-6
        assert fine['discarded_preparation_error'] > 1e-3
        assert fine['wrong_lesser_sign_error'] > .1
        assert reports[f'{label}-20']['ret_to_fine_error_ratio'] > 2
        assert reports[f'{label}-20']['less_to_fine_error_ratio'] > 2
    assert reports['actual_wrong_lesser_sign_error'] > .1
    assert reports['actual_wrong_source_diagonal_error'] > .9
    assert reports['zero_mixed_wrong_statistics_beta_difference'] < 1e-12
    assert reports['boson_position-40']['wrong_mass_reference_error'] > 1e-3
    assert reports['boson_position-40']['wrong_oscillator_scaling_error'] > 1e-3
    print(json.dumps(reports))

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])
