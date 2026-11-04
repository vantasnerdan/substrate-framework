"""Exact four-orbital Fock oracle for the NEW dual/observable boundary only."""
import json
import numpy as np
import sympy as s
from homogeneous_electronic_output import homogeneous_electronic_operators
from export_frozen_bath_comparator import car_moments, projected_factors


def dense(value):
    return np.asarray(value.tolist(),dtype=np.complex128)


def expanded(value):
    return value.applyfunc(s.expand)


def annihilator(index, orbitals):
    out = s.zeros(2**orbitals)
    for state in range(2**orbitals):
        if state & (1 << index):
            out[state ^ (1 << index),state] = (-1)**((state & ((1 << index)-1)).bit_count())
    return out


def main():
    I = s.I
    r = s.Rational
    identity = s.eye(16)
    annihilators = [annihilator(index,4) for index in range(4)]
    for a in range(4):
        for b in range(4):
            assert expanded(annihilators[a]*annihilators[b].H+
                annihilators[b].H*annihilators[a]) == (identity if a == b else s.zeros(16))
    u = s.Matrix([[r(3,5),4*I/5,0,0],[4*I/5,r(3,5),0,0],
                  [0,0,r(5,13),r(12,13)],[0,0,-r(12,13),r(5,13)]])
    assert expanded(u.H*u) == s.eye(4)
    transform = s.Matrix([[1,I/5,0],[0,2,r(1,7)],[r(1,4),0,r(3,2)]])
    inverse = transform.inv().applyfunc(s.expand_complex)
    populations = [r(1,5),r(2,5),r(3,5),s.Integer(0)]
    psi = expanded(u[:,:3]*transform)
    density = expanded(inverse*s.diag(*populations[:3])*inverse.H)
    selector = s.Matrix([[r(3,5),0],[r(4,5),0],[0,1]])
    raw_coordinates = expanded(inverse*selector)
    basis = expanded(psi*raw_coordinates)
    assert expanded(basis.H*basis) == s.eye(2)
    fock_density = identity
    for index,population in enumerate(populations):
        c = sum((s.conjugate(u[j,index])*annihilators[j] for j in range(4)),s.zeros(16))
        number = expanded(c.H*c)
        fock_density = expanded(fock_density*((1-population)*identity+(2*population-1)*number))
    assert s.expand(s.trace(fock_density)) == 1
    z = s.Matrix([[r(2,7)+I/9,r(1,3)-2*I/5],[r(4,11)+I/6,-r(1,5)+I/8]])
    other_z = s.Matrix([[2+I,r(3,2)-I/7],[-r(2,3)+I/4,3-I/5]])
    nx,ns,m = 3,2,2
    dual = np.zeros((nx+ns*m*m,1),complex)
    dual[:nx,0] = [17+3j,-5+4j,11-2j]
    dual[nx:,0] = np.stack((dense(other_z),dense(z))).reshape(-1)
    operator = homogeneous_electronic_operators(dual,dense(basis),mechanical_size=nx,
                                               source_count=ns,source_slot=1)[0]
    cactive = [sum((s.conjugate(basis[j,a])*annihilators[j] for j in range(4)),s.zeros(16))
               for a in range(2)]
    # Independent definition: coefficient z_ab multiplies c_b^dagger c_a.
    fock_homogeneous = expanded(sum((z[a,b]*cactive[b].H*cactive[a]
                                    for a in range(2) for b in range(2)),s.zeros(16)))
    extra = s.Matrix([[r(1,7),I/8,r(1,4),0],[-I/8,-r(1,9),0,r(1,6)],
                      [r(1,4),0,r(1,5),I/3],[0,r(1,6),-I/3,r(2,7)]])
    fock_extra = expanded(sum((extra[a,b]*annihilators[a].H*annihilators[b]
                               for a in range(4) for b in range(4)),s.zeros(16)))
    fock = dict(homogeneous=fock_homogeneous,omitted=fock_extra,
                combined=expanded(fock_homogeneous+fock_extra))
    truth = {}
    for a,left in fock.items():
        for b,right in fock.items():
            truth[a,b] = complex(s.expand(s.trace(fock_density*left*right)
                -s.trace(fock_density*left)*s.trace(fock_density*right)))
    frame,gram,d = dense(psi),dense(expanded(psi.H*psi)),dense(density)
    def evaluate(candidate):
        matrices = dict(homogeneous=candidate,omitted=dense(extra),combined=candidate+dense(extra))
        adjoints = {name:value.conj().T for name,value in matrices.items()}
        result = car_moments(frame,gram,d,projected_factors(matrices,adjoints,frame,
            'exact-oracle-four-orbital','original-nonorthogonal-state','bohr'),
            'original-raw-frame','exact-oracle-four-orbital','original-nonorthogonal-state')
        return {(row['left'],row['right']):row['value'] for row in result['covariance']}
    actual = evaluate(operator)
    error = max(abs(actual[key]-value) for key,value in truth.items())
    assert error < 3e-12,(error,actual,truth)
    mutations = {
        'conjugate_transpose_dual':dense(basis*z.H*basis.H),
        'wrong_saved_source_slot':homogeneous_electronic_operators(dual,dense(basis),
            mechanical_size=nx,source_count=ns,source_slot=0)[0],
        'extra_fourier_factor_on_pulledback_output':operator/8,
        'raw_columns_instead_of_saved_active_coordinates':dense(psi[:,:2]*z.T*psi[:,:2].H)}
    sensitivity = {name:max(abs(value[key]-truth[key]) for key in truth)
                   for name,candidate in mutations.items() for value in [evaluate(candidate)]}
    sensitivity['drop_both_original_output_cross_terms'] = abs(truth['combined','combined']
        -truth['homogeneous','homogeneous']-truth['omitted','omitted'])
    assert min(sensitivity.values()) > 1e-6,sensitivity
    print(json.dumps(dict(schema='P254-homogeneous-electronic-output-owning-oracle-v1',
        exact_fock_ordered_pairs=len(truth),max_complex_covariance_error=error,
        exposed_mutations=sensitivity,material_inputs_read=False,
        scope='New complex-dual/source-slot/metric embedding and original-output cross boundary only;existing filter/bath/state proofs not rerun.'),indent=2))


if __name__ == '__main__':
    main()
