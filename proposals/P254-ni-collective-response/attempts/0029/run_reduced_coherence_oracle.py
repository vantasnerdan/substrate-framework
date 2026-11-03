"""Independent finite Fock/analytic coherence oracle; no material inputs."""
from pathlib import Path
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.linalg import expm
import sympy as sp
from reduced_coherence_readout import normalized_ordered


def oracle():
    va, vb, a, b = sp.symbols('va vb a b', positive=True)
    c = sp.symbols('c', complex=True)
    assert sp.simplify(a*b*c/(sp.sqrt(a*a*va)*sp.sqrt(b*b*vb))-c/(sp.sqrt(va)*sp.sqrt(vb))) == 0
    annihilation = np.diag(np.sqrt(np.arange(1, 3)), 1).astype(complex)
    x = (annihilation+annihilation.conj().T)/np.sqrt(2)
    p = (annihilation-annihilation.conj().T)/(1j*np.sqrt(2))
    number = annihilation.conj().T@annihilation
    vacuum = np.array([1., 0., 0.], complex)
    expectation = lambda op: vacuum.conj()@op@vacuum
    time = np.pi/2
    u = expm(-1j*number*time)
    evolved = u.conj().T@x@u
    cross, vl, vr = expectation(evolved@x), expectation(evolved@evolved), expectation(x@x)
    actual = normalized_ordered([cross], vl, vr)['complex_degree'][0]
    assert abs(actual+1j) < 2e-14
    wrong_real = normalized_ordered([np.real(cross)], vl, vr)['complex_degree'][0]
    assert abs(wrong_real-actual) > .9
    squeezed = 2*x+3*p
    squeezed_cross = expectation(squeezed@x)
    squeezed_left = expectation(squeezed@squeezed)
    correct = normalized_ordered([squeezed_cross], squeezed_left, vr)['complex_degree'][0]
    assert abs(correct-(2-3j)/np.sqrt(13)) < 2e-14
    assert abs(squeezed_cross/vr-correct) > 1
    scale_defects = []
    for scale in (1e-250, 1e250):
        result = normalized_ordered([scale*cross], scale*vl, scale*vr)
        assert result['defined']
        scale_defects.append(float(abs(result['complex_degree'][0]-actual)))
    assert max(scale_defects) < 2e-14
    identity = np.eye(3)
    h = np.kron(number, identity)+2*np.kron(identity, number)
    observable = (np.kron(x, identity)+np.kron(identity, x))/np.sqrt(2)
    state = np.kron(vacuum, vacuum)
    evolution = expm(-1j*h*np.pi)
    evolved_two = evolution.conj().T@observable@evolution
    avg = lambda op: state.conj()@op@state
    two = normalized_ordered([avg(evolved_two@observable)], avg(evolved_two@evolved_two), avg(observable@observable))
    assert abs(two['complex_degree'][0]) < 2e-14
    zero = normalized_ordered([0j], 0., vr)
    assert not zero['defined']
    return dict(single_mode_ordered_degree=[float(actual.real), float(actual.imag)],
                squeezing_degree=[float(correct.real), float(correct.imag)],
                unit_scale_defects=scale_defects,
                two_mode_unitary_dephasing_magnitude=float(abs(two['complex_degree'][0])),
                real_only_mutation_defect=float(abs(wrong_real-actual)),
                fixed_initial_variance_mutation_defect=float(abs(squeezed_cross/vr-correct)),
                scope='Finite operator normalization only; dephasing is not irreversible loss or material lifetime')


if __name__ == '__main__':
    print(json.dumps(oracle(), indent=2, sort_keys=True))
