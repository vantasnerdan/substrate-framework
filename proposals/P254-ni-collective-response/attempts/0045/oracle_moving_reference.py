"""New finite oracle: noncommuting exact products, sign/order/density sensitivity."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from moving_reference_inputs import atomic_json, file_binding
from moving_reference_propagation import integrate_u, checkpoint_u


def run(output):
    output=Path(output);output.mkdir(exist_ok=False)
    hx=np.array([[.3,.7+.2j],[.7-.2j,-.1]],complex)
    hz=np.array([[.9,-.1j],[.1j,-.6]],complex)
    eye=np.eye(2,dtype=complex)
    first,w1=integrate_u(lambda t:hx,0,1,eye)
    restored,binding=checkpoint_u(output,1,first,[],w1)
    actual,w2=integrate_u(lambda t:hz,1,2,restored)
    exact=expm(-1j*hz)@expm(-1j*hx)
    wrong_order=expm(-1j*hx)@expm(-1j*hz)
    wrong_sign=expm(1j*hz)@expm(1j*hx)
    # All these U preserve the same full-rank maximally mixed density: rho
    # cannot identify transition phases or time ordering. I is one admissible
    # rho-only inferred U, and is detectably wrong as a two-time propagator.
    rho=.5*eye
    errors=dict(correct=float(np.linalg.norm(actual-exact)),
        wrongorder=float(np.linalg.norm(actual-wrong_order)),
        wrongsign=float(np.linalg.norm(actual-wrong_sign)),
        rho_inferred_identity=float(np.linalg.norm(actual-eye)),
        density_indistinguishability=float(np.linalg.norm(actual@rho@actual.conj().T-rho)),
        restart_byte_difference=float(np.linalg.norm(restored-first)))
    with (output/'exact-products.npz').open('xb') as stream:
        np.savez(stream,hx=hx,hz=hz,actual=actual,exact=exact,
            wrongorder=wrong_order,wrongsign=wrong_sign,rho=rho,rho_inferred_u=eye)
    passed=(errors['correct']<2e-7 and errors['wrongorder']>.1 and
            errors['wrongsign']>.1 and errors['rho_inferred_identity']>.1 and
            errors['density_indistinguishability']<2e-7 and errors['restart_byte_difference']==0)
    record=dict(schema='p254-a45-noncommuting-oracle-v1',passed=passed,
        errors=errors,work=[w1,w2],actual_reloaded_boundary=binding,
        arrays=file_binding(output/'exact-products.npz'),
        independence='SciPy matrix exponential exact-product reference versus adaptive DOP853 ODE; no shared time-order integration',
        scope='finite piecewise constant2level algorithm/restart sensitivity, not material accuracy certification',
        recovery='persisted U at t=1 actually reloaded; adaptive solver history lost')
    atomic_json(output/'result.json',record)
    if not passed:
        raise AssertionError(errors)
    return file_binding(output/'result.json')

if __name__=='__main__':
    import sys
    print(json.dumps(run(sys.argv[1])),flush=True)
