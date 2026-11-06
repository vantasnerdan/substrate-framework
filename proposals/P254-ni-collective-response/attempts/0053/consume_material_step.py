"""Consume actual combined material Green functions; preserve raw complex defects."""
from pathlib import Path
import json
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0029'))
from reduced_coherence_readout import normalized_ordered


def json_value(value):
    if isinstance(value, complex):
        return dict(real=value.real, imag=value.imag)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def gamma_readout(ret, less, n):
    # Physical Ni1/x Gamma contribution: e_0/sqrt(N), original Fourier gauge.
    cross = 1j*(ret[-1,0,0,0]+less[-1,0,0,0])/n
    first, last = 1j*less[0,0,0,0]/n, 1j*less[-1,-1,0,0]/n
    return dict(connected_cross=cross,
                normalized=normalized_ordered(cross,last,first))


def consume(path):
    with np.load(path, allow_pickle=False) as bank:
        times, n = bank['times'], int(bank['n'])
        e_less = bank['electron_less']
        rho = -1j*e_less[np.arange(len(times)),np.arange(len(times))]
        raw_eigenvalues = np.linalg.eigvals(rho)
        result = dict(schema='p254-actual-fixed-kernel-material-readout-v1',
            times_inverse_ry=times,
            intrinsic=gamma_readout(bank['gamma_ret'],bank['gamma_less'],n),
            represented_control=gamma_readout(bank['control_ret'],bank['control_less'],n),
            original_bare=gamma_readout(bank['source_gamma_ret'],bank['source_gamma_less'],n),
            canonical_fermion_trace=np.trace(rho,axis1=-2,axis2=-1).sum(axis=-1),
            raw_fermion_eigenvalue_real_min=raw_eigenvalues.real.min(),
            raw_fermion_eigenvalue_real_max=raw_eigenvalues.real.max(),
            raw_fermion_eigenvalue_imaginary_max=np.abs(raw_eigenvalues.imag).max(),
            raw_fermion_hermiticity_max=np.abs(rho-rho.swapaxes(-1,-2).conj()).max(),
            fermion_initial_source_maxabs=np.abs(e_less[0,0]-bank['source_electron_less0']).max(),
            gamma_initial_source_maxabs=np.abs(bank['gamma_less'][0,0]-bank['source_gamma_less'][0,0]).max(),
            all128_canonical_fermions_consumed=rho.shape[1],
            all_Gamma_position_coordinates_consumed=bank['gamma_less'].shape[-1],
            interpretation='Fixed prescribed-reference kernels, uniform prefix only; not selfconsistent, fulltime, calibrated loss or rate. No eigenvalue/state/covariance repair.')
    return result


if __name__ == '__main__':
    result = consume(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result,default=json_value,indent=2)+'\n')
    print(json.dumps(result,default=json_value))
