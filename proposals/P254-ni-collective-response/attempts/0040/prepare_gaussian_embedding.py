"""Export exact Gaussian IVP data; never infer a temperature from occupations.

The native adapter consumes orthonormal coordinates, Ry energies and inverse-Ry
real times (hbar=1). Physical e-ph vertices remain Ry/bohr for an unscaled ionic
coordinate, not oscillator-normalized couplings. This is an explicit real-time
product-Gaussian boundary, not a thermal L-contour. All time history is retained.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def binding(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())

def array_binding(value):
    """Bind actual coordinate bytes; caller retains the original bank, no copy."""
    value = np.ascontiguousarray(value)
    return dict(shape=list(value.shape), dtype=value.dtype.str, order='C',
                bytes=value.nbytes, sha256=hashlib.sha256(memoryview(value).cast('B')).hexdigest())



def _hermitian(name, value):
    value = np.asarray(value, dtype=np.complex128)
    if value.ndim != 2 or value.shape[0] != value.shape[1] or not np.isfinite(value).all():
        raise ValueError(name + ': expected a finite square matrix')
    if np.linalg.norm(value - value.conj().T) > 1e-10 * (1 + np.linalg.norm(value)):
        raise ValueError(name + ': non-Hermitian; no repair is performed')
    return value


def _uniform_grid(times, order):
    times = np.asarray(times, dtype=np.float64)
    if times.ndim != 1 or len(times) <= order or times[0] != 0 or not np.isfinite(times).all():
        raise ValueError('uniform native times must start at zero and exceed startup order')
    dt = times[1]
    if dt <= 0 or not np.allclose(times, np.arange(len(times))*dt, rtol=0, atol=8*np.finfo(float).eps*max(1, times[-1])):
        raise ValueError('NESSi requires a uniform real-time mesh; no implicit interpolation')
    return times, dt


def export_input(path, rho_system, rho_bath, h_system, h_bath, coupling, times,
                 *, source_bindings, requested_times=None, order=5, extra_metadata=None):
    """Write lossless binary64 decimal input and source-bound public manifest."""
    if order not in range(1, 6):
        raise ValueError('supported NESSi integration order is 1..5')
    times, dt = _uniform_grid(times, order)
    rs, rb = _hermitian('rho_system', rho_system), _hermitian('rho_bath', rho_bath)
    ns, nb, nt = len(rs), len(rb), len(times)-1
    for name, rho in [('system', rs), ('bath', rb)]:
        occupations = np.linalg.eigvalsh(rho)
        if occupations.min() < -1e-10 or occupations.max() > 1+1e-10:
            raise ValueError(name + ': nonfermionic Gaussian state; never clipped')
    hs, hb, v = [np.asarray(a, dtype=np.complex128) for a in (h_system, h_bath, coupling)]
    if hs.shape != (nt+1, ns, ns) or hb.shape != (nt+1, nb, nb) or v.shape != (nt+1, ns, nb):
        raise ValueError('operator nodes must cover every uniform mesh node')
    if not np.isfinite(v).all():
        raise ValueError('nonfinite embedding coupling')
    for name, nodes in [('h_system', hs), ('h_bath', hb)]:
        for node in nodes:
            _hermitian(name, node)
    requested = times if requested_times is None else np.asarray(requested_times, dtype=float)
    if requested.ndim != 1 or not np.isfinite(requested).all():
        raise ValueError('requested times must be finite and one dimensional')
    indices = np.rint(requested/dt).astype(np.int64)
    if np.any(indices < 0) or np.any(indices > nt) or not np.allclose(times[indices], requested, rtol=0, atol=8*np.finfo(float).eps*max(1, times[-1])):
        raise ValueError('original requested times must be native mesh nodes; no rounding/interpolation')
    path = Path(path).resolve()
    with path.open('x') as stream:
        stream.write(f'P254_GAUSSIAN_IVP_V1 {ns} {nb} {nt} {order} {dt:.17g}\n')
        def matrix(a):
            stream.write(' '.join(f'{z.real:.17g} {z.imag:.17g}' for z in a.ravel(order='C'))+'\n')
        matrix(rs)
        matrix(rb)
        for n in range(nt+1):
            matrix(hs[n]); matrix(hb[n]); matrix(v[n])
    result = dict(schema='P254-A40-exact-Gaussian-input-v1', input=binding(path),
        source_bindings=source_bindings, units=dict(energy='Ry', time='inverse Ry', hbar=1,
        density='dimensionless, rho_ij=<c_j^dagger c_i>', embedding_sigma='Ry^2'),
        state_policy='original matrices unchanged; no log, beta, FDT, clipping or refill',
        time_representation='H[n] held on [t[n],t[n+1]); V is nodal; no memory cutoff',
        times_inverse_ry=times.tolist(), original_requested_times_inverse_ry=requested.tolist(),
        original_requested_node_indices=indices.tolist(), order=order,
        scope='real-time product Gaussian initial value, not thermal contour or material rates',
        extra_metadata=extra_metadata or {})
    path.with_suffix(path.suffix+'.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def export_prepared_boundary(path, *, projected_frame, gram, density_coefficients,
                             h_nodes, physical_g_nodes, times, system_indices,
                             bath_indices, source_bindings, requested_times=None, order=5):
    """Map original S/D/prepared basis, literal H and physical vertices, unchanged.

    projected_frame = B^dagger Psi in an orthonormal prepared basis B;
    projected_frame^dagger projected_frame = original S; rho_B=frame D frame^dagger.
    The system/bath split must be a TRUE initial product: no initial cross block
    may be dropped, even below a tolerance. For a correlated split retain those
    modes in the system and append a separately specified product bath instead.
    physical_g_nodes has shape (n_time,n_mode,n_basis,n_basis), Ry/bohr. It is
    byte-bound to the caller's retained source bank for the next material kernel,
    NOT misused as an embedding selfenergy. No old bank is copied or republished.
    """
    frame = np.asarray(projected_frame, dtype=np.complex128)
    s, d = _hermitian('original S', gram), _hermitian('original D', density_coefficients)
    if frame.ndim != 2 or frame.shape[1] != len(s) or d.shape != s.shape:
        raise ValueError('original S/D/frame dimensions disagree')
    if not np.allclose(frame.conj().T@frame, s, rtol=1e-10, atol=1e-10):
        raise ValueError('prepared basis does not reproduce original Gram matrix')
    rho = frame@d@frame.conj().T
    hs = np.asarray(h_nodes, dtype=np.complex128)
    gs = np.asarray(physical_g_nodes, dtype=np.complex128)
    n = frame.shape[0]
    if hs.shape != (len(times), n, n) or gs.ndim != 4 or gs.shape[0] != len(times) or gs.shape[2:] != (n, n) or not np.isfinite(gs).all():
        raise ValueError('literal H or physical g time-node dimensions disagree')
    a, b = np.asarray(system_indices, dtype=int), np.asarray(bath_indices, dtype=int)
    if a.ndim != 1 or b.ndim != 1 or sorted(a.tolist()+b.tolist()) != list(range(n)):
        raise ValueError('split must partition the entire prepared basis exactly once')
    if np.count_nonzero(rho[np.ix_(a,b)]) or np.count_nonzero(rho[np.ix_(b,a)]):
        raise ValueError('unsupported correlated system/bath split: original cross density is nonzero')
    path = Path(path).resolve()
    arrays = dict(projected_frame=frame, gram=s, density_coefficients=d,
                  projected_density=rho, h_nodes=hs, physical_g_nodes=gs,
                  system_indices=a, bath_indices=b, times_inverse_ry=np.asarray(times))
    actual_arrays = {name: array_binding(value) for name,value in arrays.items()}
    return export_input(path, rho[np.ix_(a,a)], rho[np.ix_(b,b)],
        hs[:,a[:,None],a], hs[:,b[:,None],b], hs[:,a[:,None],b], times,
        source_bindings=dict(original=source_bindings, actual_array_bytes=actual_arrays),
        requested_times=requested_times, order=order, extra_metadata=dict(
            physical_g_units='Ry/bohr; original unscaled ionic coordinates',
            physical_g_status='byte-bound to retained original source; not acted on by embedding solver',
            original_state='rho_B=projected_frame D projected_frame^dagger; no H-basis refill'))


def exposing_case(path, *, mutation=False):
    """Two coherent system orbitals + occupied bath, all energies in Ry."""
    h = np.array([[0.3, 0.17+0.11j], [0.17-0.11j, -0.25]], dtype=complex)
    rho = np.array([[0.5, -0.5j], [0.5j, 0.5]], dtype=complex)  # exact eigenvalues 0,1
    if mutation:
        energies, vectors = np.linalg.eigh(h)
        rho = (vectors*(1/(1+np.exp(2*energies))))@vectors.conj().T
    times = np.arange(121)*0.025
    v = np.array([[0.31+0.09j], [-0.18+0.13j]])
    return export_input(path, rho, np.ones((1,1)),
        np.broadcast_to(h, (len(times),2,2)), np.full((len(times),1,1), 0.65),
        np.broadcast_to(v, (len(times),2,1)), times,
        source_bindings=dict(kind='explicit exposing consumer, not Ni/H material input',
            producer=binding(__file__)), extra_metadata=dict(
            mutation='physical-H thermal overwrite, beta=2 inverse Ry' if mutation else 'none',
            original_occupations=[0,1,1], initial_product=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--wrong-thermal-preparation', action='store_true')
    args = parser.parse_args()
    exposing_case(args.output, mutation=args.wrong_thermal_preparation)


if __name__ == '__main__':
    main()
