"""Original A53 observables, reduced slotwise over one canonical final bank."""
import sys
import numpy as np
from history_storage import HERE
sys.path.append(str(HERE.parent/'0053'))
from consume_material_step import gamma_readout, json_value


def raw_difference(first, second):
    if first.shape != second.shape:
        raise ValueError('iterate dimensions differ')
    return max(float(np.max(np.abs(first[index]-second[index])))
               for index in np.ndindex(first.shape[:-2]))


def consume(history, previous):
    nt = len(history.times)
    trace = np.zeros(nt, dtype=complex)
    minreal, maxreal, maximag, herm, prepared = np.inf, -np.inf, 0., 0., 0.
    delta_ret, delta_less = 0., 0.
    for slot in range(128):
        ret, less = history.electron_slot(slot)
        rho = -1j*less[np.arange(nt), np.arange(nt)]
        ev = np.linalg.eigvals(rho)
        trace += np.trace(rho, axis1=-2, axis2=-1)
        minreal, maxreal = min(minreal, ev.real.min()), max(maxreal, ev.real.max())
        maximag = max(maximag, float(np.abs(ev.imag).max()))
        herm = max(herm, float(np.abs(rho-rho.swapaxes(-1, -2).conj()).max()))
        prepared = max(prepared, float(np.abs(less[0, 0]-history.aux['source_electron_less0'][slot]).max()))
        old_ret, old_less = previous.electron_slot(slot)
        delta_ret = max(delta_ret, raw_difference(ret, old_ret))
        delta_less = max(delta_less, raw_difference(less, old_less))
    ret, less = history.gamma_slot()
    old_ret, old_less = previous.gamma_slot()
    n = int(history.aux['n'])
    return dict(schema='p254-common-history-readout-v1', times_inverse_ry=history.times,
        intrinsic=gamma_readout(ret, less, n),
        represented_control=gamma_readout(history.aux['control_ret'], history.aux['control_less'], n),
        original_bare=gamma_readout(history.aux['source_gamma_ret'], history.aux['source_gamma_less'], n),
        canonical_fermion_trace=trace,
        raw_fermion_eigenvalue_real_min=float(minreal), raw_fermion_eigenvalue_real_max=float(maxreal),
        raw_fermion_eigenvalue_imaginary_max=maximag, raw_fermion_hermiticity_max=herm,
        fermion_initial_source_maxabs=prepared,
        gamma_initial_source_maxabs=float(np.abs(less[0, 0]-history.aux['source_gamma_less'][0, 0]).max()),
        all128_canonical_fermions_consumed=128, all_Gamma_position_coordinates_consumed=114,
        delta_from_previous=dict(electron_ret_maxabs=delta_ret, electron_less_maxabs=delta_less,
            gamma_ret_raw_coordinate_maxabs=raw_difference(ret, old_ret),
            gamma_less_raw_coordinate_maxabs=raw_difference(less, old_less),
            interpretation='Raw iterate differences only; mixed Gamma units are not an error norm, conservation closure or convergence criterion.'),
        interpretation='Original unprojected A53 readout formulas. No refill, Hermitization or state repair; not temporal accuracy, calibrated loss, heat, rate or closed all-q energy conservation.')
