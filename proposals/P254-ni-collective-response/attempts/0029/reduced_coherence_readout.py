"""Selected ordered two-point coherence, not a fitted lifetime or entanglement.

For Hermitian A=u_r(t), B=u_0(s), report C_AB/sqrt(C_AA C_BB).
The homogeneous unexcited baseline makes its anchor equal-time variance the
same at every cell; every spatial q was still computed, not symmetry-dropped.
Variances and their complex defects remain unchanged. No decay cutoff or
continuous-time conclusion is invented from a few saved samples.
"""
import numpy as np


def normalized_ordered(cross, variance_left, variance_right):
    cross = np.asarray(cross, dtype=complex)
    left, right = complex(variance_left), complex(variance_right)
    if not np.isfinite(cross).all() or not np.isfinite([left, right]).all():
        raise ValueError('finite actual ordered covariance and variances required')
    result = dict(variance_left=left, variance_right=right,
                  variance_real_parts_positive=left.real > 0 and right.real > 0,
                  variance_imaginary_parts=[left.imag, right.imag])
    if left == 0 or right == 0:
        return dict(result, defined=False, reason='zero actual variance; normalized coherence is undefined')
    # Separate roots avoid squaring the variance scale into overflow/underflow.
    denominator = np.sqrt(left)*np.sqrt(right)
    if denominator == 0 or not np.isfinite(denominator):
        return dict(result, defined=False, reason='normalization outside representable numeric range')
    value = cross/denominator
    if not np.isfinite(value).all():
        raise ValueError('normalized ordered correlation is outside representable range')
    magnitude = np.abs(value)
    return dict(result, defined=True, complex_degree=value,
                magnitude=magnitude, phase_radians=np.angle(value),
                cauchy_magnitude_surplus=magnitude-1,
                definition='ordered C(r,t;s)/[sqrt(C(0,t;t))*sqrt(C(0,s;s))]; raw complex variances')


def coherence_curve(observations):
    available = {row['time_inverse_ry']: row['connected_correlation']
                 for row in observations if row['connected_correlation']['available']}
    rows = []
    defined_times = []
    for observation in observations:
        time = observation['time_inverse_ry']
        correlation = observation['connected_correlation']
        row = dict(time_inverse_ry=time, connected_available=correlation['available'])
        if not correlation['available']:
            rows.append(row)
            continue
        reference = correlation['reference_time_inverse_ry']
        row.update(reference_time_inverse_ry=reference, noise_scope=correlation['noise_scope'],
                   equal_time_ordered_by_cell_bohr2=correlation['equal_time']['ordered_ab'],
                   two_time_ordered_by_cell_bohr2=correlation['two_time']['ordered_ab'],
                   two_time_symmetrized_by_cell_bohr2=correlation['two_time']['symmetrized'],
                   two_time_commutator_by_cell_bohr2=correlation['two_time']['commutator'])
        if reference not in available:
            row['ordered_coherence'] = dict(defined=False,
                reason='no saved equal-time variance at the exact reference time')
        else:
            row['ordered_coherence'] = normalized_ordered(
                correlation['two_time']['ordered_ab'], correlation['equal_time']['ordered_ab'][0],
                available[reference]['equal_time']['ordered_ab'][0])
            if row['ordered_coherence']['defined']:
                defined_times.append(time)
        rows.append(row)
    return dict(status='computed' if defined_times else 'normalization_unavailable',
                observable='selected Hermitian Ni-displacement ordered degree of coherence',
                actual_normalized_sample_times_inverse_ry=defined_times,
                sampled_time_span_inverse_ry=None if not defined_times else [min(defined_times), max(defined_times)],
                samples=rows, continuous_time_window_certified=False,
                physical_lifetime_determined=False,
                scope='Moving unexcited active-Gaussian model; not finite-kick covariance, passive noise or entanglement',
                interpretation='Magnitude and phase distinguish ordered coherence from an oscillating real covariance. A sampled span is not a decay time or a calibrated material coherence window.')
