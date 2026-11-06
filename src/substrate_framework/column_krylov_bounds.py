"""Full-row Duhamel envelopes for independent column Arnoldi approximations.

Exact-arithmetic inequalities evaluated in ordinary floating point, NOT interval
certificates. No Hermiticity or orthogonality is assumed in the defect identity.
"""
from math import factorial

import numpy as np
from scipy.linalg import expm


def exponential(value):
    return float(np.exp(value)) if value < 709.0 else float('inf')


def column_envelope(projected, residual, coefficient, initial_error, time, *,
                    ambient_growth, seed_action_defect=0., cells=4, degree=4):
    """Integrate the actual full-row defect, not a last Ritz coefficient.

    E=HQ-shift*Q-Q*T includes all finite-precision decomposition and numerical
    breakdown directions. A thin QR compresses E for repeated small products;
    the measured reconstruction defect is added using a separate norm envelope.
    The same supplement includes the seed scaling/reconstruction defect. Neither
    QR/expm/norm roundoff nor errors in the supplied H action are interval bounded.
    """
    tau = abs(float(time))
    if tau == 0:
        return dict(absolute_error_envelope_estimate=float(initial_error),
                    cheap_envelope_estimate=float(initial_error),
                    cell_taylor_envelope_estimate=float(initial_error),
                    initial_error_term=float(initial_error), projected_forward_growth=0.,
                    projected_two_sided_growth=0., residual_compression_defect=0.,
                    seed_action_defect=float(seed_action_defect))
    b = (-1j if time >= 0 else 1j)*projected
    eigen = np.linalg.eigvalsh((b+b.conj().T)/2)
    growth = max(0., float(eigen[-1]))
    two_sided = float(np.max(np.abs(eigen)))
    q, r = np.linalg.qr(residual, mode='reduced')
    compression = float(np.linalg.norm(residual-q @ r))
    # Do not assume that computed QR columns are exactly isometric.
    r *= max(1., float(np.linalg.norm(q, 2)))
    del q
    c = np.zeros(projected.shape[0], dtype=complex)
    c[0] = coefficient
    initial = exponential(ambient_growth*tau)*initial_error
    weight = tau*exponential(max(ambient_growth, growth)*tau)*abs(coefficient)
    supplement = weight*(compression+seed_action_defect)
    cheap = initial+weight*(float(np.linalg.norm(residual))+seed_action_defect)
    remainder = float(np.linalg.norm(r @ np.linalg.matrix_power(b, degree)))
    half = tau/(2*cells)
    local = initial+supplement
    for cell in range(cells):
        left, middle = 2*cell*half, (2*cell+1)*half
        y = expm(middle*b) @ c
        power_y = y
        integral = 0.
        for order in range(degree):
            integral += 2*half**(order+1)/factorial(order+1)*float(np.linalg.norm(r @ power_y))
            power_y = b @ power_y
        integral += (2*half**(degree+1)/factorial(degree+1)*remainder
                     *float(np.linalg.norm(y))*exponential(two_sided*half))
        local += exponential(ambient_growth*(tau-left))*integral
    # Overflow/NaN cannot produce a false convergence verdict.
    candidates = [value if np.isfinite(value) else float('inf') for value in (cheap, local)]
    return dict(absolute_error_envelope_estimate=min(candidates),
                cheap_envelope_estimate=float(cheap), cell_taylor_envelope_estimate=float(local),
                initial_error_term=float(initial), projected_forward_growth=growth,
                projected_two_sided_growth=two_sided, residual_compression_defect=compression,
                seed_action_defect=float(seed_action_defect))
