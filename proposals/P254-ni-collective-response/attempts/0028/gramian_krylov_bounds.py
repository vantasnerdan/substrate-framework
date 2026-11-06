"""Source-aware full-row Duhamel estimate using a small residual Gramian.

For B=-i*T (opposite sign backwards), A=B-alpha*I and c=beta*e0,
W(t)=integral exp(A†s) R†R exp(As) ds. Cauchy-Schwarz gives
integral exp(alpha*(t-s))*||R exp(Bs)c|| ds
    <= exp(alpha*t)*sqrt(t*c†W(t)c).
The Van Loan block exp([[-A†,R†R],[0,A]]*t) has top-right Z and bottom-right V;
W=V†Z. This identity needs neither orthogonal Q nor Hermitian T. Thin QR is
used only with its measured reconstruction defect and computed ||Q|| supplement.

Only the small block exponential uses selective extra precision. The full-PW
residual and H actions remain the original floating data. This is a numerical
evaluation of an exact-arithmetic inequality, not an interval certificate.
"""
import math

import mpmath as mp
import numpy as np


def _positive_float(value):
    result = float(value)
    return float(np.nextafter(0., np.inf)) if value > 0 and result == 0 else result


def gramian_envelope(projected, residual, coefficient, initial_error, time, *,
                     ambient_growth, seed_action_defect=0., target=None, digits=50):
    tau = abs(float(time))
    if (not np.isfinite([tau, ambient_growth, initial_error, seed_action_defect]).all()
            or min(ambient_growth, initial_error, seed_action_defect) < 0):
        raise ValueError('finite nonnegative growth and defect data required')
    if tau == 0:
        return dict(absolute_error_envelope_estimate=float(initial_error), method='initial-defect')
    b = (-1j if time >= 0 else 1j)*np.asarray(projected)
    growth = max(0., float(np.linalg.eigvalsh((b+b.conj().T)/2)[-1]))
    initial = math.exp(ambient_growth*tau)*initial_error
    weight = tau*math.exp(max(ambient_growth, growth)*tau)*abs(coefficient)
    cheap = initial+weight*(float(np.linalg.norm(residual))+seed_action_defect)
    if target is not None and cheap <= target:
        return dict(absolute_error_envelope_estimate=cheap,
                    cheap_envelope_estimate=cheap, method='full-residual-norm')
    q, r = np.linalg.qr(residual, mode='reduced')
    compression = float(np.linalg.norm(residual-q@r))
    r *= max(1., float(np.linalg.norm(q, 2)))
    del q
    size = b.shape[0]
    with mp.workdps(digits):
        convert = lambda z: mp.mpc(float(np.real(z)), float(np.imag(z)))
        small = mp.matrix([[convert(z) for z in row] for row in b])
        small -= mp.mpf(ambient_growth)*mp.eye(size)
        thin = mp.matrix([[convert(z) for z in row] for row in r])
        gram = thin.H*thin
        block = mp.matrix(2*size)
        for i in range(size):
            for j in range(size):
                block[i,j] = -mp.conj(small[j,i])
                block[i,j+size] = gram[i,j]
                block[i+size,j+size] = small[i,j]
        flow = mp.expm(mp.mpf(tau)*block)
        energy = abs(convert(coefficient))**2*sum(
            mp.conj(flow[j+size,size])*flow[j,size] for j in range(size))
        if mp.re(energy) < 0:
            raise ArithmeticError('negative computed residual-energy Gramian; precision unresolved')
        gramian = _positive_float(mp.exp(mp.mpf(ambient_growth)*tau)*mp.sqrt(tau*abs(energy)))
        imaginary = str(mp.im(energy))
    supplement = weight*(compression+seed_action_defect)
    candidate = initial+gramian+supplement
    return dict(absolute_error_envelope_estimate=min(cheap, candidate),
        cheap_envelope_estimate=cheap, gramian_envelope_estimate=candidate,
        initial_error_term=initial, residual_energy_term=gramian,
        reconstruction_supplement=supplement, residual_compression_defect=compression,
        residual_energy_imaginary=imaginary, small_block_digits=digits,
        method='weighted-source-residual-Gramian',
        qualification='floating full-row data and selective-precision small exponential; not interval certified')
