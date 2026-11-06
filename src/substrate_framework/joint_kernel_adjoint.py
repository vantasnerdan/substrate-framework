"""Literal reverse of the saved Native spin-PBE grid action, not K=K† assumed.

PhysicalRhoFFT has grad†=-div on its exact admitted G+q support. Coulomb
is self-adjoint; all saved coefficient tensors below are explicitly conjugated
and transposed, including Native repeated mixed-spin sr/ss contractions.
"""
import numpy as np


def physical_kernel_adjoint(background, value, *, q_bohr_inv):
    fft = background.fft.at_q(q_bohr_inv)
    y = np.asarray(value, complex)
    if y.shape != (fft.nnr, 2) or not np.isfinite(y).all():
        raise ValueError('finite original physical spin grid required')
    co, grho = background.coefficients, background.grho
    rr, sr, ss, s = (co[name] for name in ('dvxc_rr', 'dvxc_sr', 'dvxc_ss', 'dvxc_s'))
    result = np.einsum('nij,ni->nj', (co['dmuxc']+rr).conj(), y)
    hbar = fft.gradient(y)
    gradient_bar = np.einsum('nij,cni->cnj', s.conj(), hbar)
    psbar = np.zeros((fft.nnr, 2, 2), complex)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                same3 = i == j == k
                a = sr[:, i, i] if same3 else sr[:, i, 1-i]
                c = sr[:, i, i] if same3 else sr[:, j, 1-j]
                psbar[:, j, k] += a.conj()*y[:, i]
                result[:, j] += c.conj()*np.sum(grho[:, :, k].conj()*hbar[:, :, i], axis=0)
                for l in range(2):
                    b = ss[:, i, i] if same3 and k == l else ss[:, i, 1-i]
                    psbar[:, j, k] += b.conj()*np.sum(grho[:, :, l].conj()*hbar[:, :, i], axis=0)
    gradient_bar += np.einsum('cni,nij->cnj', grho.conj(), psbar)
    result -= fft.divergence(gradient_bar)
    total = fft.to_g(y[:, :1]+y[:, 1:2])[:, 0]
    tpiba = 2*np.pi/background.inputs.reference.alat_bohr
    qg2 = np.sum(fft.wavevectors**2, axis=0)/tpiba**2
    selected = qg2 > 1e-8  # Existing Native admission, not a material tolerance.
    potential = np.zeros_like(total)
    potential[selected] = 8*np.pi*total[selected]/(tpiba**2*qg2[selected])
    return result+fft.from_g(potential[:, None])
