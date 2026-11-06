"""Bare Ry Coulomb factors with original directed-edge and Nyquist closure maps.

For normalized cell orbitals D has units Omega^-1/2. The physical Fourier
multiplier is 8*pi/Q^2, without another volume. Canonical supercell vertices
therefore contain exactly (1/N)*D*D. H_int=1/2 V_ij;kl c†i c†j c_l c_k,
A_ij;kl=V_ij;kl-V_ij;lk and H_int=1/4 A_ij;kl c†i c†j c_l c_k.
No saved full PBE K, instantaneous mean or screened interaction is consumed.
"""
import numpy as np


def coulomb_multiplier(wavevectors, tpiba):
    q2 = np.sum(np.asarray(wavevectors)**2, axis=-1)
    admitted = q2/tpiba**2 > 1.e-8
    value = np.zeros(q2.shape)
    value[admitted] = 8*np.pi/q2[admitted]
    return value, admitted


class BareCoulombFactors:
    def __init__(self, inputs):
        self.c = inputs
        self.weights, self.wavevectors, self.admitted = [], [], []
        for q, gs in enumerate(inputs.feature_g):
            frac = inputs.q_indices[q]/4.
            frac -= np.floor(frac+.5)
            waves = (inputs.reciprocal@(gs+frac).T).T
            weight, admitted = coulomb_multiplier(waves, inputs.tpiba)
            self.weights.append(weight)
            self.wavevectors.append(waves)
            self.admitted.append(admitted)
        for q in range(inputs.nq):
            qm = int(inputs.minus[q])
            partner = inputs.feature_partner[q][:len(inputs.feature_g[q])]
            if not np.allclose(self.wavevectors[q], -self.wavevectors[qm][partner], rtol=0, atol=1e-13):
                raise ValueError('opposite physical Fourier feature differs')

    def transfer(self, target, source):
        c = self.c
        delta = (c.k_indices[c.slot_k[target]]-c.k_indices[c.slot_k[source]])%4
        return int(np.ravel_multi_index(tuple(delta), (4,4,4)))

    def direct(self, i, a, b, c):
        """Return X[r,i_band,b_band],Y[r,a_band,c_band] for V_ia;bc.

        Opposite physical multiplication legs come from their original archives,
        not an assumed adjoint. Explicit spin deltas can make a term absent.
        The first factor includes v/N, so two vertices give precisely 1/N^2.
        """
        data = self.c
        if data.slot_spin[i] != data.slot_spin[b] or data.slot_spin[a] != data.slot_spin[c]:
            return None
        q = self.transfer(i,b)
        qm = int(data.minus[q])
        if int(data.target[q,b]) != i or int(data.target[qm,c]) != a:
            return None
        ng = len(data.feature_g[q])
        partner = data.feature_partner[q][:ng]
        x = data.density[q][b].conj()*(self.weights[q]/data.n)[:,None,None]
        y = data.density[qm][c,partner].conj()
        return x,y

    def quartets(self, external):
        """All spin/momentum conserving (a,b,c); no internal band truncation.

        Enumerate every a and every b momentum. Allowed b spin is external or a;
        the corresponding c spin is determined. Equal spins are enumerated once.
        """
        data = self.c
        si = external//data.n
        ki = data.k_indices[data.slot_k[external]]
        for a in range(data.ns):
            sa = a//data.n
            ka = data.k_indices[data.slot_k[a]]
            for sb in sorted({si,sa}):
                sc = sa if sb == si else si
                for kb in range(data.n):
                    kc = (ki+ka-data.k_indices[kb])%4
                    c = sc*data.n+int(np.ravel_multi_index(tuple(kc), (4,4,4)))
                    yield a, sb*data.n+kb, c

    def term_factors(self, i, a, b, c):
        """Both independent vertices: (ia;bc) at t and (bc;ia) at s."""
        left = (self.direct(i,a,b,c), self.direct(i,a,c,b))
        right = (self.direct(b,c,i,a), self.direct(b,c,a,i))
        return left,right
