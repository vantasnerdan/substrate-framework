"""Complex trace-bilinear Gamma gradients, positive dH/dR in Ry/bohr.

This is an ordinary-NC, no-NLCC reduction, not a Native real-force driver.
For X=L R^dagger both differentiated projector sides are retained separately;
no real part, spin factor, force minus sign or symmetry processing is applied.
The density is physical bohr^-3. Its single inverse FFT below implements the
BILINEAR Fourier transpose, not the conjugate density/potential pairing.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg.blas import zgemm


class GammaGradientReduction:
    """Share species projectors and block overlaps over every ionic coordinate.

    The owning GammaDyadCoupling validates the common PW/grid/geometry identity.
    No global block-diagonal D, derivative projector array or PW-square matrix
    is formed. Overlaps live only for the current factor block and species.
    """

    def __init__(self, factors, ionic):
        self.factors = factors
        self.ionic = ionic
        self.nat = len(ionic.species)
        self.omega = factors.fft.omega
        self.local_values = {}
        moduli = np.linalg.norm(ionic.g, axis=1)
        for symbol in dict.fromkeys(ionic.species):
            if symbol not in ionic.vloc:
                ionic.vloc[symbol] = ionic.tables[symbol].vloc(moduli)
            self.local_values[symbol] = ionic.vloc[symbol]

    def local(self, density):
        """All atoms from ONE transform: Omega sum_G n(-G) dV(G)/dR."""
        ionic = self.ionic
        # ifft (not fft, and no conjugation) gives sum_r n(r)e^(+iGr)/N.
        transposed = np.fft.ifftn(density.reshape(ionic.shape, order="F"))
        sampled = transposed[ionic.local_indices]
        out = np.empty((self.nat, 3), dtype=np.complex128)
        for atom, symbol in enumerate(ionic.species):
            weighted = (sampled * self.local_values[symbol]
                        * np.exp(-1j * (ionic.g @ ionic.positions[:, atom])))
            out[atom] = -1j * self.omega * (weighted @ ionic.g)
        return out.reshape(-1)

    def nonlocal_block(self, left, right):
        """Return tr[(dB D B† + B D dB†) L R†] for all atoms.

        Base overlaps B†L/B†R are computed once per species and reused for all
        three axes; each derivative overlap B†(p_alpha L/R) is computed once.
        Species D is applied locally in its nh channels, never globally.
        """
        out = np.zeros((self.nat, 3), dtype=np.complex128)
        width = left.shape[1]
        for species in self.factors.species:
            nh = species.active_d.shape[0]
            atoms = np.asarray(species.atom_indices, dtype=int) - 1
            shape = (len(atoms), nh, width)
            beta = species.beta
            bl = zgemm(1., beta, left, trans_a=2).reshape(shape)
            br = zgemm(1., beta, right, trans_a=2).reshape(shape)
            dbl = np.einsum("ij,ajb->aib", species.active_d, bl, optimize=True)
            for axis in range(3):
                momentum = self.ionic.p[:, axis, None]
                dleft = zgemm(1., beta, momentum * left, trans_a=2).reshape(shape)
                dright = zgemm(1., beta, momentum * right, trans_a=2).reshape(shape)
                ddleft = np.einsum("ij,ajb->aib", species.active_d, dleft, optimize=True)
                # dB†=+i B†p and dB=-i pB, independently of X's adjoint.
                bra = np.einsum("aib,aib->a", br.conj(), ddleft, optimize=True)
                ket = np.einsum("aib,aib->a", dright.conj(), dbl, optimize=True)
                out[atoms, axis] = 1j * (bra - ket)
        return out.reshape(-1)
