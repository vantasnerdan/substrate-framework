"""Full physical ionic support closed between canonical opposite q edges.

No projected coefficient is averaged or repaired. The original density Miller
support is retained and extended by its exact physical-wavevector partner.
Aliased grid modes, if present after extension, are ADDITIVELY accumulated;
they remain distinct physical |q+G| radial/derivative contributions.
"""
from dataclasses import dataclass
import numpy as np
from substrate_framework.qe_nc_radial import build_radial_tables


@dataclass(frozen=True)
class PairedLocalSupport:
    physical: tuple
    miller: np.ndarray
    reciprocal_bohr_inv: np.ndarray
    q_bohr_inv: np.ndarray
    omega: float
    nl: np.ndarray
    wavevectors: np.ndarray
    nnr: int
    original_count: int

    def to_g(self, field):
        field = np.asarray(field)
        spectrum = np.fft.fftn(field.reshape((*self.physical, *field.shape[1:]), order='F'),
                               axes=(0, 1, 2), norm='forward')
        return spectrum.reshape((self.nnr, *field.shape[1:]), order='F')[self.nl]

    def from_g(self, coefficients):
        coefficients = np.asarray(coefficients)
        if coefficients.shape[0] != len(self.nl):
            raise ValueError('paired physical support coefficient count differs')
        spectrum = np.zeros((self.nnr, *coefficients.shape[1:]), complex)
        np.add.at(spectrum, self.nl, coefficients)
        field = np.fft.ifftn(spectrum.reshape((*self.physical, *coefficients.shape[1:]), order='F'),
                             axes=(0, 1, 2), norm='forward')
        return field.reshape(spectrum.shape, order='F')


def paired_millers(original, q_fractional, minus_fractional):
    original = np.asarray(original)
    if original.ndim != 2 or original.shape[0] != 3 or original.dtype.kind not in 'iu':
        raise ValueError('original integer physical Miller columns required')
    shift = np.asarray(q_fractional)+np.asarray(minus_fractional)
    integer = np.rint(shift).astype(int)
    if not np.allclose(shift, integer, rtol=0, atol=2e-10):
        raise ValueError('canonical opposite q representatives must sum to an integer')
    columns = [tuple(g) for g in original.T]
    seen = set(columns)
    for g in original.T:
        partner = tuple(-g-integer)
        if partner not in seen:
            seen.add(partner); columns.append(partner)
    return np.asarray(columns, dtype=int).T


def create_support(original_fft, q_fractional, minus_fractional):
    miller = paired_millers(original_fft.miller, q_fractional, minus_fractional)
    shape = original_fft.physical
    q_cart = original_fft.reciprocal_bohr_inv @ np.asarray(q_fractional)
    nl = np.ravel_multi_index(tuple(miller[a] % shape[a] for a in range(3)), shape, order='F')
    return PairedLocalSupport(shape, miller, original_fft.reciprocal_bohr_inv, q_cart,
        original_fft.omega, nl, original_fft.reciprocal_bohr_inv @ miller+q_cart[:, None],
        original_fft.nnr, original_fft.miller.shape[1])


def attach_paired_local_support(ionic, minus_fractional):
    """Extend the physical action itself before any action/projection is made."""
    support = create_support(ionic.bridge.density_fft, ionic.edge.q_fractional, minus_fractional)
    ionic.bridge.density_fft = support
    ionic.gq = support.wavevectors.T
    moduli = np.linalg.norm(ionic.gq, axis=1)
    required = float(np.max(moduli, initial=0.))
    tables, vloc = {}, {}
    for name in dict.fromkeys(ionic.inputs.reference.atom_names):
        original = ionic.inputs.tables[name]
        current = ionic.tables[name]
        if current.qmax_bohr_inv < required:
            current = build_radial_tables(original.upf, original.omega_bohr3, required)
        if current.upf is not original.upf or current.omega_bohr3 != original.omega_bohr3:
            raise ValueError('paired support must retain exact original UPF and cell volume')
        tables[name] = current
        vloc[name] = current.vloc(moduli)
    ionic.tables = tables; ionic.vloc = vloc
    return dict(original_count=support.original_count, paired_count=len(support.nl),
                aliased_additive_modes=int(len(support.nl)-len(np.unique(support.nl))),
                support_rule='original union (-original-q-qminus); original UPF radial extension')
