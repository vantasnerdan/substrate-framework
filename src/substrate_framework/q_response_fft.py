"""Execution-only physical FFT with one owned multidimensional work buffer.

Preserves PhysicalRhoFFT support, flattening and normalization. SciPy's single
multidimensional transform can overwrite this private buffer, unlike NumPy's
sequence of axis transforms. Inputs and source/model authority are unchanged.
The response owner opts in explicitly; saved constructors are not modified.
"""
from dataclasses import fields

import numpy as np
from scipy import fft as spfft

from .qe_material_grid import PhysicalRhoFFT


class ResponseRhoFFT(PhysicalRhoFFT):
    @classmethod
    def from_physical(cls, original):
        if not isinstance(original, PhysicalRhoFFT):
            raise TypeError('an existing checked physical FFT is required')
        return cls(**{field.name: getattr(original, field.name) for field in fields(original)})

    def to_g(self, field):
        field = np.asarray(field)
        # Never permit overwrite_x to modify the caller's actual density/frame.
        work = np.array(self.physical_view(field), dtype=np.complex128,
                        order='F', copy=True)
        spectrum = spfft.fftn(work, axes=(0, 1, 2), norm='forward',
                             overwrite_x=True, workers=1)
        return spectrum.reshape((self.nnr, *field.shape[1:]), order='F')[self.nl]

    def from_g(self, coefficients):
        coefficients = np.asarray(coefficients)
        if coefficients.shape[0] != self.nl.size:
            raise ValueError('coefficient rows must match the explicit Miller map')
        # F order makes physical_view a view, not another full-grid allocation.
        work = np.zeros((self.nnr, *coefficients.shape[1:]), dtype=np.complex128,
                        order='F')
        work[self.nl] = coefficients
        result = spfft.ifftn(self.physical_view(work), axes=(0, 1, 2),
                            norm='forward', overwrite_x=True, workers=1)
        return result.reshape(work.shape, order='F')
