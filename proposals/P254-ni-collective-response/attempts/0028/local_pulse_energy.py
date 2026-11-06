"""Initial canonical local-pulse input, not a propagated energy partition.

For U=exp(+i J x/hbar), U^dagger p U=p+J.  Only the struck
Cartesian kinetic term changes: J*(<p>-Z*<a>)/M + J^2/(2*M).
This identity needs only finite first moments and kinetic energy, not a
Gaussian or equilibrium.  It leaves covariances and coordinate-only terms
unchanged.  Fourier rows below are the centered preparation's quadratic
input shares, not energies of independently propagated complex q states.
"""
from cmath import exp
from itertools import product
from math import pi, prod, sqrt


def pulse_components(mass, charge, amplitude, canonical_mean=0, vector_potential_mean=0):
    """Scalar struck-axis formula; supports exact symbolic scalars as well."""
    mechanical_mean = canonical_mean-charge*vector_potential_mean
    cross = amplitude*mechanical_mean/mass
    quadratic = amplitude**2/(2*mass)
    return dict(mechanical_mean=mechanical_mean, linear_cross=cross,
                quadratic=quadratic, total=cross+quadratic)


def fourier_pulse(mesh, cell, amplitude):
    """Complete orthonormal forward DFT of J delta_(R,cell), without an FFT.

    q is fractional reciprocal-cell coordinates.  All q and -q are retained;
    the inverse is sum_q exp(+2pi i q.R) Pq/sqrt(N), without a second dagger.
    """
    mesh, cell = tuple(mesh), tuple(cell)
    if (len(mesh) != 3 or len(cell) != 3
            or any(type(n) is not int or n < 1 for n in mesh)
            or any(type(r) is not int for r in cell)):
        raise ValueError('three positive integer mesh sizes and integer cell indices required')
    ncells = prod(mesh)
    rows = []
    for step in product(*(range(n) for n in mesh)):
        q = tuple(k/n for k, n in zip(step, mesh))
        coefficient = amplitude*exp(-2j*pi*sum(k*r/n for k, r, n in zip(step, cell, mesh)))/sqrt(ncells)
        rows.append(dict(mesh_index=step, q_fractional=q, coefficient=coefficient,
                         quadratic_weight=(coefficient.conjugate()*coefficient).real))
    return rows
