"""Exact initial jet of the existing canonical total-force observable.

For a canonical momentum-only seed, all initial electronic, displacement and
field variations vanish. The coupled RHS then has u'=p/M, electronic'=0 and
field-coordinate'=0. Thus F(0)=0 and F'(0)=-C(q)p/M, even on a moving mean.
This is an initial derivative, NOT a finite-time frozen-electron propagator.
"""
import numpy as np
from substrate_framework.joint_spatial_readout import complete_grid


def force_rate(curvature, masses, q, *, ncells, cell, amplitude, coordinate):
    """Original collective Fourier force derivative, retaining complex C(q)."""
    c = np.asarray(curvature, complex)
    m = np.asarray(masses)
    q, cell = np.asarray(q, float), np.asarray(cell)
    if (c.ndim != 2 or c.shape[0] != c.shape[1] or m.shape != (len(c),)
            or not np.isfinite(c).all() or not np.isfinite(m).all()
            or np.any(m.imag != 0) or np.any(m.real <= 0)):
        raise ValueError('finite original square curvature and positive real Cartesian masses required')
    if (type(ncells) is not int or ncells < 1 or q.shape != (3,)
            or cell.shape != (3,) or not np.issubdtype(cell.dtype, np.integer)
            or not np.isfinite(q).all() or not np.isfinite(amplitude)
            or not 0 <= coordinate < len(m)):
        raise ValueError('original q/cell, positive N, finite seed and Cartesian coordinate required')
    phase = np.exp(-2j*np.pi*q @ cell)
    return -c[:, coordinate]*(amplitude*phase/(np.sqrt(ncells)*m[coordinate].real))


def full_grid(curvatures, masses, *, mesh, cell, amplitude, coordinate):
    """Complete original character grid; return only the initial derivative."""
    mesh = np.asarray(mesh)
    q = np.indices(tuple(mesh)).reshape(3, -1).T/mesh
    q, cells = complete_grid(q, mesh)
    c = np.asarray(curvatures)
    d = len(masses)
    if c.shape != (*mesh, d, d):
        raise ValueError('complete original q curvature grid required')
    values = np.asarray([force_rate(matrix, masses, point, ncells=len(q), cell=cell,
        amplitude=amplitude, coordinate=coordinate) for matrix, point in zip(c.reshape(-1, d, d), q)])
    real_space = np.fft.ifftn(values.reshape(*mesh, d), axes=(0, 1, 2), norm='ortho').reshape(len(q), d)
    return dict(q_fractional=q, cells=cells, force_rate_q=values, force_rate_cells=real_space,
        units='Ry^2/bohr: derivative of probe-defined total force (Ry/bohr) with respect to time in Ry^-1',
        scope='Exact initial jet of original coupled scalar and momentum-only seed; no finite-time trajectory, remainder bound, quantum covariance, coherence window, heat or lattice ranking.')
