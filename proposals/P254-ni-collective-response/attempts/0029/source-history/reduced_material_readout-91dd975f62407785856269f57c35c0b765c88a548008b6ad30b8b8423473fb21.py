"""Spatial pulse and finite kinetic exchange from reduced-model actual states.

Reuse the existing species/displacement metric. Complex energy defects survive;
participation is not an entanglement witness and reversible exchange is not heat.
"""
import sys
from pathlib import Path

import numpy as np

A28 = Path(__file__).resolve().parent.parent/'0028'
if str(A28) not in sys.path:
    sys.path.insert(0, str(A28))
from clamped_scientific_readout import displacement_readout
from substrate_framework.spatial_readout_metrics import periodic_distances


def original_geometry(arrays, mesh, anchor=(0, 0, 0)):
    mesh = np.asarray(mesh, dtype=int)
    cells = np.indices(tuple(mesh)).reshape(3, -1).T
    ids = np.asarray(arrays['atom_indices'], dtype=int)
    names = np.asarray(arrays['atom_names'])
    cell = np.asarray(arrays['cell_bohr'], dtype=float)
    anchor = np.asarray(anchor, dtype=int)
    if (mesh.shape != (3,) or np.any(mesh <= 0) or anchor.shape != (3,)
            or np.any(anchor < 0) or np.any(anchor >= mesh)
            or names.shape != ids.shape or len(set(ids)) != len(ids)):
        raise ValueError('original atom identities and periodic geometry required')
    return dict(mesh=mesh, cells=cells, ids=ids, names=names, cell=cell,
                anchor=anchor, n=len(cells), d=3*len(ids))


def spatial_pulse(delta_q, geometry):
    """Read all original coordinates; zero signal remains undefined, not clipped."""
    d, n = geometry['d'], geometry['n']
    if len(delta_q) != n:
        raise ValueError('every original q block is required for spatial readout')
    ions_q = np.stack([np.asarray(block)[:2*d] for block in delta_q])
    ions = np.fft.ifftn(ions_q.reshape(*geometry['mesh'], 2*d),
                       axes=(0, 1, 2), norm='ortho')
    distances = periodic_distances(geometry['cells'], geometry['anchor'],
                                   geometry['mesh'], geometry['cell'])
    return dict(ions=ions, displacement=displacement_readout(ions, geometry, distances),
                definition='actual kicked-minus-unexcited state; orthonormal inverse FFT',
                scope='Conditional active-electron model, not full-PW material accuracy.')


def finite_kinetic_exchange(baseline_gamma, delta_q, geometry, masses, charges,
                            field_mass):
    """Exact baseline-delta kinetic differences, never two large-energy subtraction.

    Gamma ionic coordinates are orthonormal Fourier amplitudes. Maxwell a,Pi
    are one physical global field, not per-cell or Fourier-rescaled copies.
    """
    d, n = geometry['d'], geometry['n']
    masses = np.asarray(masses).reshape(d)
    charges = np.asarray(charges).reshape(-1)
    if len(charges)*3 != d or np.any(masses <= 0) or field_mass <= 0:
        raise ValueError('original positive masses/charges and Maxwell mass required')
    baseline_gamma = np.asarray(baseline_gamma)
    delta_gamma = np.asarray(delta_q[0])
    if baseline_gamma.shape != (2*d+6,) or delta_gamma.shape != (2*d+6,):
        raise ValueError('original Gamma (u,p,a,Pi) order with three field modes required')
    ions_q = np.stack([np.asarray(block)[:2*d] for block in delta_q])
    delta_p = np.fft.ifftn(ions_q.reshape(*geometry['mesh'], 2*d),
                          axes=(0, 1, 2), norm='ortho').reshape(n, 2*d)[:, d:]
    z = np.repeat(charges, 3)
    a0, da = baseline_gamma[2*d:2*d+3], delta_gamma[2*d:2*d+3]
    baseline_mechanical = baseline_gamma[d:2*d]/np.sqrt(n)-z*np.tile(a0, d//3)
    delta_mechanical = delta_p-z[None, :]*np.tile(da, d//3)[None, :]
    energy = (baseline_mechanical[None, :]*delta_mechanical
              +delta_mechanical*delta_mechanical/2)/masses[None, :]
    by_cell_atom = energy.reshape(n, -1, 3).sum(axis=2)
    by_atom = by_cell_atom.sum(axis=0)
    pi0, dpi = baseline_gamma[2*d+3:], delta_gamma[2*d+3:]
    field = (pi0@dpi+dpi@dpi/2)/field_mass
    return dict(ionic_by_cell_atom_ry=by_cell_atom,
                ionic_by_original_atom_ry={str(i): v for i, v in zip(geometry['ids'], by_atom)},
                ionic_by_species_ry={name: by_atom[geometry['names'] == name].sum()
                                     for name in dict.fromkeys(geometry['names'])},
                ionic_total_ry=by_atom.sum(), maxwell_free_mean_change_ry=field,
                scope='Finite mean kinetic exchange; no covariance energy, dissipation or heat inferred.')
