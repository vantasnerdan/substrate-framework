"""Selected-Gamma force/number action factors, with explicit physical measure.

Reuses full saved gPsi; never recomputes ionic actions. Number probes are
integrated spin-resolved particle numbers weighted by cos(b1.r), sin(b1.r),
not local densities or charge. PW states are volume-normalized: the grid
pairing Omega/Ngrid cancels their 1/Omega normalization. No volume or electron
charge multiplies these dimensionless operators. Both collinear spin sectors
are independent inputs, never a degeneracy factor. Only k=q=Gamma is admitted.

Whole-model obligations remain: common physical reference/preparation,
full k/q and region maps, same-background H0/K/curvature dynamics, ionic/joint
state, boundaries/environment and justified dissipation, matter-field inputs,
coherence observables/calibration/design, and nuclear energy-disposal bridge.
This factor interface supplies none of those claims implicitly.
"""
from dataclasses import dataclass
import math
import numpy as np

from .physical_factor_preparation import ActionFactors
from .qe_material_ionic_action import fft_indices


@dataclass(frozen=True)
class NumberProbeGrid:
    """Physical grid/PW identity for projected multiplication by cos/sin(b1.r).

cell_bohr and reciprocal_bohr_inv have lattice vectors as columns, with
cell.T @ reciprocal = 2*pi I. Miller columns and physical FFT bin ordering
match the saved wavefunction. Projection is onto ALL original finite-PW rows,
not saved bands. This is a discrete physical-grid observable, not a new
continuum-cutoff or quadrature convergence assertion.
"""
    fft_shape: tuple
    pw_miller: np.ndarray
    cell_bohr: np.ndarray
    reciprocal_bohr_inv: np.ndarray
    pw_basis_identity: str

    def __post_init__(self):
        shape = tuple(self.fft_shape)
        if (len(shape) != 3 or any(isinstance(n, (bool, np.bool_)) or
                not isinstance(n, (int, np.integer)) or n < 3 for n in shape)):
            raise ValueError('three physical FFT dimensions >=3 required')
        cell = np.asarray(self.cell_bohr, dtype=float)
        reciprocal = np.asarray(self.reciprocal_bohr_inv, dtype=float)
        if (cell.shape != (3,3) or reciprocal.shape != (3,3)
                or not np.isfinite(cell).all() or not np.isfinite(reciprocal).all()
                or not self.pw_basis_identity.strip()):
            raise ValueError('finite physical lattice matrices and PW identity required')
        if not np.allclose(cell.T@reciprocal, 2*np.pi*np.eye(3), rtol=2e-13, atol=2e-13):
            raise ValueError('physical cell and reciprocal columns are not dual')
        indices = fft_indices(self.pw_miller, shape)
        object.__setattr__(self, 'fft_shape', shape)
        object.__setattr__(self, '_indices', indices)
        object.__setattr__(self, 'volume_bohr3', abs(float(np.linalg.det(cell))))
        # Fractional physical coordinates r=cell*(i/n); b1.r=2*pi*i1/n1.
        phase = 2*np.pi*np.arange(shape[0])/shape[0]
        object.__setattr__(self, '_cos', np.cos(phase)[:, None, None])
        object.__setattr__(self, '_sin', np.sin(phase)[:, None, None])

    def actions(self, psi):
        """Return (cos Psi,sin Psi), sharing each band's unitary inverse FFT."""
        psi = np.asarray(psi, dtype=np.complex128)
        if (psi.ndim != 2 or psi.shape[0] != len(self._indices[0])
                or not np.isfinite(psi).all()):
            raise ValueError('finite factors in original full-PW row order required')
        cosine, sine = np.empty_like(psi), np.empty_like(psi)
        grid = np.zeros(self.fft_shape, dtype=complex)
        for band in range(psi.shape[1]):
            grid.fill(0)
            grid[self._indices] = psi[:, band]
            wave = np.fft.ifftn(grid, norm='ortho')
            cosine[:, band] = np.fft.fftn(wave*self._cos, norm='ortho')[self._indices]
            sine[:, band] = np.fft.fftn(wave*self._sin, norm='ortho')[self._indices]
        return cosine, sine

    def metadata(self):
        return {'fft_shape': list(self.fft_shape), 'pw_basis_identity': self.pw_basis_identity,
                'volume_bohr3': self.volume_bohr3,
                'reciprocal_b1_bohr_inv': np.asarray(self.reciprocal_bohr_inv)[:,0].tolist(),
                'grid_measure': 'Omega/Ngrid; PW normalization 1/sqrt(Omega)',
                'operator_measure': 'dimensionless integrated number; volume cancels',
                'charge': 'no -e factor', 'density': 'not number/volume or charge/volume',
                'projection': 'full original finite-PW rows; no saved-band closure'}


def gamma_probe_actions(psi, *, spin_one_based, pw_basis_identity, source_authority,
                        gpsi_h32_x, gpsi_ni1_x, number_grid):
    """Build four Hermitian probes per spin from bound original array inputs.

The two gPsi inputs are the archive's positive ionic gradients in Ry/bohr
for full37atom111-coordinate H32-x/Ni1-x, q=0. They must retain every source
band and full PW row; Root validates immutable archive/row/band receipts.
For these real directions g=g^dagger, so force=-g has the same adjoint action.
Saved MV occupation/energy metadata are deliberately absent from this API.
"""
    psi = np.asarray(psi, dtype=np.complex128)
    if spin_one_based not in (1, 2) or isinstance(spin_one_based, bool):
        raise ValueError('explicit original collinear spin 1 or 2 required')
    if (psi.ndim != 2 or not psi.shape[1] or not np.isfinite(psi).all()
            or not source_authority.strip() or not pw_basis_identity.strip()
            or number_grid.pw_basis_identity != pw_basis_identity):
        raise ValueError('raw factors and matching source/grid/PW authority required')
    result = {}
    scope = '; spin='+str(spin_one_based)+'; k=q=Gamma; full original finite-PW'
    for name, atom, raw in (('force_H32_x', 'H32', gpsi_h32_x),
                            ('force_Ni1_x', 'Ni1', gpsi_ni1_x)):
        raw = np.asarray(raw, dtype=np.complex128)
        if raw.shape != psi.shape or not np.isfinite(raw).all():
            raise ValueError('full gPsi must match every raw factor row and column')
        force = -raw
        identity = '-g('+atom+'-x); full37atom displacement, not normal mode'+scope
        result[name] = ActionFactors(force, force, identity, identity, 'Ry/bohr',
                                     pw_basis_identity, source_authority)
    for kind, action in zip(('cos', 'sin'), number_grid.actions(psi)):
        identity = 'integral d^3r '+kind+'(b1.r) n_spin(r); no charge/volume factor'+scope
        result['number_'+kind+'_b1'] = ActionFactors(action, action, identity, identity,
            '1 (integrated number)', pw_basis_identity, source_authority+'; physical-grid Fourier probe')
    return result


@dataclass(frozen=True)
class CompleteGridNormalization:
    """A chosen complete-grid convention, not a claim that its modes exist.

For N cells, collective centered coordinates sqrt(N)*delta(c^dagger c)
and probe coefficient 1/N give connected sector weight 1/N, not 1/N^2.
Uncentered per-cell means separately have weight 1/N. Raw means remain raw;
the centered collective variable itself has zero mean. No IBZ weight enters.
"""
    shape: tuple
    authority: str

    def __post_init__(self):
        if (len(self.shape) != 3 or any(isinstance(n, (bool, np.bool_)) or
                not isinstance(n, (int, np.integer)) or n < 1 for n in self.shape)
                or not self.authority.strip()):
            raise ValueError('explicit positive complete-grid shape and authority required')

    @property
    def ncells(self):
        return math.prod(int(n) for n in self.shape)

    def metadata(self):
        return {'complete_grid_shape': [int(n) for n in self.shape], 'ncells': self.ncells,
                'authority': self.authority, 'connected_sector_weight': 1/self.ncells,
                'mean_convention': 'uncentered per-cell contribution, weight 1/N',
                'connected_convention': 'collective Gamma: N*(1/N)^2=1/N',
                'actual_coverage': 'selected Gamma only, both spins; NOT complete BvK'}


def gamma_correlations(preparations, left, right, normalization, *, independence_authority):
    """Keep raw sectors and separately weighted PARTIAL electronic results.

The explicit independence authority licenses summing spin covariances with
zero spin cross-cumulants. It is not implied by LSDA input files. No ionic,
electron-ion cross, or total covariance is fabricated. BA is computed on its
own; symmetrized=(AB+BA)/2, commutator=AB-BA (not divided by 2i).
"""
    if set(preparations) != {1, 2} or not independence_authority.strip():
        raise ValueError('both independent spin preparations and their authority required')
    raw = {spin: prep.correlations(left, right) for spin, prep in preparations.items()}
    units = {(item['units_A'], item['units_B']) for item in raw.values()}
    if len(units) != 1:
        raise ValueError('spin sector observable units differ')
    weight = 1/normalization.ncells
    keys = ('mean_A', 'mean_B', 'ordered_AB', 'ordered_BA', 'symmetrized', 'commutator')
    electronic = {key: weight*sum(item[key] for item in raw.values()) for key in keys}
    electronic['units_A'], electronic['units_B'] = next(iter(units))
    return {'raw_sectors': raw, 'electronic': electronic, 'ionic': None, 'cross': None,
            'total': None, 'complete': False, 'scope': 'PARTIAL selected-Gamma equal-time electronic contributions',
            'normalization': normalization.metadata(), 'independence_authority': independence_authority,
            'operator_keys': [left, right], 'stationarity': None, 'dynamics': None, 'noise_spectrum': None}
