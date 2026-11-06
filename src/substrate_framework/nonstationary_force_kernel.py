"""Ordered nonstationary quasifree force kernels with the full PW complement.

For prescribed number-conserving quadratic evolution U(t,0), use Heisenberg
factors X_A=U† A(t) U Psi, X_Aadj=U† A(t)† U Psi. Both arrays retain every
physical PW row and raw column in the INITIAL frame's Hilbert space. Merely
supplying A(t) Psi(t) is wrong at unequal times: it omits the return U†.
C_AB(t,s)=Tr(D X_Aadj(t)† X_B(s))-Tr(D M_A(t) D M_B(s)).
No stationarity, FDT, bath temperature, Markov approximation or fitted loss.
This is exact for the declared quasifree quadratic path, not an assertion
that interacting higher cumulants or self-consistent screened response equal
the bare propagated kernel. PreparedBlochResponse owns moving-mean screening.
For delta H=+f B, chi_AB=-i theta(t-s)<[A(t),B(s)]>; for -u F flip sign.
"""
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np

from .physical_factor_preparation import FactorFunctional, FactorPreparation
from .prepared_reference import array, authority, frozen


def _cross(left, right, row_block_size):
    out = np.zeros((left.shape[1], right.shape[1]), complex)
    for start in range(0, left.shape[0], row_block_size):
        end = min(start+row_block_size, left.shape[0])
        out += left[start:end].conj().T @ right[start:end]
    return out


@dataclass(frozen=True)
class HeisenbergForceFactors:
    time_inverse_ry: float
    action: np.ndarray
    adjoint_action: np.ndarray
    operator_identity: str
    units: str
    initial_factor_identity: str
    pw_basis_identity: str
    propagation_authority: str

    def validate(self, functional, row_block_size):
        for name in ('operator_identity', 'units', 'initial_factor_identity',
                     'pw_basis_identity', 'propagation_authority'):
            authority(getattr(self, name), name)
        if not np.isfinite(self.time_inverse_ry):
            raise ValueError('finite absolute time in inverse Ry required')
        if (self.initial_factor_identity != functional.factor_identity
                or self.pw_basis_identity != functional.pw_basis_identity):
            raise ValueError('Heisenberg factors must return to the admitted INITIAL frame/basis')
        shape = (functional.full_pw_dimension, len(functional.gram))
        for value in (self.action, self.adjoint_action):
            if not isinstance(value, np.ndarray) or value.shape != shape or value.dtype.kind != 'c':
                raise ValueError('complex full-PW actions retaining ALL raw columns required')
            for start in range(0, shape[0], row_block_size):
                if not np.isfinite(value[start:start+row_block_size]).all():
                    raise ValueError('nonfinite propagated action')


@dataclass(frozen=True)
class OrderedForcePair:
    time_a: float
    time_b: float
    mean_a: complex
    mean_b: complex
    ordered_ab: complex
    ordered_ba: complex

    @property
    def symmetrized(self):
        return (self.ordered_ab+self.ordered_ba)/2

    @property
    def commutator(self):
        return self.ordered_ab-self.ordered_ba

    def retarded(self, *, perturbation_sign, equal_time_theta=0.5):
        """Response to delta H=sign*f*B, sign explicitly +1 or -1.

        A force-displacement perturbation -u F uses -1, not the +a L
        optical sign. Equal-time theta is an explicit distribution convention;
        instantaneous contact terms are separate derivatives of the scalar.
        """
        if perturbation_sign not in (-1, 1) or equal_time_theta not in (0, 0.5, 1):
            raise ValueError('explicit +/-1 perturbation sign and theta(0)=0,1/2,1 required')
        theta = (1. if self.time_a > self.time_b else
                 0. if self.time_a < self.time_b else equal_time_theta)
        return -1j*perturbation_sign*theta*self.commutator


class NonstationaryForceKernel:
    def __init__(self, preparation, *, row_block_size=2048):
        if not isinstance(preparation, FactorPreparation):
            raise TypeError('admitted FactorPreparation required; D is not inferred')
        if type(row_block_size) is not int or row_block_size < 1:
            raise ValueError('positive row block required')
        self.preparation = preparation
        self.row_block_size = row_block_size
        self.initial_covariance = None
        self.initial_means = None
        self.initial_names = None
        self.checkpoint_authority = None

    @classmethod
    def from_initial_checkpoint(cls, arrays, metadata, *, checkpoint_authority,
                                tolerance=1e-10):
        """Consume actual stored equal-time covariance; never redo force actions.

        Receipt/byte authentication belongs to caller. Validate the stored
        quantum covariance against the archived exact S,D,M,T interface and
        all directed/symmetrized/commutator companions. Return ORIGINAL stored
        covariance after validation, not a substituted recomputation or repair.
        """
        authority(checkpoint_authority, 'authenticated initial checkpoint')
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('positive finite comparison tolerance required')
        if (metadata.get('schema') != 'material-equal-time-force-noise-v1'
                or metadata.get('all_raw_columns') is not True
                or metadata.get('full_physical_pw_rows') is not True
                or metadata.get('source_slot_measure_applied') is not False):
            raise ValueError('actual unweighted full-row full-column initial material noise required')
        layout = metadata['array_layout']
        names = tuple(metadata['operator_order'])
        if not names or len(set(names)) != len(names) or set(names) != set(layout['overlaps']):
            raise ValueError('unique complete archived operator order required')
        overlaps = {name: arrays[layout['overlaps'][name]] for name in names}
        products = {(entry['left'], entry['right']): arrays[entry['array']]
                    for entry in layout['ordered_products']}
        if set(products) != {(a, b) for a in names for b in names}:
            raise ValueError('both ordered products for every pair required')
        gram = np.asarray(arrays[layout['gram']])
        shape = gram.shape
        if shape != (metadata['factor_count'], metadata['factor_count']):
            raise ValueError('archived Gram/raw column count mismatch')
        for name, value in list(overlaps.items())+list(products.items()):
            array(value, shape, 'archived overlap/product '+str(name))
        f = FactorFunctional(gram, MappingProxyType(overlaps), MappingProxyType(products),
            MappingProxyType(metadata['operators']), metadata['factor_identity'],
            metadata['pw_basis_identity'], metadata['full_pw_dimension'])
        prep = FactorPreparation(f, arrays['density_coefficients'], metadata['state_authority'])
        kernel = cls(prep)
        c = array(arrays['ordered_connected'], (len(names), len(names)), 'stored C_AB')
        means = array(arrays['means'], (len(names),), 'stored means')
        expected = np.asarray([[prep.ordered(a, b) for b in names] for a in names])
        def require(actual, reference, label):
            actual = array(actual, np.shape(reference), label)
            if not np.allclose(actual, reference, rtol=tolerance, atol=tolerance):
                raise ValueError('stored quantum kernel mismatch: '+label)
        require(c, expected, 'C_AB versus Tr(D T)-Tr(D M D M)')
        require(means, np.asarray([prep.mean(a) for a in names]), 'means')
        require(arrays['reverse_ordered_connected'], c.T, 'BA transpose without conjugation')
        require(arrays['symmetrized_connected'], (c+c.T)/2, 'anticommutator')
        require(arrays['commutator'], c-c.T, 'commutator')
        require(arrays['ordered_second_moment'], c+means[:, None]*means[None, :], 'second moment')
        require(arrays['reverse_ordered_second_moment'],
                (c+means[:, None]*means[None, :]).T, 'reverse second moment')
        # The actual checkpoint declares real Hermitian force/number operators.
        for name in names:
            op = metadata['operators'][name]
            if op['operator_identity'] != op['adjoint_identity']:
                raise ValueError('initial force checkpoint must declare Hermitian operators')
        require(c, c.conj().T, 'Hermitian covariance; no Hermitization')
        if np.linalg.eigvalsh(c).min() < -tolerance*(1+np.linalg.norm(c)):
            raise ValueError('stored Hermitian force covariance is not positive semidefinite')
        kernel.initial_covariance = frozen(c)
        kernel.initial_means = frozen(means)
        kernel.initial_names = names
        kernel.checkpoint_authority = checkpoint_authority
        return kernel

    def initial_pair(self, left, right):
        if self.initial_covariance is None:
            raise ValueError('actual authenticated equal-time checkpoint has not been supplied')
        i, j = self.initial_names.index(left), self.initial_names.index(right)
        return OrderedForcePair(0., 0., self.initial_means[i], self.initial_means[j],
                                self.initial_covariance[i, j], self.initial_covariance[j, i])

    def pair(self, initial_frame, left, right):
        """Evaluate supplied unequal-time full-row factors, independent AB and BA.

        Cost per pair O(P*r^2+r^3), storage O(r^2+row_block*r) beyond supplied
        factors. Storing n times and m probes costs up to 2*n*m*P*r complex
        values; stream pairs instead if the propagation service can replay.
        """
        prep, block = self.preparation, self.row_block_size
        f, d = prep.functional, prep.density_coefficients
        left.validate(f, block)
        right.validate(f, block)
        if left.propagation_authority != right.propagation_authority:
            raise ValueError('both endpoints must belong to the same supplied unitary history')
        psi = np.asarray(initial_frame)
        shape = (f.full_pw_dimension, len(f.gram))
        if psi.shape != shape or psi.dtype.kind != 'c':
            raise ValueError('unchanged initial full-PW frame with ALL raw columns required')
        for start in range(0, len(psi), block):
            if not np.isfinite(psi[start:start+block]).all():
                raise ValueError('nonfinite initial frame')
        if not np.allclose(_cross(psi, psi, block), f.gram, rtol=1e-12, atol=1e-12):
            raise ValueError('initial raw frame differs from the admitted Gram; no normalization repair')
        ma, mb = (_cross(psi, endpoint.action, block) for endpoint in (left, right))
        tab = _cross(left.adjoint_action, right.action, block)
        tba = _cross(right.adjoint_action, left.action, block)
        ab = np.trace(d @ tab)-np.trace(d @ ma @ d @ mb)
        ba = np.trace(d @ tba)-np.trace(d @ mb @ d @ ma)
        return OrderedForcePair(left.time_inverse_ry, right.time_inverse_ry,
            complex(np.trace(d @ ma)), complex(np.trace(d @ mb)), complex(ab), complex(ba))
