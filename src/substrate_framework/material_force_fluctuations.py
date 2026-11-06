"""Initial full-PW electronic force noise in an explicitly supplied raw state.

For real dimensionless atom-major directions u, F_u=-u.dH/dR (Ry/bohr).
The state is exactly rho=Psi D Psi† with empty initial orthogonal complement.
T_AB=(A†Psi)†(B Psi) retains transitions through that complement. No saved-band
projection, frame normalization, state selection, thermal or dynamical model
occurs here. Results are ONE UNWEIGHTED spin/k sector; no implicit spin factor.
"""
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from types import MappingProxyType

import numpy as np
from scipy.linalg.blas import zgemm

from .physical_factor_preparation import FactorFunctional, FactorPreparation


def _readonly(value):
    value.flags.writeable = False
    return value


def _cross_gram(left, right, row_block_size):
    """Full-row contraction with only row-block-sized BLAS layout copies."""
    out = np.zeros((left.shape[1], right.shape[1]), dtype=np.complex128, order='F')
    for start in range(0, left.shape[0], row_block_size):
        stop = min(start + row_block_size, left.shape[0])
        a = np.asfortranarray(left[start:stop], dtype=np.complex128)
        b = a if left is right else np.asfortranarray(right[start:stop], dtype=np.complex128)
        out = zgemm(1., a, b, trans_a=2, beta=1., c=out, overwrite_c=1)
    return _readonly(out)


@dataclass(frozen=True)
class ForceFluctuations:
    preparation: FactorPreparation
    directions: object
    coordinate_identity: str
    gram_authority: str

    def archive_payload(self):
        """Small matrices plus all directed connected and unconnected moments.

        Matrix indices are operator_order, including number. C[i,j] is AB,
        never BA or A†B. reverse_ordered_connected[i,j] is C[j,i] WITHOUT
        conjugation. Raw complex values and diagnostic imaginary parts remain.
        """
        prep = self.preparation
        arrays, metadata = prep.functional.archive_payload()
        names = tuple(prep.functional.overlaps)
        means = np.asarray([prep.mean(name) for name in names])
        connected = np.asarray([[prep.ordered(a, b) for b in names] for a in names])
        second = connected + means[:, None] * means[None, :]
        arrays.update(density_coefficients=prep.density_coefficients,
                      means=means, ordered_connected=connected,
                      reverse_ordered_connected=connected.T,
                      ordered_second_moment=second,
                      reverse_ordered_second_moment=second.T,
                      symmetrized_connected=(connected + connected.T)/2,
                      commutator=connected-connected.T,
                      directions=np.stack(tuple(self.directions.values())))
        # This diagnostic is the term lost by projection, not a replacement
        # for the retained product. Only raw-factor-sized solves are formed.
        f, d = prep.functional, prep.density_coefficients
        solved = {name: np.linalg.solve(f.gram, f.overlaps[name]) for name in names}
        discarded = np.asarray([[np.trace(d @ (f.ordered_products[a, b]
                                - f.overlaps[a] @ solved[b])) for b in names] for a in names])
        arrays['empty_complement_contribution'] = discarded
        metadata.update(state='explicit gauge-invariant quasifree source checkpoint',
            state_authority=prep.state_authority, admissibility=prep.admissibility,
            density_coefficients='density_coefficients', operator_order=list(names),
            direction_order=list(self.directions), coordinate_identity=self.coordinate_identity,
            gram_authority=self.gram_authority,
            force_sign='F_u=-u.dH/dR; u real dimensionless Cartesian atom-major direction',
            normalization='one unweighted spin/k sector; no spin degeneracy or k measure applied',
            moment_layout={'means': '<A>', 'ordered_connected': '<delta A delta B>',
                'reverse_ordered_connected': '<delta B delta A>',
                'ordered_second_moment': '<AB>', 'reverse_ordered_second_moment': '<BA>',
                'symmetrized_connected': '<{delta A,delta B}>/2',
                'commutator': '<[A,B]>',
                'empty_complement_contribution': 'Tr[D(T_AB-M_A S^-1 M_B)]; retained in C_AB'},
            scope='initial electronic equal-time noise only; no equilibrium, FDT, damping or two-time kernel')
        return arrays, metadata


def initial_force_fluctuations(psi, density_coefficients, coupling, directions, *,
        factor_identity, pw_basis_identity, state_authority, operator_authority,
        coordinate_identity, gram=None, gram_authority=None, scratch_directory=None,
        block_size=8, row_block_size=2048):
    """Compose physical ionic actions with existing FactorFunctional/Preparation.

    coupling supplies npw, ndof and ionic.prepare(direction[nat,3]).apply(block),
    the existing POSITIVE ionic gradient API. Prepare each direction only once.
    All source columns are consumed; the caller binds the exact source/frame/
    coupling and supplies D without refilling it.
    Optional gram must be the authenticated Gram of THIS unchanged raw frame;
    gram_authority must identify that receipt, not merely assert orthogonality.
    The existing public FactorFunctional constructor is the cached-Gram path.

    Disk-backed actions hold one full frame per force direction. FFT/action
    temporaries are column-bounded, contraction temporaries row-bounded; no
    PW-square object or full-frame conjugate/layout copy is allocated here.
    Temporary action files are removed on success or exception.
    """
    if (type(block_size) is not int or not 1 <= block_size <= 16
            or type(row_block_size) is not int or row_block_size < 1):
        raise ValueError('block_size in 1..16 and positive integer row_block_size required')
    for value in (factor_identity, pw_basis_identity, state_authority,
                  operator_authority, coordinate_identity):
        if not isinstance(value, str) or not value.strip():
            raise ValueError('explicit source, basis, state, operator and coordinate identities required')
    psi = np.asarray(psi)
    if (psi.ndim != 2 or psi.shape[0] != coupling.npw
            or not 0 < psi.shape[1] <= psi.shape[0] or psi.dtype.kind != 'c'):
        raise ValueError('full-PW complex raw frame with every source column required')
    for start in range(0, psi.shape[0], row_block_size):
        if not np.isfinite(psi[start:start+row_block_size]).all():
            raise ValueError('nonfinite source frame')
    selected = {}
    if not directions:
        raise ValueError('at least one explicitly selected physical force direction required')
    for name, value in directions.items():
        if not isinstance(name, str) or not name.strip() or name == 'number':
            raise ValueError('nonempty distinct force names excluding reserved number required')
        value = np.asarray(value)
        if (np.iscomplexobj(value) or value.shape != (coupling.ndof,)
                or not np.isfinite(value).all() or not np.any(value)):
            raise ValueError('finite nonzero REAL atom-major physical direction required')
        selected[name] = _readonly(np.array(value, dtype=float, copy=True))
    if gram is None:
        if gram_authority is not None:
            raise ValueError('Gram authority supplied without Gram')
        s = _cross_gram(psi, psi, row_block_size)
        gram_authority = 'computed full-row Gram from '+factor_identity
    else:
        if not isinstance(gram_authority, str) or not gram_authority.strip():
            raise ValueError('authenticated same-frame Gram receipt authority required')
        s = np.array(gram, dtype=np.complex128, copy=True)
        if s.shape != (psi.shape[1], psi.shape[1]):
            raise ValueError('authenticated Gram must retain every source column')
        _readonly(s)
    metadata = {name: MappingProxyType(dict(operator_identity='force:'+name,
        adjoint_identity='force:'+name, units='Ry/bohr',
        pw_basis_identity=pw_basis_identity, source_authority=operator_authority)) for name in selected}
    metadata['number'] = MappingProxyType(dict(operator_identity='full-PW-identity',
        adjoint_identity='full-PW-identity', units='dimensionless electron count',
        pw_basis_identity=pw_basis_identity, source_authority=operator_authority))
    overlaps, products = {}, {}
    actions = {}
    with TemporaryDirectory(prefix='force-actions-', dir=scratch_directory) as directory:
        try:
            for index, (name, direction) in enumerate(selected.items()):
                out = np.memmap(Path(directory)/('action-'+str(index)+'.bin'),
                                mode='w+', dtype=np.complex128, shape=psi.shape, order='F')
                actions[name] = out
                prepared = coupling.ionic.prepare(direction.reshape(-1, 3))
                for start in range(0, psi.shape[1], block_size):
                    stop = min(start+block_size, psi.shape[1])
                    positive = np.asarray(prepared.apply(psi[:, start:stop]))
                    if positive.shape != (psi.shape[0], stop-start) or not np.isfinite(positive).all():
                        raise ValueError('ionic gradient action lost rows/columns or returned nonfinite values')
                    np.negative(positive, out=out[:, start:stop])
                out.flush()
                overlaps[name] = _cross_gram(psi, out, row_block_size)
                for previous, other in actions.items():
                    products[previous, name] = _cross_gram(other, out, row_block_size)
                    if previous != name:
                        # Exact adjoint identity for declared real Hermitian
                        # directions, not numerical Hermitization.
                        products[name, previous] = _readonly(products[previous, name].conj().T)
                products['number', name] = overlaps[name]
                products[name, 'number'] = _readonly(overlaps[name].conj().T)
            overlaps['number'] = s
            products['number', 'number'] = s
        finally:
            for action in actions.values():
                action._mmap.close()
    functional = FactorFunctional(s, MappingProxyType(overlaps), MappingProxyType(products),
        MappingProxyType(metadata), factor_identity, pw_basis_identity, psi.shape[0])
    preparation = FactorPreparation(functional, density_coefficients, state_authority)
    return ForceFluctuations(preparation, MappingProxyType(selected), coordinate_identity, gram_authority)
