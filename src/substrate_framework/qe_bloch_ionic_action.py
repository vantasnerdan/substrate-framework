"""General-k q=0 bare ordinary-NC ionic derivative, Ry/bohr.

Reuses checked k+G projector factors and the accepted two-sided complex-linear
DirectionalAction algebra. GammaIonicAction remains unchanged. Recompute each
physical atom direction on target geometry; no Gamma archive rotation occurs.
"""
from dataclasses import dataclass
import numpy as np
from .qe_bloch_qzero import BlochQZeroBridge
from .qe_material_hamiltonian import build_pw_factors_for_basis
from .qe_material_ionic_action import DirectionalAction


@dataclass(frozen=True)
class BlochDirectionalAction(DirectionalAction):
    def adjoint_apply(self, psi):
        """Hilbert adjoint at q=0: conjugate the direction, not its apply input."""
        terms = tuple((pu.conj(), bd, adjoint) for pu, bd, adjoint in self.terms)
        adjoint = BlochDirectionalAction(self.owner, self.direction.conj(), self.field.conj(), terms)
        return adjoint.apply(psi)


class BlochIonicAction:
    def __init__(self, inputs, basis, *, factors=None):
        self.bridge = BlochQZeroBridge(inputs, basis)
        self.inputs, self.basis = inputs, basis
        self.factors = build_pw_factors_for_basis(inputs, basis) if factors is None else factors
        provenance = self.factors.sources[-1]
        if (provenance.get('identity') != basis.identity or self.factors.sources[:-1] != inputs.sources
                or not np.array_equal(self.factors.fft.miller, basis.miller)
                or not np.array_equal(self.factors.k_bohr_inv, basis.k_cart_bohr_inv)
                or not np.array_equal(self.factors.fft.reciprocal_bohr_inv, basis.reciprocal_bohr_inv)):
            raise ValueError('ionic action requires factors of the same checked basis/inputs')
        self.p = (basis.reciprocal_bohr_inv @ basis.miller + basis.k_cart_bohr_inv[:, None]).T
        self.g = self.bridge.density_fft.wavevectors.T  # q=0 density G, not k+G
        self.atoms, self.vloc = {}, {}

    def local_action(self, field, psi):
        return self.bridge.local_apply(field, psi)

    def _atom(self, atom):
        if atom not in self.atoms:
            for species in self.factors.species:
                if atom+1 in species.atom_indices:
                    index = species.atom_indices.index(atom+1)
                    nh = species.active_d.shape[0]
                    beta = species.beta[:, index*nh:(index+1)*nh]
                    self.atoms[atom] = (beta @ species.active_d, beta.conj().T)
                    break
            else:
                raise ValueError('atom is absent from checked NC projector factors')
        return self.atoms[atom]

    def prepare(self, direction):
        """Complex-linear dimensionless direction[nat,3]; epsilon has units bohr."""
        ref = self.inputs.reference
        u = np.array(direction, dtype=np.complex128, copy=True)
        if u.shape != (ref.nat, 3) or not np.isfinite(u).all():
            raise ValueError('finite full original atom-order direction required')
        local = np.zeros(len(self.g), dtype=np.complex128)
        terms = []
        for atom in np.flatnonzero(np.any(u != 0, axis=1)):
            name = ref.atom_names[atom]
            if name not in self.vloc:
                self.vloc[name] = self.inputs.tables[name].vloc(np.linalg.norm(self.g, axis=1))
            local += (-1j*(self.g @ u[atom])*self.vloc[name]
                      * np.exp(-1j*(self.g @ ref.positions_bohr[:, atom])))
            bd, adjoint = self._atom(int(atom))
            terms.append((self.p @ u[atom], bd, adjoint))
        field = self.bridge.density_fft.from_g(local)
        u.setflags(write=False)
        field.setflags(write=False)
        return BlochDirectionalAction(self, u, field, tuple(terms))
