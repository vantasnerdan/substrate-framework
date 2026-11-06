"""Full-PW electronic time generator in Ry time (hbar=1).

For exp(-i*z*tau), plus amplitudes have denominator z-(eps_plus-eps_s).
The stored minus amplitudes are conjugated ket coefficients: their denominator
is -z-(eps_minus-eps_s). Thus i*a_dot=gap_plus*a+drive_plus and
 i*b_dot=-gap_minus*b-conjugate(drive_minus[v*,u*]). Both state and drives
enter complex-linearly. There is no numerical eta, screening, mass, fitted
counterterm, or population-relaxation term in this generator.

State packing is channel order, plus then minus for each channel, each a
(target_eigenstate, source_eigenstate) array flattened in C order. All supplied
nonzero-FD sources and complete target spectra are retained without a cutoff.
Models and arrays must remain immutable while this object is in use.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .pw_response_data import PWResponseModel, PWResponseValue
from .pw_transition_space import PWTransitionSpace, finite_array


@dataclass(frozen=True)
class PWTransitionBlock:
    """Public packing descriptor; state[start:stop].reshape(shape) is a view."""

    channel_index: int
    branch: str
    start: int
    stop: int
    shape: tuple[int, int]


@dataclass(frozen=True)
class _Branch:
    block: PWTransitionBlock
    coordinates: tuple[NDArray, ...]
    vectors: NDArray
    adjoint: NDArray
    gap: NDArray
    ionic_drive: NDArray
    ionic_observable: NDArray


class PWElectronicDynamics:
    """Reusable complex-linear electronic generator with physical observables.

    ``state_layout`` documents stable flat packing. ``derivative`` returns
    d/dtau, not i*d/dtau; observables return physical spin density and positive
    ionic gradient delta<g_-q>, not its negative mechanical force. Occupation
    weights belong to observables, not amplitude evolution.
    """

    def __init__(self, model: PWResponseModel):
        self.model = model
        self._space = PWTransitionSpace(model)
        branches = []
        offset = 0
        for index, (channel, maps) in enumerate(zip(model.channels, self._space.maps)):
            for name, coordinates in zip(("plus", "minus"), maps):
                vectors = getattr(channel, f"{name}_vectors")
                adjoint = vectors.conj().T
                gap = getattr(channel, f"{name}_energies")[:, None] - channel.source_energies[None, :]
                shape = gap.shape
                stop = offset + gap.size
                block = PWTransitionBlock(index, name, offset, stop, shape)
                spectral_g = adjoint @ getattr(channel, f"g_{name}")
                drive = spectral_g if name == "plus" else spectral_g.conj()
                observable = (spectral_g.conj() if name == "plus" else spectral_g)
                observable = observable * channel.weighted_occupations[None, None, :]
                branches.append(_Branch(block, coordinates, vectors, adjoint,
                                        gap, drive, observable))
                offset = stop
        self._branches = tuple(branches)
        self.state_layout = tuple(branch.block for branch in branches)
        self.state_size = offset
        self._source_conjugates = tuple(c.source_grid.conj() for c in model.channels)

    def zero_state(self) -> NDArray[np.complex128]:
        return np.zeros(self.state_size, dtype=np.complex128)

    def derivative(self, state, local_potential, displacement) -> NDArray[np.complex128]:
        """Evaluate the undamped isolated electronic EOM at arbitrary complex state."""
        state = finite_array(state, (self.state_size,), "electronic state")
        potential = finite_array(local_potential, (self.model.fft.nnr, 2), "local potential")
        u = finite_array(displacement, (6,), "displacement")
        result = np.empty(self.state_size, dtype=np.complex128)
        for branch in self._branches:
            block = branch.block
            channel = self.model.channels[block.channel_index]
            amplitude = state[block.start:block.stop].reshape(block.shape)
            v = potential[:, channel.spin]
            if block.branch == "minus":
                v = v.conj()
            local = branch.adjoint @ self._space.to_pw(
                v[:, None] * channel.source_grid, branch.coordinates)
            if block.branch == "minus":
                local = local.conj()
            drive = np.einsum("atn,a->tn", branch.ionic_drive, u)
            drive += local
            sign = -1j if block.branch == "plus" else 1j
            result[block.start:block.stop] = (sign * (branch.gap * amplitude + drive)).reshape(-1)
        return result

    def observables(self, state) -> PWResponseValue:
        """Contract the state with full physical FFT/volume/spin/k normalization."""
        state = finite_array(state, (self.state_size,), "electronic state")
        fft = self.model.fft
        density = np.zeros((fft.nnr, 2), dtype=np.complex128)
        gradient = np.zeros(6, dtype=np.complex128)
        for branch in self._branches:
            block = branch.block
            channel = self.model.channels[block.channel_index]
            amplitude = state[block.start:block.stop].reshape(block.shape)
            if block.branch == "plus":
                target_grid = self._space.to_grid(branch.vectors @ amplitude, branch.coordinates)
                source_grid = self._source_conjugates[block.channel_index]
            else:
                target_grid = self._space.to_grid(
                    branch.vectors @ amplitude.conj(), branch.coordinates).conj()
                source_grid = channel.source_grid
            density[:, channel.spin] += np.einsum(
                "rn,rn,n->r", source_grid, target_grid, channel.weighted_occupations) / fft.omega
            gradient += np.einsum("atn,tn->a", branch.ionic_observable, amplitude)
        return PWResponseValue(density, gradient)
