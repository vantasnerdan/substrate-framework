"""Physical boundary readouts, separate from RK stages and covariance means.

Connected parts are Cov(A(t), B(0)) (or its reverse order), NOT equal-time
Cov(A(t), B(t)). No unspecified source part is silently promoted to zero.
"""
from __future__ import annotations

import numpy as np
from .prepared_covariance_source import PreparedProbe, ordered_readout, probe_mean
from .prepared_reference import array, frozen


def regional_number_probe(model, mask, *, identity):
    """Per-cell number in an explicit full-grid real region (same for both spins).

    The physical quadrature is supplied by each coupling's potential action.
    Collective connected readouts retain the source's existing sqrt(N) scaling.
    """
    region = frozen(array(mask, (model.reference.nbar.shape[0],), 'region mask', real=True))
    if np.any(region < 0) or np.any(region > 1):
        raise ValueError('regional mask must lie in [0,1]')
    actions = tuple((lambda vectors, c=s.coupling: c.potential_action(region, vectors)) for s in model.sectors)
    zero = np.zeros(model.reference.ndof)
    held = model.measure*np.sum(region[:, None]*model.reference.held_density)
    return PreparedProbe(actions, actions, zero, zero, complex(held), identity, 'electrons/cell')


def boundary_diagnostics(model, state):
    """Raw energy/trace/Gram measurements; no additional mean/tangent RHS call.

    Rank-sized Gram arrays are emitted only on explicit request to this helper.
    No clipping, Hermitization or conservation guarantee is made.
    """
    r = model.reference
    density = np.zeros(r.nbar.shape, complex)
    gradient = np.zeros(r.ndof, complex)
    energy = 0j
    grams = []
    for sector, psi in zip(model.sectors, state.frames):
        n, g = sector.coupling.moments(psi, psi, sector.d)
        density[:, sector.density_component] += sector.weight*n
        gradient += sector.weight*g
        energy += sector.weight*np.trace(sector.d @ (psi.conj().T @ sector.h0(psi)))
        grams.append(psi.conj().T @ psi)
    delta = density + r.held_density - r.nbar
    field = r.kernel(delta)
    energy += (model.measure*np.sum(delta*field)/2
               + state.u @ (gradient+r.held_gradient-r.ionic_force-r.support_force)
               + state.u @ r.curvature @ state.u/2 + np.sum(state.p*state.p/r.masses)/2)
    return dict(energy=complex(energy), traces=model.traces(state), grams=tuple(grams))


def physical_readout(model, boundary, probes):
    """Compact means and separately labeled initial-origin ordered connections."""
    return dict(time=boundary.time, completed_steps=boundary.completed_steps,
                probes={name: dict(identity=probe.identity, units=probe.units,
                                   mean=probe_mean(model, boundary.mean, probe),
                                   connected=[dict(right_identity=s.right_identity,
                                                   preparation_authority=s.preparation_authority,
                                                   **ordered_readout(model, boundary.mean, s, probe))
                                              for s in boundary.sources])
                        for name, probe in probes.items()},
                connection_scope='A(t) against B(initial); not both observables evolved')


def stream_readouts(evolution, probes, *, dt, steps):
    """Yield completed-boundary physical records without retaining states/history."""
    for boundary in evolution.stream(dt=dt, steps=steps):
        yield physical_readout(evolution.model, boundary, probes)
