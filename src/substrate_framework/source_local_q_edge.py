"""Original dual-edge equations with moments and free generators kept separate.

Context publication consumes only n/g/ell; source replay consumes only ket/bra.
The original physical coupling, field actions, Hamiltonians and weights are reused.
"""
import numpy as np


def dual_moments(executor, index, time, state):
    r = executor.r
    n = np.zeros(r.density_shape, complex)
    g, ell = np.zeros(r.ndof, complex), np.zeros(3, complex)
    if executor.dynamics != 'screened-joint':
        return dict(n=n, g=g, ell=ell)
    edge = r.edges[index]
    source, target = r.sectors[edge.source], r.sectors[edge.target]
    with executor.edge_stage(index, time) as (ps, pt, hs, hsa, ht, hta, coupling, axes):
        left, right = 1j*state.ket[index] @ source.d, ps
        left2, right2 = pt, 1j*state.bra[index].conj().T @ target.d.conj().T
        direction = -state.p
        field = -r.sqrt_n*r.polarizations.T @ state.momentum if r.gamma else np.zeros(3)
        paired = getattr(executor.edge_stage, 'paired_moments_and_drive', None)
        if paired is None:
            n1, g1 = coupling.moments(left, right, weight=1.)
            n2, g2 = coupling.moments(left2, right2, weight=1.)
            n, g = n1+n2, g1+g2
        else:
            n, g, unused_drive, unused_adjoint = paired(
                coupling, left, right, left2, right2, direction.reshape(-1, 3), weight=1.)
            del unused_drive, unused_adjoint
        if r.gamma:
            paired_field = getattr(executor.edge_stage, 'paired_field_moments_and_drive', None)
            if paired_field is None:
                ell = (executor._field_moments(axes, left, right)
                       + executor._field_moments(axes, left2, right2))
            else:
                ell, unused_drive, unused_adjoint = paired_field(
                    axes, left, right, left2, right2, field, weight=1.)
                del unused_drive, unused_adjoint
    return dict(n=n, g=g, ell=ell)


def dual_free(executor, index, time, state):
    r = executor.r
    source = r.sectors[r.edges[index].source]
    with executor.edge_stage(index, time) as (ps, pt, hs, hsa, ht, hta, coupling, axes):
        # Exact zeros are omitted algebraically, never through a size threshold.
        if executor.dynamics == 'screened-joint' and np.any(state.p):
            operation = coupling.prepare((-state.p).reshape(-1, 3))
            drive, adjoint = operation.apply(ps), operation.adjoint_apply(pt)
        else:
            drive = np.zeros_like(state.ket[index])
            adjoint = np.zeros(state.bra[index].shape[::-1], dtype=np.complex128)
        if executor.dynamics == 'screened-joint' and r.gamma:
            field = -r.sqrt_n*r.polarizations.T @ state.momentum
            for axis_index, coefficient in enumerate(field):
                if coefficient != 0:
                    axis = axes[axis_index]
                    drive += coefficient*axis.apply(ps)
                    adjoint += coefficient.conjugate()*axis.adjoint_apply(pt)
        ket = 1j*hta(state.ket[index])+source.weight*drive
        bra = (-1j*hs(state.bra[index].conj().T).conj().T
               + source.weight*adjoint.conj().T)
    return ket, bra
