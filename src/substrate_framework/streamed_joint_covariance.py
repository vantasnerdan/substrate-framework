"""Two actual pullbacks contracted at the original product preparation.

No operator is fabricated from incomplete columns. For a physical cotangent,
alpha_s=w A Psi_s, beta_t=w Psi_t† A, so Psi_t†alpha=beta Psi_s.
AB=N[tr(alpha_A†alpha_B Ds)-tr(J_A† Dt J_B Ds)];
BA=N[tr(beta_A† Dt beta_B)-tr(J_A† Dt J_B Ds)].
These identities use rho=Psi D Psi† literally: raw nonorthogonal frames and
singular D require neither orthogonalization nor D/Gram inverses. All empty PW
rows in alpha/beta participate in the first terms. Discrete propagation can
leave a compatibility defect; report it rather than silently project it away.
"""
import numpy as np

from .streamed_joint_response import MappedResponse


def initial_contraction(preparation, left, right):
    r, stage = preparation.response, preparation.initial_stage
    ab, ba = 0j, 0j
    defects = []
    for i, e in enumerate(r.edges):
        ps, pt = stage.frames[e.source], stage.frames[e.target]
        ds, dt = r.sectors[e.source].d, r.sectors[e.target].d
        ja, jb = pt.conj().T @ left.ket[i], pt.conj().T @ right.ket[i]
        defects.append(dict(edge=i, left=float(np.linalg.norm(ja-left.bra[i] @ ps)),
                            right=float(np.linalg.norm(jb-right.bra[i] @ ps))))
        occupied = np.trace(ja.conj().T @ dt @ jb @ ds)
        ab += r.ncells*(np.trace(left.ket[i].conj().T @ right.ket[i] @ ds)-occupied)
        ba += r.ncells*(np.einsum('ij,ij->', left.bra[i].conj(), dt @ right.bra[i])-occupied)
    ai, bi = np.concatenate((left.u, left.p)), np.concatenate((right.u, right.p))
    af, bf = np.concatenate((left.coordinate, left.momentum)), np.concatenate((right.coordinate, right.momentum))
    first = dict(electronic=complex(ab), ionic=complex(np.vdot(ai, preparation.ionic @ bi)),
                 field=complex(np.vdot(af, preparation.field @ bf)))
    second = dict(electronic=complex(ba), ionic=complex(np.vdot(ai, preparation.ionic_minus.T @ bi)),
                  field=complex(np.vdot(af, preparation.field.T @ bf)))
    first['total'], second['total'] = sum(first.values()), sum(second.values())
    return dict(AB=first, BA=second, commutator=first['total']-second['total'],
                symmetrized=(first['total']+second['total'])/2,
                raw_frame_compatibility=defects, preparation_identity=preparation.identity,
                convention='<O_A(q,t) O_B(q,s)^dagger> connected; origins refer to original t0')


class StreamedTwoSidedCovariance:
    def __init__(self, preparation, executor):
        if preparation.response is not executor.r:
            raise ValueError('preparation and propagated q block differ')
        self.preparation, self.executor = preparation, executor

    def pullback(self, probe, times, directory, *, checkpoint=None):
        from pathlib import Path
        directory = Path(directory)
        times = np.asarray(times, float)
        if not len(times) or times[0] != 0 or np.any(np.diff(times) <= 0):
            raise ValueError('pullback grid must start at original t0')
        if hasattr(self.executor.edge_stage, 'clear_interval'):
            self.executor.edge_stage.clear_interval()
        executor = self.executor
        previous = getattr(executor, 'provenance', None)
        provenance = dict(previous or {}, pullback=dict(probe=probe.identity, units=probe.units,
                          terminal_time=float(times[-1]), grid=times.tolist()))
        executor.provenance = provenance
        try:
            if checkpoint is None:
                terminal = MappedResponse(executor.r, directory/'terminal', executor.budget, create=True)
                executor.terminal_probe(times[-1], probe, terminal)
                terminal.publish(time=times[-1], identity=executor.identity, budget=executor.budget,
                                 kind='terminal-cotangent', provenance=provenance)
            else:
                terminal, record = MappedResponse.restore(executor.r, checkpoint, executor.budget, executor.identity)
                if (record['kind'] not in ('dual', 'terminal-cotangent')
                        or record['provenance'] != provenance or record['time'] not in times):
                    raise ValueError('dual restart probe/grid/history differs')
                times = times[times <= record['time']]
            if len(times) == 1:
                return terminal, [checkpoint] if checkpoint else []
            return executor.propagate(terminal, times, directory/'pullback', dual=True)
        finally:
            executor.provenance = previous

    def pair(self, left, right, times, second_times, directory, *, left_checkpoint=None, right_checkpoint=None):
        from pathlib import Path
        directory = Path(directory)
        a, ar = self.pullback(left, times, directory/'left', checkpoint=left_checkpoint)
        b, br = self.pullback(right, second_times, directory/'right', checkpoint=right_checkpoint)
        result = initial_contraction(self.preparation, a, b)
        result.update(time=float(times[-1]), second_time=float(second_times[-1]),
                      dynamics=self.executor.dynamics, propagation_identity=self.executor.identity,
                      left_checkpoints=ar, right_checkpoints=br,
                      whole_spatial_reconstruction=False)
        return result
