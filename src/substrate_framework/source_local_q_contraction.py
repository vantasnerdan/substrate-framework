"""Original product contraction, partitioned into independent electronic sources.

Each worker reconstructs both full-PW cotangents for one original edge. Only
its two ordered scalars and compatibility norms are retained. Reduction follows
original edge order; canonical preparation terms are added once, not per source.
The equations are those of streamed_joint_covariance.initial_contraction.
"""
from pathlib import Path

import numpy as np

from .partitioned_q_storage import release_mapped_pages
from .partitioned_q_work import QWork, arrays
from .streamed_material_state import file_binding


def _electronic_source(preparation, left, right, index):
    r, stage = preparation.response, preparation.initial_stage
    edge = r.edges[index]
    ps, pt = stage.frames[edge.source], stage.frames[edge.target]
    ds, dt = r.sectors[edge.source].d, r.sectors[edge.target].d
    try:
        ak, ab = left.edge(index)
        bk, bb = right.edge(index)
        ja, jb = pt.conj().T @ ak, pt.conj().T @ bk
        occupied = np.trace(ja.conj().T @ dt @ jb @ ds)
        return dict(
            ordered=np.asarray([
                r.ncells*(np.trace(ak.conj().T @ bk @ ds)-occupied),
                r.ncells*(np.einsum('ij,ij->', ab.conj(), dt @ bb)-occupied)], complex),
            compatibility=np.asarray([np.linalg.norm(ja-ab @ ps),
                                      np.linalg.norm(jb-bb @ ps)], float))
    finally:
        left.release()
        right.release()


def source_contraction(preparation, left, right, executor, directory):
    directory = Path(directory)
    work = QWork(directory, executor.budget, executor.edge_work.cpus, executor.identity)
    request = dict(preparation=preparation.identity, left=left.binding, right=right.binding)
    records = work.execute('source-local-initial-contraction', request,
        len(executor.r.edges), lambda index: _electronic_source(preparation, left, right, index))
    ab, ba, defects = 0j, 0j, []
    for index, record in enumerate(records):
        values = arrays(record)
        ab += complex(values['ordered'][0])
        ba += complex(values['ordered'][1])
        defects.append(dict(edge=index, left=float(values['compatibility'][0]),
                            right=float(values['compatibility'][1])))
        for value in values.values():
            release_mapped_pages(value)
    ai, bi = np.concatenate((left.u, left.p)), np.concatenate((right.u, right.p))
    af, bf = np.concatenate((left.coordinate, left.momentum)), np.concatenate((right.coordinate, right.momentum))
    first = dict(electronic=ab, ionic=complex(np.vdot(ai, preparation.ionic @ bi)),
                 field=complex(np.vdot(af, preparation.field @ bf)))
    second = dict(electronic=ba, ionic=complex(np.vdot(ai, preparation.ionic_minus.T @ bi)),
                  field=complex(np.vdot(af, preparation.field.T @ bf)))
    first['total'], second['total'] = sum(first.values()), sum(second.values())
    root = work.last[1]
    return dict(AB=first, BA=second, commutator=first['total']-second['total'],
        symmetrized=(first['total']+second['total'])/2,
        raw_frame_compatibility=defects, preparation_identity=preparation.identity,
        convention='<O_A(q,t) O_B(q,s)^dagger> connected; origins refer to original t0',
        source_contraction=dict(request=file_binding(root/'request.json'), count=len(records),
            receipts=[file_binding(root/('item-%04d' % index)/'complete.json')
                      for index in range(len(records))],
            reduction='original source order; canonical terms added once'))
