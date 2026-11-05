"""Same-source residual substitution, including preparation and all cross orders.

A30's initial sigma dual ALREADY contains represented feedback. A31/A38
'retained_vertex' means P O P propagated by the SAME frozen full-bath H;
it is NOT the original moving active response and cannot replace prepdual.
We add only full-minus-retained force/current/density residual operators to
that original prepdual. This declares an external-residual driven model,
not the selfconsistent A32 hybrid mean or a whole-E covariance correction.
"""
import numpy as np


def compose_same_source(preparation, preparation_adjoint, ionic_residual,
                        ionic_adjoint, field_residual, field_adjoint,
                        density_residual, density_adjoint):
    matrices = dict(preparation=preparation, ionic=ionic_residual,
                    field=field_residual, density=density_residual)
    adjoints = dict(preparation=preparation_adjoint, ionic=ionic_adjoint,
                    field=field_adjoint, density=density_adjoint)
    matrices['matched'] = sum(matrices.values())
    adjoints['matched'] = sum(adjoints.values())
    if not all(np.isfinite(a).all() for a in (*matrices.values(), *adjoints.values())):
        raise ValueError('nonfinite same-source operator')
    return matrices, adjoints


def substitution_readout(moments, original_ordered, represented_ordered):
    """Replace source0's old preparation contribution; NEVER add a variance.

    The original total must be from the SAME backward dual, baseline, terminal
    probe and state used to embed preparation. It is not an A32 hybrid total.
    All16 ordered cross terms of the four component operators are retained.
    represented_ordered is the actual source0 active Wick term in that total,
    not an inferred P O P variance. Its full-frame CAR compression defect is
    reported separately; numerical state defects are not silently equated.
    """
    values = {(r['left'], r['right']): r['value'] for r in moments['covariance']}
    parts = ('preparation', 'ionic', 'field', 'density')
    required = {(a, b) for a in (*parts, 'matched') for b in (*parts, 'matched')}
    if set(values) != required:
        raise ValueError('complete25 ordered same-source pairs required')
    old = represented_ordered
    new = values['matched', 'matched']
    return dict(old_source_preparation=old, matched_source=new,
        same_source_substitution=new-old,
        conditional_matched_active_total=original_ordered-old+new,
        raw_preparation_compression_difference=values['preparation', 'preparation']-old,
        all_cross_orders=[dict(left=a, right=b, value=values[a, b])
                          for a in parts for b in parts if a != b],
        raw_all16_reconstruction=new-sum(values[a, b] for a in parts for b in parts),
        forbidden_independent_variance_addition_error=
            new-old-sum(values[a, a] for a in parts[1:]),
        scope='external frozen-bath residual driven original moving active model; not A32 hybrid or full physical E')
