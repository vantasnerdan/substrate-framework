"""Budgeted matrix-function action of one actual physical full-PW H0.

The endpoint/basis is owned by the caller. Independent rectangular legs and
original mean columns use the same arbitrary-full-row operation. No eigenband
projection or assumption of exact floating Hermiticity is used.
"""
from pathlib import Path
from uuid import uuid4

import numpy as np

from joint_column_mean import ColumnLimits
from moving_exponential_panels import bounded_panel_columns
from run_column_krylov_material import operator_norm
from run_block_krylov_material import ambient_growth
from substrate_framework.material_frozen_action import FrozenMaterialAction


class PhysicalH0Exponential:
    def __init__(self, h0, *, budget, limits, identity, receipt_directory,
                 fft_block_size=8):
        if not isinstance(limits, ColumnLimits) or not identity:
            raise ValueError('explicit full-column limits and physical operator identity required')
        self.budget, self.limits, self.identity = budget, limits, identity
        self.frozen = FrozenMaterialAction(h0, fft_block_size=fft_block_size)
        self.norm = operator_norm(self.frozen)[0]
        self.growth = ambient_growth(self.frozen)[0]
        self.directory = Path(receipt_directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.actions = self.columns = 0

    def propagate(self, vectors, interval, *, adjoint=False):
        self.budget.check()
        if not np.isfinite(interval):
            raise ValueError('finite physical interval required')
        if interval == 0 or not np.any(vectors):
            return vectors
        action = self.frozen.adjoint_apply if adjoint else self.frozen.apply

        def counted(value):
            self.budget.action(value)
            self.actions += 1
            self.columns += 1 if value.ndim == 1 else value.shape[1]
            return action(value)

        evolved, receipt = bounded_panel_columns(counted, vectors, interval,
            limits=self.limits, norm_bound=self.norm, growth_bound=self.growth,
            identity=self.identity+':'+str(adjoint)+':'+repr(float(interval)))
        self.budget.json(self.directory/(uuid4().hex+'.json'),
            dict(kind='physical-H0-exponential', operator_identity=self.identity,
                 interval_ry_inverse=float(interval), adjoint=adjoint, **receipt))
        return evolved
