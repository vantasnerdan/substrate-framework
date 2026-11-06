"""A29 coupled Gaussian tangent at literal A78 mean/rho nodes.

Only the coherent Gamma curvature is added. Original electronic density/K,
ionic/field vertices, contacts and product preparation remain owned by A29.
"""
import numpy as np
from reduced_electronic_correlations import ActiveBandCorrelations
from produce_mean_feedback import _node_value
import mean_control_blocks as blocks




class HeldMeanBaseline:
    def __init__(self, functional, times, rho, mean):
        self.functional, self.times = functional, times
        self.rho, self.mean = rho, mean

    def __call__(self, time):
        rho = _node_value(time, self.times, self.rho)
        mean = _node_value(time, self.times, self.mean)
        return rho, self.functional.hamiltonian(rho, mean)


class CoupledMeanControl(ActiveBandCorrelations):
    """Original full coupled tangent, with the common Gamma Phi000*mean block."""
    def __init__(self, coefficients, functional, times, rho, mean, **settings):
        super().__init__(coefficients,
            HeldMeanBaseline(functional, times, rho, mean), **settings)
        self.times = times
        self.nx = len(coefficients.mechanical.generators[0])
        self.j_columns = coefficients.mechanical.j[0][:, :coefficients.d]
        # Curvature is linear in mean: interpolate the precontracted shifts,
        # never subtract two large Hessians or repeat Phi000 contractions/RHS.
        self.curvature_nodes = np.einsum('ijk,tk->tij', functional.phi,
            mean[:, :coefficients.d], optimize=True)

    def _curvature(self, q, time):
        if q != 0:
            raise ValueError('this projected control owns Gamma only')
        return _node_value(time, self.times, self.curvature_nodes)

    def _tangent(self, q, time, value, block):
        curvature = self._curvature(q, time)
        x, _, ndim = self._split(q, value)
        result = super()._tangent(q, time, value, block)
        increment = blocks.gamma_tangent_increment(self.j_columns, curvature,
            x[:self.c.d])
        result[:self.nx] += increment[:, 0] if ndim == 1 else increment
        return result

    def _transpose(self, q, time, value, block):
        curvature = self._curvature(q, time)
        dual, _, ndim = self._split(q, value)
        result = super()._transpose(q, time, value, block)
        increment = blocks.gamma_transpose_increment(self.j_columns, curvature, dual)
        result[:self.c.d] += increment[:, 0] if ndim == 1 else increment
        return result
