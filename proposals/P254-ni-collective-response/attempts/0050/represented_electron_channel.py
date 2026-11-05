"""Eliminate ONLY A29 active electron fluctuations with prescribed mechanics.

Retains the moving mean and density->K->delta_h tangent. This is NOT the
undressed A41 bubble. Noise and response use the SAME original active sector.
The public maps also permit prescribed mechanical forcing without elimination.
"""
import numpy as np
from scipy.integrate import solve_ivp
from joint_gamma_adapter import vertices
from reduced_correlation_transpose import commutator_tangent, commutator_transpose, wick_edge_pair


class RepresentedElectronChannel:
    def __init__(self, coefficients, baseline, *, rtol=2e-8, atol=2e-10):
        self.c, self.baseline = coefficients, baseline
        self.block = coefficients.q(0)
        self.forward, self.reverse = vertices(coefficients, self.block)
        self.shape = (coefficients.ns, coefficients.m, coefficients.m)
        self.rtol, self.atol = rtol, atol
        self.rhs_evaluations = 0

    def forcing(self, time, coordinates):
        rho, h = self.baseline(time)
        dh = np.einsum('siab,ip->sabp', self.forward, coordinates, optimize=True)
        return commutator_tangent(h, h, rho, rho, np.zeros_like(dh), dh)

    def gradient(self, sigma):
        return np.einsum('siab,sabp->ip', self.reverse, sigma, optimize=True)

    def tangent(self, time, sigma, coordinates=None):
        c, b = self.c, self.block
        rho, h = self.baseline(time)
        density = np.einsum('srab,sabp->rp', b.density, sigma, optimize=True)/c.sqrt_n
        potential = b.kernel @ density
        dh = np.einsum('sabr,rp->sabp', b.potential, potential, optimize=True)/c.sqrt_n
        if coordinates is not None:
            dh += np.einsum('siab,ip->sabp', self.forward, coordinates, optimize=True)
        return commutator_tangent(h, h, rho, rho, sigma, dh)

    def transpose(self, time, dual):
        c, b = self.c, self.block
        rho, h = self.baseline(time)
        result, dh_dual = commutator_transpose(h, h, rho, rho, dual)
        feature = np.einsum('sabr,sabp->rp', b.potential, dh_dual, optimize=True)/c.sqrt_n
        result += np.einsum('srab,rp->sabp', b.density, b.kernel.T@feature,
                            optimize=True)/c.sqrt_n
        return result

    def response_from_dual(self, dual, time):
        rho, h = self.baseline(time)
        _, dh_dual = commutator_transpose(h, h, rho, rho, dual)
        # Output index first, input coordinate second. Bilinear, not vdot.
        return np.einsum('sabp,siab->pi', dh_dual, self.forward, optimize=True)

    def pullback(self, terminal, sample_times, coordinates=None):
        """C U_e(terminal,s), sampled only at the explicitly requested times.

        solve_ivp evaluates the original saved degree7 moving baseline; no
        interpolation between the four contour samples or inferred U is used.
        """
        times = sorted({float(t) for t in sample_times if t <= terminal}, reverse=True)
        if not times or times[-1] != 0 or times[0] != terminal:
            raise ValueError('include zero and terminal in pullback sample times')
        indices = np.arange(self.c.d+self.c.f) if coordinates is None else np.asarray(coordinates)
        if indices.ndim != 1 or not len(indices) or np.any(indices<0) or np.any(indices>=self.c.d+self.c.f):
            raise ValueError('nonempty original joint coordinate indices required')
        initial = self.reverse[:,indices].transpose(0, 2, 3, 1).copy()
        if terminal == 0:
            return {0.: initial}, dict(nfev=0, terminal=0., success=True,coordinates=indices.tolist())
        def rhs(t, flat):
            self.rhs_evaluations += 1
            return -self.transpose(t, flat.reshape(initial.shape)).ravel()
        result = solve_ivp(rhs, (terminal, 0.), initial.ravel(), t_eval=times,
                           method='DOP853', rtol=self.rtol, atol=self.atol)
        if not result.success:
            raise RuntimeError(result.message)
        values = {float(t): result.y[:, j].reshape(initial.shape)
                  for j, t in enumerate(result.t)}
        return values, dict(nfev=int(result.nfev), terminal=float(terminal), success=True,coordinates=indices.tolist())

    def ordered_noise(self, left_at_zero, right_at_zero):
        c = self.c
        return wick_edge_pair(left_at_zero, right_at_zero[c.target[0]],
                              c.rho0, c.rho0[c.target[0]])


def polarization_difference(new, represented):
    """Pi-only new-minus-represented in identical units, NEVER subtract Sigma.

    new and represented contain pi_retarded/pi_greater/pi_lesser. Represented
    retarded is the actual driven tangent response, not inferred from its noise.
    """
    names = ('pi_retarded', 'pi_greater', 'pi_lesser')
    result = {}
    for name in names:
        if np.shape(new[name]) != np.shape(represented[name]):
            raise ValueError('same coordinate/time/component axes required')
        result[name] = new[name] - represented[name]
    return result


CUTOVER = {
    'remove_response': 'A29 _tangent sigma -> ionic_reverse/sqrtN + field_reverse_modes -> mechanical.j gradient, for this active sector ONLY',
    'remove_noise': 'The SAME active sector Wick propagation in initial_ordered; keep original mechanical w0 and electronic rho0, do not add it again as independent force/current noise',
    'retain': 'moving mean, rho-rho0 forces/tadpoles, passive initial totals, complete mechanical generators/j/w0, mixed/contact/minimal-coupling/Maxwell blocks',
    'alternative': 'retain explicit tangent and insert new-minus-represented Pi in the identical joint coordinate components; Sigma is not a represented mechanical channel',
    'unresolved': 'No ab-initio PBE/XC correlation-residual double-count prescription is supplied',
}
