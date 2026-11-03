"""Finite-amplitude baseline-plus-delta evolution of the shifted active scalar.

Times are inverse Ry. Densities are global canonical densities, with explicit
spins and no sector weights. Complex residuals are never discarded or repaired.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from reduced_coefficients import ReducedCoefficients


@dataclass(frozen=True)
class BaselineState:
    rho: np.ndarray
    x: np.ndarray


@dataclass(frozen=True)
class DeltaState:
    rho: np.ndarray
    x: tuple


@dataclass(frozen=True)
class EnergyDifference:
    mechanical: complex
    band: complex
    interaction: complex
    screening: complex

    @property
    def total(self):
        return self.mechanical+self.band+self.interaction+self.screening

    def as_dict(self):
        return dict(mechanical=self.mechanical, band=self.band,
                    interaction=self.interaction, screening=self.screening,
                    total=self.total)


class BaselineSolution:
    def __init__(self, model, result):
        self.model, self.result = model, result

    def state(self, time):
        if not self.result.sol.t_min <= time <= self.result.sol.t_max:
            raise ValueError('baseline interpolation outside integrated interval')
        return self.model.unpack_baseline(self.result.sol(time))

    def __call__(self, time):
        state = self.state(time)
        return state.rho, self.model.baseline_hamiltonian(state)


class FiniteDoseSolution:
    """Requested finite-dose states, with optional explicitly retained interpolation."""
    def __init__(self, model, baseline, result):
        self.model, self.baseline, self.result = model, baseline, result

    def state(self, time):
        retained = np.flatnonzero(self.result.t == time)
        if len(retained):
            return self.model.unpack_delta(self.result.y[:,retained[0]])
        if self.result.sol is None:
            raise ValueError('finite-dose time was not retained; request it in t_eval '
                             'or explicitly enable dense_output=True')
        if not self.result.sol.t_min <= time <= self.result.sol.t_max:
            raise ValueError('delta interpolation outside integrated interval')
        return self.model.unpack_delta(self.result.sol(time))

    def energy(self, time):
        return self.model.energy_difference(self.baseline.state(time), self.state(time))


class ActiveBandDynamics:
    def __init__(self, coefficients: ReducedCoefficients):
        self.c = coefficients
        self.length = coefficients.nk*coefficients.m
        self.rho_shape = (coefficients.nspin,self.length,self.length)
        self.rho_size = int(np.prod(self.rho_shape))
        self.baseline_rho_size = coefficients.ns*coefficients.m**2
        self.x_offsets = np.cumsum((0,)+coefficients.mechanical_sizes)
        self.delta_size = self.rho_size+int(self.x_offsets[-1])

    def initial_baseline(self):
        return BaselineState(self.c.rho0.copy(), np.zeros(self.c.mechanical_sizes[0],complex))

    def initial_delta(self, coordinate=0, cell=(0,0,0), amplitude=1.):
        """Unit canonical momentum impulse by default: Ni1/x, original cell zero.

        This is delta p=1, NOT the clamped helper's one-trap-quantum kick.
        All original q sectors receive the orthonormal Fourier impulse.
        """
        c = self.c
        cell = np.asarray(cell)
        if (cell.shape != (3,) or not np.issubdtype(cell.dtype,np.integer)
                or np.any(cell < 0) or np.any(cell >= c.mesh)):
            raise ValueError('original integer cell inside periodic mesh required')
        if not 0 <= coordinate < c.d or not np.isfinite(amplitude):
            raise ValueError('finite impulse and original Cartesian coordinate required')
        phase = np.exp(-2j*np.pi*(c.q_indices/np.asarray(c.mesh))@cell)
        x = tuple(np.zeros(size,complex) for size in c.mechanical_sizes)
        for q in range(c.nq):
            x[q][c.d+coordinate] = amplitude*phase[q]/c.sqrt_n
        return DeltaState(np.zeros(self.rho_shape,complex),x)

    def pack_baseline(self, state):
        c = self.c
        if state.rho.shape != c.rho0.shape or state.x.shape != (c.mechanical_sizes[0],):
            raise ValueError('baseline shape mismatch')
        return np.concatenate((state.rho.ravel(),state.x)).astype(complex,copy=False)

    def unpack_baseline(self, vector):
        vector = np.asarray(vector)
        c = self.c
        if vector.shape != (self.baseline_rho_size+c.mechanical_sizes[0],):
            raise ValueError('baseline packed shape mismatch')
        return BaselineState(vector[:self.baseline_rho_size].reshape(c.rho0.shape),
                             vector[self.baseline_rho_size:])

    def pack_delta(self, state):
        if state.rho.shape != self.rho_shape or len(state.x) != self.c.nq:
            raise ValueError('finite delta shape mismatch')
        for x,size in zip(state.x,self.c.mechanical_sizes):
            if x.shape != (size,):
                raise ValueError('mechanical delta block shape mismatch')
        return np.concatenate((state.rho.ravel(),*state.x)).astype(complex,copy=False)

    def unpack_delta(self, vector):
        vector = np.asarray(vector)
        if vector.shape != (self.delta_size,):
            raise ValueError('finite delta packed shape mismatch')
        mechanical = vector[self.rho_size:]
        return DeltaState(vector[:self.rho_size].reshape(self.rho_shape),
                          tuple(mechanical[a:b] for a,b in zip(self.x_offsets[:-1],self.x_offsets[1:])))

    def q_blocks(self, rho, q):
        """Global rho_(k+q,k), in explicit source-slot order."""
        c = self.c
        blocks = rho.reshape(c.nspin,c.nk,c.m,c.nk,c.m)
        spin = np.arange(c.ns)//c.nk
        return blocks[spin,c.slot_k[c.target[q]],:,c.slot_k,:]

    def density(self, data, blocks):
        return np.einsum('srab,sab->r',data.density,blocks)/self.c.sqrt_n

    def _gradient(self, data, blocks, gamma=False):
        c = self.c
        gradient = np.zeros(2*c.d+(2*c.f if gamma else 0),complex)
        gradient[:c.d] = np.einsum('siab,sab->i',data.ionic_reverse,blocks)/c.sqrt_n
        if gamma:
            gradient[2*c.d:2*c.d+c.f] = np.einsum('sfab,sab->f',c.field_reverse_modes,blocks)
        return gradient

    def baseline_hamiltonian(self, state, data=None):
        c = self.c
        data = c.q(0) if data is None else data
        density = self.density(data,state.rho-c.rho0)
        out = c.h0 + np.einsum('siab,i->sab',data.ionic,state.x[:c.d])/c.sqrt_n
        out = out + np.einsum('sabr,r->sab',data.potential,data.kernel@density)/c.sqrt_n
        return out + np.einsum('sfab,f->sab',c.field_modes,state.x[2*c.d:2*c.d+c.f])

    def baseline_rhs(self, time, vector):
        del time
        c = self.c
        state = self.unpack_baseline(vector)
        data = c.q(0)
        h = self.baseline_hamiltonian(state,data)
        rho_dot = -1j*(h@state.rho-state.rho@h)
        gradient = self._gradient(data,state.rho-c.rho0,True)
        x_dot = c.mechanical.generators[0]@state.x+c.mechanical.j[0]@(c.mechanical.linear+gradient)
        return self.pack_baseline(BaselineState(rho_dot,x_dot))

    def delta_hamiltonian(self, state):
        """Assemble the finite global delta H; no equilibrium or tangent shortcut."""
        return self._delta_terms(state, feedback=False)[0]

    def _delta_terms(self, state, *, feedback):
        c = self.c
        h = np.zeros(self.rho_shape,complex)
        view = h.reshape(c.nspin,c.nk,c.m,c.nk,c.m)
        spin = np.arange(c.ns)//c.nk
        x_dot = []
        for q in range(c.nq):
            data = c.q(q)
            blocks = self.q_blocks(state.rho,q)
            density = self.density(data,blocks)
            block = (np.einsum('siab,i->sab',data.ionic,state.x[q][:c.d])
                     + np.einsum('sabr,r->sab',data.potential,data.kernel@density))/c.sqrt_n
            if q == 0:
                block += np.einsum('sfab,f->sab',c.field_modes,state.x[0][2*c.d:2*c.d+c.f])
            view[spin,c.slot_k[c.target[q]],:,c.slot_k,:] = block
            if feedback:
                gradient = self._gradient(data,blocks,q == 0)
                x_dot.append(c.mechanical.generators[q]@state.x[q]+c.mechanical.j[q]@gradient)
        return h,tuple(x_dot)

    def delta_rhs_state(self, baseline, state, baseline_h=None):
        c = self.c
        hb = self.baseline_hamiltonian(baseline) if baseline_h is None else baseline_h
        dh,x_dot = self._delta_terms(state,feedback=True)
        rho_dot = np.empty_like(state.rho,dtype=complex)
        for spin, slots in enumerate(c.slots_by_spin):
            h = hb[slots]
            rho_b = baseline.rho[slots]
            delta = state.rho[spin].reshape(c.nk,c.m,c.nk,c.m)
            perturbation = dh[spin].reshape(c.nk,c.m,c.nk,c.m)
            # Block-diagonal homogeneous products; never assemble a dense H_b.
            comm = np.einsum('kab,kblc->kalc',h,delta)-np.einsum('kalb,lbc->kalc',delta,h)
            comm += np.einsum('kalb,lbc->kalc',perturbation,rho_b)-np.einsum('kab,kblc->kalc',rho_b,perturbation)
            # The indispensable finite-amplitude term, absent in a tangent model.
            rho_dot[spin] = -1j*(comm.reshape(self.length,self.length)
                                +dh[spin]@state.rho[spin]-state.rho[spin]@dh[spin])
        return DeltaState(rho_dot,x_dot)

    def delta_rhs(self, time, vector, baseline):
        state_b = baseline.state(time)
        return self.pack_delta(self.delta_rhs_state(state_b,self.unpack_delta(vector)))

    def canonical_velocity(self, state):
        """q ionic velocities, including the one global field-charge term."""
        c = self.c
        return tuple((c.mechanical.generators[q]@state.x[q])[:c.d] for q in range(c.nq))

    def energy_difference(self, baseline, state):
        """Exact E(b+delta)-E(b), without subtracting two large total scalars."""
        c = self.c
        xb = tuple(baseline.x if q == 0 else np.zeros(size,complex)
                   for q,size in enumerate(c.mechanical_sizes))
        mechanical = c.mechanical.linear@state.x[0]
        for q,h in enumerate(c.mechanical.h):
            minus = c.minus[q]
            mechanical += (xb[minus]@h@state.x[q]+state.x[minus]@h@xb[q]
                           +state.x[minus]@h@state.x[q])/2
        diagonal = self.q_blocks(state.rho,0)
        band = np.einsum('sab,sba->',c.h0,diagonal)
        interaction = 0j
        screening = 0j
        blocks = state.rho.reshape(c.nspin,c.nk,c.m,c.nk,c.m)
        spin = np.arange(c.ns)//c.nk
        # Only small feature vectors are retained, not a second coefficient bank.
        delta_density = []
        screened_delta = []
        partners = []
        for q in range(c.nq):
            data = c.q(q)
            delta_density.append(self.density(data,self.q_blocks(state.rho,q)))
            screened_delta.append(data.kernel@delta_density[q])
            partners.append(data.feature_partner)
            potential = np.einsum('siab,i->sab',data.ionic,state.x[q][:c.d])/c.sqrt_n
            if q == 0:
                potential += np.einsum('sfab,f->sab',c.field_modes,state.x[0][2*c.d:2*c.d+c.f])
                base_potential = np.einsum('siab,i->sab',data.ionic,baseline.x[:c.d])/c.sqrt_n
                base_potential += np.einsum('sfab,f->sab',c.field_modes,baseline.x[2*c.d:2*c.d+c.f])
                interaction += np.einsum('sab,sba->',base_potential,diagonal)
                interaction += np.einsum('sab,sba->',potential,baseline.rho-c.rho0)
                base_density = self.density(data,baseline.rho-c.rho0)
                screening += (base_density[data.feature_partner]@screened_delta[0]
                              +delta_density[0][data.feature_partner]@data.kernel@base_density)/2
            # Tr V_(target,source) delta-rho_(source,target), never a dagger.
            opposite = blocks[spin,c.slot_k,:,c.slot_k[c.target[q]],:]
            interaction += np.einsum('sab,sba->',potential,opposite)
        for q in range(c.nq):
            screening += delta_density[c.minus[q]][partners[q]]@screened_delta[q]/2
        return EnergyDifference(mechanical,band,interaction,screening)

    def integrate_baseline(self, interval, initial=None, *, rtol=2e-8, atol=2e-10,
                           method='DOP853', max_step=np.inf, t_eval=None):
        initial = self.initial_baseline() if initial is None else initial
        result = solve_ivp(self.baseline_rhs,interval,self.pack_baseline(initial),
                           method=method,rtol=rtol,atol=atol,max_step=max_step,
                           t_eval=t_eval,dense_output=True)
        if not result.success:
            raise RuntimeError(f'baseline integration failed: {result.message}')
        return BaselineSolution(self,result)

    def integrate_finite(self, interval, initial=None, *, baseline=None,
                         baseline_initial=None, rtol=2e-8, atol=2e-10,
                         method='DOP853', max_step=np.inf, t_eval=None,
                         dense_output=False):
        """Integrate while retaining only requested global-density states.

        t_eval=None retains the final time only, not every adaptive step.
        Production consumers should request a bounded physical output grid;
        state(t) and energy(t) then accept exact retained times only.
        dense_output=True is an explicit small-model/debug option: its SciPy
        interpolants retain several full density arrays for every solver step.
        Baseline interpolation remains dense independently for correlations.
        """
        output_times = np.asarray([interval[1]] if t_eval is None else t_eval,dtype=float)
        if output_times.ndim != 1 or not len(output_times):
            raise ValueError('a nonempty one-dimensional t_eval output grid is required')
        if baseline is None:
            baseline = self.integrate_baseline(interval,baseline_initial,rtol=rtol,
                atol=atol,method=method,max_step=max_step)
        if baseline.result.sol.t_min > interval[0] or baseline.result.sol.t_max < interval[1]:
            raise ValueError('baseline must cover entire finite-dose interval')
        initial = self.initial_delta() if initial is None else initial
        result = solve_ivp(lambda t,y:self.delta_rhs(t,y,baseline),interval,self.pack_delta(initial),
                           method=method,rtol=rtol,atol=atol,max_step=max_step,
                           t_eval=output_times,dense_output=dense_output)
        if not result.success:
            raise RuntimeError(f'finite-dose integration failed: {result.message}')
        return FiniteDoseSolution(self,baseline,result)


def compare_observables(coarse, fine, times, readout):
    """Consumer-selected physical observables, not an opaque state-norm gate."""
    times = np.asarray(times,dtype=float)
    left = np.asarray([readout(coarse,float(t)) for t in times])
    right = np.asarray([readout(fine,float(t)) for t in times])
    difference = np.abs(left-right)
    return dict(times=times,coarse=left,fine=right,absolute_error=difference,
                relative_error=difference/np.maximum(1.,np.abs(right)))
