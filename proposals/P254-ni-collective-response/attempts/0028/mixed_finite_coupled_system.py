"""Finite Gamma coupled dynamics and independent nonlinear density comparator.

Only the tiny independent oracle may form PW-square matrices. The dynamics
provider uses the admitted borrowed density descriptors and endpoint frame RHS.
The comparator never calls either: it evolves full density commutators and
canonical derivatives of a SymPy scalar, with physical parameter continuations.
The fixture retains its original single canonical photon mode, including actual
one-mode checkpoints; no padding or artificial polarizations are introduced.
No Native/material/nonGamma result is implied by this fixture.
"""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import sympy as sp

from run_joint_response_oracle import FiniteSystem
from mixed_frame_response import density_dyads_item, frame_rhs_item
from mixed_coordinate_response import coordinate_rhs
from mixed_excitation_energy import mixed_components
from mixed_response_storage import COMPONENTS, SMALL_NAMES, MappedMixedState, publish_checkpoint
from substrate_framework.partitioned_q_storage import durable_publish, sync_file
from substrate_framework.partitioned_q_work import arrays
from substrate_framework.streamed_material_state import authenticate, file_binding

ORDERS = ('mean', 'plus', 'minus', 'mixed')
ENERGIES = ('bare_electronic', 'hartree', 'electron_ionic',
            'ionic_reference_linear', 'ionic_C', 'ionic_kinetic_minimal',
            'delta_C', 'electron_field', 'field_contact', 'ion_field_mixed', 'maxwell')


def small_vector(small):
    return np.concatenate([small.u, small.p, small.coordinate[:1], small.momentum[:1]])


def small_dict(vector):
    return dict(u=np.asarray(vector[:3]), p=np.asarray(vector[3:6]),
                coordinate=np.asarray(vector[6:7]), momentum=np.asarray(vector[7:8]))


def model_metadata(s):
    reference = NS(ncells=s.n, masses=s.mass, curvature=s.c, nbar=.04,
                   ionic_force=np.zeros(3), support_force=np.zeros(3))
    coefficients = NS(delta_curvature=s.dc, contact=s.contact, mixed=s.mixed)
    return reference, coefficients, NS(mass=s.mf)


def contractions(block, operator):
    # Weight and original D occur exactly once; no density or Gram inverse.
    return sum(d.weight*np.trace(d.right.conj().T @ operator @ d.left @ d.D)
               for d in block.dyads)


def local_moments(s, index, descriptors):
    return np.asarray([[contractions(descriptors[c], s.v[index]),
                        *(contractions(descriptors[c], g) for g in s.g[index]),
                        *(contractions(descriptors[c], l) for l in s.ell[index])]
                       for c in ORDERS], complex)


class PreparedFinite:
    def __init__(self, binding, state, budget):
        self.binding = dict(binding)
        record = json.loads(authenticate(binding).read_text())
        if (record['schema'] != 'P254-finite-coupled-context-v1'
                or record['state'] != state.binding):
            raise ValueError('finite context is not bound to the exact stage')
        self._record, self._budget = record, budget
        self._arrays = {}
        for name, entry in record['arrays'].items():
            self._arrays[name] = np.load(authenticate(entry), allow_pickle=False)
        self.small_rhs = {c: {name: self._arrays[c+'-'+name]
                              for name in SMALL_NAMES} for c in COMPONENTS}

    def edge_rhs(self, index, state):
        if state.binding != self._record['state']:
            raise ValueError('finite prepared context used with another stage')
        values = self._arrays
        psi = values['psi-'+str(index)]
        h = {c: values[c+'-h-'+str(index)] for c in ORDERS}
        if self._record['mutation'] == 'missing_mixed_source':
            h['mixed'] = np.zeros_like(h['mixed'])
        def apply(component, x):
            self._budget.action(x)
            return h[component] @ x
        rhs = frame_rhs_item(index, psi, *(state.leg(c, index) for c in COMPONENTS),
            plus_endpoint=index, minus_endpoint=index,
            h0=lambda i, x: apply('mean', x),
            hplus=lambda i, x: apply('plus', x),
            hminus=lambda i, x: apply('minus', x),
            hmixed=lambda i, x: apply('mixed', x))
        if self._record['mutation'] == 'half_mixed_cross_source':
            cross = -1j*(apply('plus', state.leg('minus', index))
                          +apply('minus', state.leg('plus', index)))
            rhs = (rhs[0], rhs[1], rhs[2]-.5*cross)
        return dict(zip(COMPONENTS, rhs))


class FiniteDynamics:
    """Durable prepared interface; global reduction completes before retirement."""
    def __init__(self, system, directory, budget, sources, mutation='none'):
        self.s, self.directory, self.budget = system, Path(directory), budget
        self.sources, self.mutation = sources, mutation
        self.directory.mkdir(parents=True, exist_ok=True)

    def prepare(self, state, work):
        s = self.s
        psi, u, p, a, pi, rho, hs = s.mean_values(s.mean.sol(state.time))
        # This finite baseline is the exact moving fixture, not a frozen mean.
        readers = {c: _LegReader(state, c, len(psi)) for c in COMPONENTS}
        weights = [s.w]*len(psi)
        if self.mutation == 'wrong_source_weight':
            weights[1] *= 1.7
        identity = dict(stage=state.binding, sources=self.sources, mutation=self.mutation)
        def reduce_edge(index):
            desc = density_dyads_item(index, psi, readers['plus'], readers['minus'],
                readers['mixed'], s.ds, target_plus=state.target_plus,
                target_minus=state.target_minus, source_weights=weights)
            result = local_moments(s, index, desc)
            state.release_edge(index)
            return dict(moments=result)
        records = work.execute('finite-coupled-global-moments', identity, len(psi), reduce_edge)
        total = np.zeros((4, 7), complex)
        for record in records:
            mapped = arrays(record)
            total += mapped['moments']
            for value in mapped.values():
                value._mmap.close()
        reference, coefficients, field = model_metadata(s)
        payloads = {}
        for j, c in enumerate(COMPONENTS, 1):
            raw = state.small[c]
            active = NS(u=raw.u, p=raw.p, coordinate=raw.coordinate[:1],
                        momentum=raw.momentum[:1])
            rhs = coordinate_rhs(reference=reference, coefficients=coefficients,
                field=field, charges=[s.z], polarizations=s.e, curvature=s.c,
                gamma=True, small=active, gradient=total[j, 1:4],
                field_gradient=total[j, 4:7])
            if self.mutation == 'wrong_canonical_N':
                rhs['momentum'] = rhs['momentum']/s.n
            if self.mutation == 'wrong_canonical_charge':
                rhs['u'] = raw.p/s.mass
                rhs['momentum'] = -s.n*(s.e @ total[j, 4:7]
                    +s.e @ s.mixed.T @ raw.u+s.e @ s.contact @ s.e.T @ active.coordinate)
            for name in SMALL_NAMES:
                payloads[c+'-'+name] = rhs[name]
        for i, frame in enumerate(psi):
            payloads['psi-'+str(i)] = frame
            payloads['mean-h-'+str(i)] = hs[i]
            for j, c in enumerate(COMPONENTS, 1):
                raw = state.small[c]
                payloads[c+'-h-'+str(i)] = (s.k*total[j, 0]*s.v[i]
                    +np.einsum('a,aij->ij', raw.u, s.g[i])
                    +np.einsum('a,aij->ij', s.e.T @ raw.coordinate[:1], s.ell[i]))
        key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        root = self.directory/key
        root.mkdir(exist_ok=True)
        target = root/'context.json'
        if target.exists():
            return self.restore(file_binding(target), state, work)
        bindings = {}
        for name, value in payloads.items():
            bindings[name] = self.budget.array(root/(name+'.npy'), np.asarray(value, complex))
            sync_file(bindings[name]['path'])
        binding = durable_publish(target, dict(schema='P254-finite-coupled-context-v1',
            state=state.binding, mutation=self.mutation, sources=self.sources,
            arrays=bindings, global_moment_receipts=records), self.budget)
        return self.restore(binding, state, work)

    def restore(self, binding, state, work):
        prepared = PreparedFinite(binding, state, self.budget)
        if (prepared._record['sources'] != self.sources
                or prepared._record['mutation'] != self.mutation):
            raise ValueError('finite context scientific source identity changed')
        return prepared


class _LegReader:
    def __init__(self, state, component, count):
        self.state, self.component, self.count = state, component, count
    def __len__(self):
        return self.count
    def __getitem__(self, index):
        return self.state.leg(self.component, index)


class CoupledFixture:
    def __init__(self):
        self.system = FiniteSystem()
        self.start, self.stop = .073, .193
        s = self.system
        self.frames, self.u, self.p, self.a, self.pi, self.rho, _ = s.mean_values(s.mean.sol(self.start))
        rng = np.random.default_rng(2540028)
        self.A, self.B, self.C = [], [], []
        for n in s.sizes:
            a = .18*(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
            c = .11*(rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)))
            self.A.append(a)
            self.B.append(-a.conj().T)
            self.C.append(c-c.conj().T)
        self.first = .13*(rng.normal(size=8)+1j*rng.normal(size=8))
        self.cross = .17*rng.normal(size=8)
        self.mean_small = np.concatenate([self.u, self.p, self.a, self.pi])
        self.symbolic = SymbolicScalar(s)

    def initial(self, directory, budget, sources):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=False)
        vals = {c: [] for c in COMPONENTS}
        for psi, a, b, c in zip(self.frames, self.A, self.B, self.C):
            vals['plus'].append(a @ psi)
            vals['minus'].append(b @ psi)
            vals['mixed'].append((c+(a @ b+b @ a)/2) @ psi)
        legs = {c: [budget.array(directory/f'{c}-{i}.npy', value)
                    for i, value in enumerate(vals[c])] for c in COMPONENTS}
        small = {c: small_dict(value) for c, value in zip(COMPONENTS,
                 (self.first, self.first.conj(), self.cross))}
        return publish_checkpoint(directory/'checkpoint.json', legs=legs, small=small,
            q=[0., 0., 0.], target_plus=[0, 1, 2], target_minus=[0, 1, 2],
            time=self.start, identity='finite-self-consistent-coupled-Gamma',
            kind='physical-exponential-initial-parameter-jet',
            provenance=dict(scope='finite only; not Native/material', sources=sources,
                            physical_curve='exp(lambda A+eta B+lambda eta C); B=-A^dagger; C anti-Hermitian'),
            budget=budget)

    def pack_density(self, rho, small):
        return np.concatenate([*(r.ravel() for r in rho), small])

    def unpack_density(self, vector):
        rho, start = [], 0
        for n in self.system.sizes:
            rho.append(vector[start:start+n*n].reshape(n, n))
            start += n*n
        return rho, vector[start:]

    def continuation_initial(self, lam):
        eta = lam.conjugate()
        rho = []
        for r, a, b, c in zip(self.rho, self.A, self.B, self.C):
            g = lam*a+eta*b+lam*eta*c
            # Independent analytic density continuation; no frame jets or D inverse.
            rho.append(expm(g) @ r @ expm(-g))
        small = self.mean_small+lam*self.first+eta*self.first.conj()+lam*eta*self.cross
        return self.pack_density(rho, small)

    def density_rhs(self, time, vector):
        s = self.system
        rho, small = self.unpack_density(vector)
        moments = self.symbolic.moments(rho)
        u, p, a, pi = small[:3], small[3:6], small[6], small[7]
        h = [h0+s.k*(moments[1]-.04)*v+np.einsum('a,aij->ij', u, g)
             +np.einsum('a,aij->ij', s.e.ravel()*a, ell)
             for h0, v, g, ell in zip(s.h0, s.v, s.g, s.ell)]
        drho = [-1j*(hi @ r-r @ hi) for hi, r in zip(h, rho)]
        dsmall = self.symbolic.canonical(small, moments)
        return self.pack_density(drho, dsmall)

    def integrate(self, lam):
        sol = solve_ivp(self.density_rhs, (self.start, self.stop),
            self.continuation_initial(lam), method='DOP853', rtol=3e-13,
            atol=3e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1]

    def readout(self, binding, budget, *, omit_mean_mixed=False):
        s = self.system
        with MappedMixedState(binding, budget=budget) as state:
            psi, u, p, a, pi, _, _ = s.mean_values(s.mean.sol(state.time))
            readers = {c: _LegReader(state, c, 3) for c in COMPONENTS}
            moments = np.zeros((4, 7), complex)
            one = np.zeros(4, complex)
            density_mixed = []
            for i in range(3):
                desc = density_dyads_item(i, psi, readers['plus'], readers['minus'],
                    readers['mixed'], s.ds, target_plus=state.target_plus,
                    target_minus=state.target_minus, source_weights=[s.w]*3)
                moments += local_moments(s, i, desc)
                one += [contractions(desc[c], s.h0[i]) for c in ORDERS]
                # Comparator-visible complete density, not merely n/g/ell moments.
                density_mixed.append(sum(d.left @ d.D @ d.right.conj().T
                                         for d in desc['mixed'].dyads))
            mean_small = NS(**small_dict(np.concatenate([u, p, a, pi])))
            stages = {}
            for j, c in enumerate(ORDERS):
                small = mean_small if c == 'mean' else state.small[c]
                stages[c] = NS(density=np.array([moments[j, 0]]),
                    potential=np.array([s.k*(moments[j, 0]-(.04 if c == 'mean' else 0))]),
                    u=small.u, p=small.p, shift=(s.e.T @ small.coordinate[:1]).ravel(),
                    momentum=small.momentum[:1], gradient=moments[j, 1:4],
                    field_gradient_electronic=moments[j, 4:7])
            if omit_mean_mixed:
                # Consumer-visible stationary shortcut mutation, not a source-text test.
                zero = stages['mean']
                for name in ('u', 'p', 'shift', 'momentum', 'gradient', 'field_gradient_electronic'):
                    setattr(zero, name, np.zeros_like(getattr(zero, name)))
                zero.density = np.array([.04+0j])
                zero.potential = np.zeros(1, complex)
            reference, coefficients, field = model_metadata(s)
            energy = mixed_components(reference=reference, coefficients=coefficients,
                field=field, charges=[s.z], stages=stages,
                one_body={'mixed': one[-1]}, grid_measure=1.,
                curvature_plus=s.c, curvature_minus=s.c, gamma=True)
            density = self.pack_density(density_mixed, small_vector(state.small['mixed']))
            return density, np.asarray([energy[k] for k in ENERGIES])


class SymbolicScalar:
    """One independent eleven-component scalar; SymPy derives ALL forces."""
    def __init__(self, s):
        self.s = s
        u = sp.symbols('u0:3'); p = sp.symbols('p0:3')
        a, pi = sp.symbols('a pi')
        b, density = sp.symbols('bare density')
        g = sp.symbols('g0:3'); ell = sp.symbols('ell0:3')
        self.variables = (*u, *p, a, pi, b, density, *g, *ell)
        U, P, G, L = map(sp.Matrix, (u, p, g, ell))
        shift = sp.Matrix(s.e.ravel())*a
        velocity = P-sp.Matrix(s.e.ravel())*s.z*a
        mat = lambda x: sp.Matrix(np.asarray(x).real)
        n = sp.Integer(s.n)
        self.expressions = (
            n*b, n*sp.Float(s.k)*(density-sp.Rational(1, 25))**2/2,
            n*(U.dot(G)), sp.Integer(0), n*(U.dot(mat(s.c)*U))/2,
            n*sum(velocity[i]**2/sp.Float(s.mass[i]) for i in range(3))/2,
            n*(U.dot(mat(s.dc)*U))/2, n*shift.dot(L),
            n*shift.dot(mat(s.contact)*shift)/2,
            n*U.dot(mat(s.mixed)*shift), pi**2/(2*sp.Float(s.mf)))
        energy = sum(self.expressions)
        # Per-cell ionic p,u have symplectic factor 1/N; global a,Pi do not.
        canonical = [*(sp.diff(energy, x)/n for x in p),
                     *(-sp.diff(energy, x)/n for x in u),
                     sp.diff(energy, pi), -sp.diff(energy, a)]
        self._energy = sp.lambdify(self.variables, self.expressions, 'numpy', cse=True)
        self._canonical = sp.lambdify(self.variables, canonical, 'numpy', cse=True)
        self.force_derivatives = [str(x) for x in canonical]

    def moments(self, rho):
        s = self.s
        return np.asarray([s.w*sum(np.trace(h @ r) for h, r in zip(s.h0, rho)),
            s.w*sum(np.trace(v @ r) for v, r in zip(s.v, rho)),
            *(s.w*sum(np.trace(g[i] @ r) for g, r in zip(s.g, rho)) for i in range(3)),
            *(s.w*sum(np.trace(l[i] @ r) for l, r in zip(s.ell, rho)) for i in range(3))], complex)

    def canonical(self, small, moments):
        return np.asarray(self._canonical(*small, *moments), complex)

    def energy(self, fixture, vector):
        rho, small = fixture.unpack_density(vector)
        return np.asarray(self._energy(*small, *self.moments(rho)), complex)
