"""Independent finite oracle; parent runs this file, children do not execute it.

Run from repository root:
  python proposals/P254-ni-collective-response/attempts/0029/oracle_reduced_electronic_correlations.py

The reference embeds each q edge in a global 12-mode one-body space and forms
its dense finite linearization by applying global commutators to basis vectors.
It independently evolves full finite transition matrices.  Initial electronic
noise is also checked by explicit sparse Fock operators in the eigenbasis of
the quasifree product state, NOT by forwarding the production Wick routine.
Only this oracle uses a full transition matrix or a finite Fock space.
"""
from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, expm
from scipy.sparse import coo_matrix, csr_matrix

from reduced_electronic_correlations import ActiveBandCorrelations


def symbolic_transpose():
    """Exact complex-polynomial check; no assumption of real/Hermitian entries."""
    def matrix(name):
        return sp.Matrix(2, 2, sp.symbols(f'{name}0:4'))
    ht, hs, rt, rs, x, v, z = [matrix(s) for s in ('ht', 'hs', 'rt', 'rs', 'x', 'v', 'z')]
    tangent = -sp.I * (ht*x-x*hs+v*rs-rt*v)
    zx = -sp.I * (ht.T*z-z*hs.T)
    zv = -sp.I * (z*rs.T-rt.T*z)
    residue = sum(z[a,b]*tangent[a,b]-zx[a,b]*x[a,b]-zv[a,b]*v[a,b]
                  for a in range(2) for b in range(2))
    if sp.expand(residue) != 0:
        raise AssertionError('symbolic bilinear transpose identity failed')


def canonical(n):
    return np.block([[np.zeros((n,n)), np.eye(n)], [-np.eye(n), np.zeros((n,n))]])


def fixture(n=3):
    rng = np.random.default_rng(290031)
    ns, m, d, f = 2*n, 2, 3, 1
    c = SimpleNamespace(mesh=(n,1,1), n=n, nk=n, nq=n, ns=ns, m=m, d=d, f=f, r0=2)
    c.q_indices = np.column_stack((np.arange(n), np.zeros((n,2), int)))
    c.k_indices = c.q_indices.copy()
    c.slot_spin = np.repeat(np.arange(2), n)
    c.slot_k = np.tile(np.arange(n), 2)
    c.target = np.array([[spin*n+(k+q)%n for spin,k in zip(c.slot_spin,c.slot_k)]
                         for q in range(n)])
    c.minus = (-np.arange(n)) % n
    c.rho0, c.h0 = [], []
    for s in range(ns):
        a = rng.normal(size=(m,m))+1j*rng.normal(size=(m,m))
        u, _ = np.linalg.qr(a)
        c.rho0.append(u@np.diag([.16+.025*s, .79-.021*s])@u.conj().T)
        a = rng.normal(size=(m,m))+1j*rng.normal(size=(m,m))
        c.h0.append(.18*(a+a.conj().T)+np.diag([-.2,.3]))
    c.rho0, c.h0 = np.array(c.rho0), np.array(c.h0)
    field = rng.normal(size=(ns,f,m,m))+1j*rng.normal(size=(ns,f,m,m))
    c.field_modes = .16*(field+field.conj().swapaxes(-1,-2))
    c.field_reverse_modes = c.field_modes.swapaxes(-1,-2)
    generators, js, w0, stiffness = [], [], [], []
    for q in range(n):
        nx = 2*d+(2*f if q == 0 else 0)
        j = block_diag(canonical(d), canonical(f)) if q == 0 else canonical(d)
        a = rng.normal(size=(nx,nx))
        h = .025*(a+a.T)+np.diag(np.linspace(.7,1.1,nx))
        if q > c.minus[q]:
            h = stiffness[c.minus[q]].T
        elif 0 < q < c.minus[q]:
            b = rng.normal(size=(nx,nx))
            h = h+.014j*(b-b.T)
        stiffness.append(h)
        w = np.diag(np.linspace(.8,1.2,nx))+.5j*j
        generators.append(j@h)
        js.append(j.astype(complex))
        w0.append(w)
    c.mechanical = SimpleNamespace(generators=tuple(generators), j=tuple(js), w0=tuple(w0))
    blocks = []
    for q in range(n):
        ionic = .2*(rng.normal(size=(ns,d,m,m))+1j*rng.normal(size=(ns,d,m,m)))
        potential = np.zeros((ns,m,m,2), complex)
        for s in range(ns):
            v = .2*(rng.normal(size=(m,m))+1j*rng.normal(size=(m,m)))
            if q == 0:
                v = v+v.conj().T
                ionic[s] = ionic[s]+ionic[s].conj().swapaxes(-1,-2)
            potential[s,:,:,c.slot_spin[s]] = v
        if q > c.minus[q]:
            for s in range(ns):
                reverse_source = c.target[q,s]
                ionic[s] = blocks[c.minus[q]].ionic[reverse_source].conj().swapaxes(-1,-2)
                potential[s] = blocks[c.minus[q]].potential[reverse_source].conj().swapaxes(0,1)
        elif q != 0 and q == c.minus[q]:
            ionic = ionic+ionic[c.target[q]].conj().swapaxes(-1,-2)
            potential = potential+potential[c.target[q]].conj().swapaxes(1,2)
        density = potential.conj().transpose(0,3,1,2)
        blocks.append(SimpleNamespace(ionic=ionic, ionic_reverse=ionic.conj(),
                                      potential=potential, density=density,
                                      kernel=np.array([[.24,.073],[.073,.19]], complex),
                                      feature_partner=np.arange(2)))
    c.q = lambda q: blocks[q]
    c.feature_count = lambda q: blocks[q].density.shape[1]

    def baseline(time):
        phase = time+.23*np.sin(time)
        rho = np.array([expm(-1j*h*phase)@r@expm(1j*h*phase)
                        for h,r in zip(c.h0,c.rho0)])
        return rho, (1+.23*np.cos(time))*c.h0
    return c, baseline, rng


class GlobalLinearization:
    """Independent dense embedding, never calling consumer tangent/transpose."""

    def __init__(self, c, baseline, q, *, omit_root_n=False, omit_field=False):
        self.c, self.baseline, self.q = c, baseline, q
        self.nx = len(c.mechanical.generators[q])
        self.size = self.nx+c.ns*c.m*c.m
        self.total = c.ns*c.m
        self.g = []
        self.gr = []
        self.p = []
        self.dr = []
        b = c.q(q)
        norm = 1.0 if omit_root_n else np.sqrt(c.n)
        for i in range(c.d):
            self.g.append(self.embed(b.ionic[:,i])/norm)
            self.gr.append(self.readout(b.ionic_reverse[:,i])/norm)
        for r in range(b.density.shape[1]):
            self.p.append(self.embed(b.potential[:,:,:,r])/norm)
            self.dr.append(self.readout(b.density[:,r])/norm)
        self.field = []
        self.field_reverse = []
        if q == 0:
            for i in range(c.f):
                factor = 0.0 if omit_field else 1.0
                self.field.append(self.embed(c.field_modes[:,i])*factor)
                self.field_reverse.append(self.readout(c.field_reverse_modes[:,i])*factor)

    def embed(self, sigma):
        c = self.c
        full = np.zeros((self.total,self.total), complex)
        for s in range(c.ns):
            t = c.target[self.q,s]
            full[t*c.m:(t+1)*c.m, s*c.m:(s+1)*c.m] = sigma[s]
        return full

    def readout(self, coefficients):
        """Operator M such that trace(M sigma_global) is elementwise readout."""
        return self.embed(coefficients).T

    def extract(self, matrix):
        c = self.c
        return np.array([matrix[t*c.m:(t+1)*c.m,s*c.m:(s+1)*c.m]
                         for s,t in enumerate(c.target[self.q])])

    def matrix(self, time):
        c = self.c
        rho, h = self.baseline(time)
        rho, h = block_diag(*rho), block_diag(*h)
        result = np.zeros((self.size,self.size), complex)
        for column in range(self.size):
            state = np.zeros(self.size, complex)
            state[column] = 1
            x = state[:self.nx]
            sigma = self.embed(state[self.nx:].reshape(c.ns,c.m,c.m))
            n = np.array([np.trace(d@sigma) for d in self.dr])
            potential_coefficients = c.q(self.q).kernel@n
            dh = np.zeros_like(h)
            for i,g in enumerate(self.g):
                dh += x[i]*g
            for r,p in enumerate(self.p):
                dh += potential_coefficients[r]*p
            gradient = np.zeros(self.nx, complex)
            for i,gr in enumerate(self.gr):
                gradient[i] = np.trace(gr@sigma)
            for i,field in enumerate(self.field):
                a = 2*c.d+i
                dh += x[a]*field
                gradient[a] = np.trace(self.field_reverse[i]@sigma)
            dx = c.mechanical.generators[self.q]@x+c.mechanical.j[self.q]@gradient
            ds = -1j*(h@sigma-sigma@h+dh@rho-rho@dh)
            result[:,column] = np.concatenate((dx,self.extract(ds).ravel()))
        return result

    def transition(self, time):
        identity = np.eye(self.size, dtype=complex)
        if time == 0:
            return identity
        solution = solve_ivp(lambda t,y: (self.matrix(t)@y.reshape(self.size,self.size)).ravel(),
                             (0.,time), identity.ravel(), method='DOP853',
                             rtol=5e-11, atol=5e-13, t_eval=[time])
        if not solution.success:
            raise AssertionError(solution.message)
        return solution.y[:,-1].reshape(self.size,self.size)

    def observable(self, dual):
        return self.readout(dual[self.nx:].reshape(self.c.ns,self.c.m,self.c.m))


def direct_quasifree(c, q, left, right):
    """Global one-body operator products; no directed-index Wick routine."""
    other = int(c.minus[q])
    # Baseline is irrelevant for this initial operator embedding.
    left_embed = GlobalLinearization(c, None, q)
    right_embed = GlobalLinearization(c, None, other)
    rho = block_diag(*c.rho0)
    hole = np.eye(len(rho))-rho
    lm = [left_embed.observable(left[:,p]) for p in range(left.shape[1])]
    rm = [right_embed.observable(right[:,r]) for r in range(right.shape[1])]
    electronic_ab = np.array([[np.trace(a@hole@b@rho) for b in rm] for a in lm])
    electronic_ba = np.array([[np.trace(b@hole@a@rho) for b in rm] for a in lm])
    lx, rx = left[:left_embed.nx], right[:right_embed.nx]
    ab = lx.T@c.mechanical.w0[q]@rx+electronic_ab
    ba = (rx.T@c.mechanical.w0[other]@lx).T+electronic_ba
    return ab, ba, electronic_ab, lm, rm


class DirectFockProduct:
    """Explicit CAR matrices and probability-weighted operator products."""

    def __init__(self, rho):
        occupations, self.u = np.linalg.eigh(rho)
        modes = len(occupations)
        dimension = 1 << modes
        states = np.arange(dimension, dtype=np.int64)
        probability = np.ones(dimension)
        self.annihilation = []
        for mode in range(modes):
            occupied = (states >> mode)&1
            probability *= np.where(occupied, occupations[mode], 1-occupations[mode])
            columns = states[occupied.astype(bool)]
            rows = columns ^ (1 << mode)
            signs = np.array([(-1)**int(int(s)&((1 << mode)-1)).bit_count() for s in columns])
            self.annihilation.append(coo_matrix((signs,(rows,columns)),
                                               shape=(dimension,dimension)).tocsr())
        self.probability = probability
        self.dimension = dimension

    def operator(self, one_body):
        matrix = self.u.conj().T@one_body@self.u
        out = csr_matrix((self.dimension,self.dimension), dtype=complex)
        for i,ci in enumerate(self.annihilation):
            for j,cj in enumerate(self.annihilation):
                if matrix[i,j] != 0:
                    out += matrix[i,j]*(ci.conj().T@cj)
        return out

    def connected(self, a, b):
        oa, ob = self.operator(a), self.operator(b)
        mean_a = self.probability@oa.diagonal()
        mean_b = self.probability@ob.diagonal()
        ordered = np.asarray(oa.multiply(ob.T).sum(axis=1)).ravel()
        return self.probability@ordered-mean_a*mean_b


def relative(value, reference):
    return float(np.linalg.norm(value-reference)/max(1.,np.linalg.norm(reference)))


def short_interval_baseline_domain(consumer_type=ActiveBandCorrelations):
    """A valid short pullback must not query unavailable baseline history."""
    c, baseline, _ = fixture()
    time = 1e-6
    rho, h = baseline(0.)

    def bounded(time_at):
        if not 0. <= time_at <= time:
            raise ValueError('baseline outside its actual integrated interval')
        return rho, h

    consumer = consumer_type(c, bounded, rtol=1e-10, atol=1e-12)
    probes = consumer.ionic_probes(0, (0,))
    expected = expm(time*GlobalLinearization(c, bounded, 0).matrix(0.).T) @ probes
    actual = consumer.backward_probes(0, time, (0,))
    error = relative(actual, expected)
    if error > 2e-12:
        raise AssertionError('short-domain pullback differs from independent transition')
    return error, relative(probes, expected)


def spatial_batch_regressions(c, consumer, duals_left, duals_right, t, s, cells):
    """Independent dense/global products expose time, coordinate and q collisions."""
    errors = {}
    phases = np.exp(2j*np.pi*(cells/np.asarray(c.mesh))@c.q_indices.T)
    # Reordering AND repeating coordinates must survive the cache and both axes.
    coordinates, columns = (2,0,2), (1,0,1)
    left = [dual[:,columns] for dual in duals_left]
    right = [dual[:,columns] for dual in duals_right]

    def expected(left_duals, right_duals):
        ab, ba = [], []
        for q in range(c.nq):
            a, b, _, _, _ = direct_quasifree(c, q, left_duals[q],
                                            right_duals[int(c.minus[q])])
            ab.append(a)
            ba.append(b)
        return (np.tensordot(phases, np.array(ab), axes=(1,0))/c.n,
                np.tensordot(phases, np.array(ba), axes=(1,0))/c.n)

    def compare(name, actual, reference):
        ab, ba = reference
        for component, value in (
                ('ordered_ab', ab), ('ordered_ba', ba),
                ('symmetrized', (ab+ba)/2), ('commutator', ab-ba)):
            errors[f'{name}_{component}'] = relative(getattr(actual, component), value)

    equal_reference, two_reference = expected(left, left), expected(left, right)
    equal, two = consumer.spatial_pair(t, s, cells, coordinates)
    compare('batch_equal', equal, equal_reference)
    compare('batch_two', two, two_reference)
    same_equal, same_two = consumer.spatial_pair(t, t, cells, coordinates)
    compare('batch_same_time_equal', same_equal, equal_reference)
    compare('batch_same_time_two', same_two, equal_reference)
    # Mutating any returned component must not alter either the other component
    # or a future call, even when all pullbacks were identical.
    for name in ('ordered_ab', 'ordered_ba', 'symmetrized', 'commutator'):
        getattr(same_equal, name)[0,0,0] += 17+3j
    compare('batch_return_independence', same_two, equal_reference)
    fresh_equal, fresh_two = consumer.spatial_pair(t, s, cells, coordinates)
    compare('batch_fresh_equal', fresh_equal, equal_reference)
    compare('batch_fresh_two', fresh_two, two_reference)
    # Distinct left/right lengths and orders collide at self-minus q if the
    # ordered coordinate tuple is absent from the pullback key.
    short_right = [dual[:,(0,)] for dual in duals_left]
    distinct_reference = expected(left, short_right)
    distinct = consumer.spatial(t, t, cells, coordinates, (0,))
    compare('spatial_distinct_axes', distinct, distinct_reference)
    for q in np.flatnonzero(c.minus == np.arange(c.nq)):
        ab, ba, _, _, _ = direct_quasifree(c, q, left[q], short_right[q])
        actual = consumer.two_time_q(int(q), t, t, coordinates, (0,))
        compare(f'direct_self_minus_q{q}', actual, (ab, ba))
    zero = []
    for q in range(c.nq):
        probes = np.zeros_like(left[q])
        probes[np.asarray(coordinates), np.arange(len(coordinates))] = 1
        zero.append(probes)
    zero_equal, zero_two = consumer.spatial_pair(t, 0., cells, coordinates)
    compare('batch_zero_reference_equal', zero_equal, equal_reference)
    compare('batch_zero_reference_two', zero_two, expected(left, zero))
    # A separate reversed-coordinate request cannot reuse stale earlier columns.
    reversed_equal, reversed_two = consumer.spatial_pair(t, s, cells, (0,2))
    compare('batch_changed_coordinates_equal', reversed_equal, expected(duals_left, duals_left))
    compare('batch_changed_coordinates_two', reversed_two, expected(duals_left, duals_right))
    # The same q/time/coordinates on the same consumer must observe a changed
    # baseline on the next call, rather than retrieving a persistent dual.
    original_baseline = consumer.baseline
    fixed_rho, fixed_h = original_baseline(0.)
    frozen = lambda _: (fixed_rho, fixed_h)
    frozen_left, frozen_right = [], []
    for q in range(c.nq):
        dense = GlobalLinearization(c, frozen, q).matrix(0.)
        frozen_left.append(expm(t*dense).T@zero[q])
        frozen_right.append(expm(s*dense).T@zero[q])
    consumer.baseline = frozen
    try:
        changed_equal, changed_two = consumer.spatial_pair(t, s, cells, coordinates)
    finally:
        consumer.baseline = original_baseline
    compare('batch_changed_baseline_equal', changed_equal, expected(frozen_left, frozen_left))
    compare('batch_changed_baseline_two', changed_two, expected(frozen_left, frozen_right))
    return errors


def run_oracle():
    symbolic_transpose()
    c, baseline, rng = fixture()
    consumer = ActiveBandCorrelations(c, baseline, rtol=5e-11, atol=5e-13)
    errors, mutations = {}, {}
    errors['short_interval_pullback'], mutations['unpropagated_short_probe'] = short_interval_baseline_domain()
    t, s = .43, .29
    coordinates = (0,2)
    refs = [GlobalLinearization(c,baseline,q) for q in range(c.nq)]
    duals_left, duals_right = [], []
    for q,reference in enumerate(refs):
        dense = reference.matrix(.31)
        state = rng.normal(size=reference.size)+1j*rng.normal(size=reference.size)
        z = rng.normal(size=reference.size)+1j*rng.normal(size=reference.size)
        actual = consumer.tangent(q,.31,state)
        transpose = consumer.transpose(q,.31,z)
        errors[f'global_tangent_q{q}'] = relative(actual,dense@state)
        errors[f'global_transpose_q{q}'] = relative(transpose,dense.T@z)
        errors[f'bilinear_identity_q{q}'] = float(abs(z@actual-transpose@state)/max(1.,abs(z@actual)))
        mutations[f'dagger_not_transpose_q{q}'] = relative(transpose,dense.conj().T@z)
        mutations[f'missing_physical_root_n_q{q}'] = relative(
            actual,GlobalLinearization(c,baseline,q,omit_root_n=True).matrix(.31)@state)
        probes = consumer.ionic_probes(q,coordinates)
        dl, dr = reference.transition(t).T@probes, reference.transition(s).T@probes
        duals_left.append(dl)
        duals_right.append(dr)
        errors[f'backward_nonstationary_q{q}'] = relative(consumer.backward_probes(q,t,coordinates),dl)
        mutations[f'frozen_baseline_q{q}'] = relative(dl,expm(t*reference.matrix(0)).T@probes)
    fock = DirectFockProduct(block_diag(*c.rho0))
    reference_ab, reference_ba = [], []
    for q in range(c.nq):
        minus = int(c.minus[q])
        left, right = duals_left[q], duals_right[minus]
        ab, ba, electronic, lm, rm = direct_quasifree(c,q,left,right)
        reference_ab.append(ab)
        reference_ba.append(ba)
        measured = consumer.two_time_q(q,t,s,coordinates)
        errors[f'ordered_ab_q{q}'] = relative(measured.ordered_ab,ab)
        errors[f'ordered_ba_q{q}'] = relative(measured.ordered_ba,ba)
        errors[f'symmetrized_q{q}'] = relative(measured.symmetrized,(ab+ba)/2)
        errors[f'commutator_q{q}'] = relative(measured.commutator,ab-ba)
        # Use generic complex electronic probes so Wick mutations cannot hide
        # behind a displacement probe's short-time electronic smallness.
        electronic_left = np.zeros_like(left)
        electronic_right = np.zeros_like(right)
        nx = refs[q].nx
        nrm = refs[minus].nx
        electronic_left[nx:] = rng.normal(size=left[nx:].shape)+1j*rng.normal(size=left[nx:].shape)
        electronic_right[nrm:] = rng.normal(size=right[nrm:].shape)+1j*rng.normal(size=right[nrm:].shape)
        wab, wba, wel, olm, orm = direct_quasifree(c,q,electronic_left,electronic_right)
        errors[f'global_wick_q{q}'] = relative(consumer.initial_ordered(q,electronic_left,electronic_right),wab)
        if q in (0,1):
            explicit = np.array([[fock.connected(a,b) for b in orm] for a in olm])
            reverse = np.array([[fock.connected(b,a) for b in orm] for a in olm])
            errors[f'explicit_fock_ab_q{q}'] = relative(wab,explicit)
            errors[f'explicit_fock_ba_q{q}'] = relative(wba,reverse)
        rho = block_diag(*c.rho0)
        no_hole = np.array([[np.trace(a@b@rho) for b in orm] for a in olm])
        mutations[f'missing_pauli_hole_q{q}'] = relative(wel,no_hole)
        mutations[f'extra_sector_weight_q{q}'] = relative(wel,wel/c.n)
        mutations[f'extra_spin_weight_q{q}'] = relative(wel,wel/2)
        if q == 1:
            ls = electronic_left[nx:].reshape(c.ns,c.m,c.m,-1)
            rs = electronic_right[nrm:].reshape(c.ns,c.m,c.m,-1)
            wrong = np.zeros_like(wel)
            for slot in range(c.ns):
                target = c.target[q,slot]
                for a in range(c.m):
                    for b in range(c.m):
                        for cc in range(c.m):
                            for dd in range(c.m):
                                wrong += np.outer(ls[slot,a,b],rs[slot,cc,dd])*c.rho0[slot,cc,b]*(
                                    (1 if a == dd else 0)-c.rho0[target,a,dd])
            mutations['unreversed_opposite_source'] = relative(wel,wrong)
    cells = np.array([[0,0,0],[1,0,0],[2,0,0]])
    phases = np.exp(2j*np.pi*(cells/np.array(c.mesh))@c.q_indices.T)
    expected_ab = np.tensordot(phases,np.array(reference_ab),axes=(1,0))/c.n
    expected_ba = np.tensordot(phases,np.array(reference_ba),axes=(1,0))/c.n
    measured = consumer.spatial(t,s,cells,coordinates)
    errors['spatial_ab'] = relative(measured.ordered_ab,expected_ab)
    errors['spatial_ba'] = relative(measured.ordered_ba,expected_ba)
    errors['spatial_commutator'] = relative(measured.commutator,expected_ab-expected_ba)
    mutations['missing_fourier_inverse_n'] = relative(expected_ab,c.n*expected_ab)
    no_field = GlobalLinearization(c,baseline,0,omit_field=True)
    probes = consumer.ionic_probes(0,coordinates)
    mutations['missing_gamma_field'] = relative(duals_left[0],no_field.transition(t).T@probes)
    # Equal-time entrypoint against independent transitions, with nonzero field
    # preparation and both actual spin populations present.
    eq_ab, eq_ba, _, _, _ = direct_quasifree(c,1,duals_left[1],duals_left[2])
    equal = consumer.equal_time_q(1,t,coordinates)
    errors['equal_time_ab'] = relative(equal.ordered_ab,eq_ab)
    errors['equal_time_ba'] = relative(equal.ordered_ba,eq_ba)
    errors.update(spatial_batch_regressions(
        c, consumer, duals_left, duals_right, t, s, cells))
    # An even mesh adds a non-Gamma self-minus (Nyquist) sector; do not
    # mistake its two ordered coordinate requests for the same probe bank.
    even, even_baseline, _ = fixture(n=4)
    even_consumer = ActiveBandCorrelations(even, even_baseline,
                                          rtol=5e-11, atol=5e-13)
    even_left, even_right = [], []
    for q in range(even.nq):
        reference = GlobalLinearization(even, even_baseline, q)
        probes = np.zeros((reference.size, len(coordinates)), complex)
        probes[np.asarray(coordinates), np.arange(len(coordinates))] = 1
        even_left.append(reference.transition(t).T@probes)
        even_right.append(reference.transition(s).T@probes)
    even_errors = spatial_batch_regressions(
        even, even_consumer, even_left, even_right, t, s, cells)
    errors.update({f'even_mesh_{name}': value for name, value in even_errors.items()})
    failures = {name:error for name,error in errors.items() if not np.isfinite(error) or error > 2e-8}
    invisible = {name:gap for name,gap in mutations.items() if not np.isfinite(gap) or gap < 1e-7}
    report = dict(schema='P254-active-correlations-independent-finite-oracle-v1',
                  errors=errors, mutation_gaps=mutations, symbolic_transpose='exact',
                  fixture=dict(mesh=list(c.mesh), spins=2, bands=c.m,
                               global_fermion_modes=c.ns*c.m, gamma_field_modes=c.f,
                               nonstationary_baseline=True),
                  noise_scope=measured.noise_scope)
    if failures or invisible:
        raise AssertionError(json.dumps(dict(errors=failures,invisible_mutations=invisible),sort_keys=True))
    return report


if __name__ == '__main__':
    print(json.dumps(run_oracle(),indent=2,sort_keys=True))
