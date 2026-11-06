"""Non-Abelian local composite candidate. Parent law: SU3 Wilson links + CAR
u,d,anti-u + local quantum trace-matrix deformation. Approximation: covariant
1+3+bar3 Peter-Weyl link space and conditional local vibronic states. No scalar
clock substitution, fixed epsilon-only colour sector or whole-baryon hopping.
"""
from __future__ import annotations
import itertools
import math
from collections import defaultdict
from typing import Any
import numpy as np
from scipy.linalg import eigh
from scipy import sparse

Wave = dict[int, complex]
FERMIONS = 72
FERMION_MASK = (1 << FERMIONS)-1
CLOUD_SHIFT = FERMIONS+30


def epsilon(a: int, b: int, c: int) -> int:
    if len({a, b, c}) < 3:
        return 0
    return 1 if (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] else -1


def generators() -> list[Any]:
    result = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        x = np.zeros((3, 3), dtype=complex); x[i, j] = x[j, i] = .5; result.append(x)
        y = np.zeros((3, 3), dtype=complex); y[i, j], y[j, i] = -.5j, .5j; result.append(y)
    result += [np.diag([.5, -.5, 0]), np.diag([1, 1, -2])/np.sqrt(12)]
    return result


def link_operators() -> tuple[dict[tuple[int, int, bool], Any], list[Any], list[Any]]:
    transport: dict[tuple[int, int, bool], Any] = {}
    for a, b in itertools.product(range(3), repeat=2):
        u = np.zeros((19, 19))
        u[1+3*a+b, 0] = 1/np.sqrt(3)
        u[0, 10+3*a+b] = 1/np.sqrt(3)
        for c, d in itertools.product(range(3), repeat=2):
            for e, f in itertools.product(range(3), repeat=2):
                u[10+3*e+f, 1+3*c+d] = epsilon(a, c, e)*epsilon(b, d, f)/2
        transport[a, b, False] = u
        transport[a, b, True] = u.T
    left, right = [], []
    for t in generators():
        l, r = np.zeros((19, 19), dtype=complex), np.zeros((19, 19), dtype=complex)
        l[1:10, 1:10] = -np.kron(t.T, np.eye(3)); l[10:, 10:] = np.kron(t, np.eye(3))
        r[1:10, 1:10] = np.kron(np.eye(3), t); r[10:, 10:] = -np.kron(np.eye(3), t.T)
        left.append(l); right.append(r)
    return transport, left, right


def fermion(bits: int, operators: list[tuple[bool, int]]) -> tuple[int, int] | None:
    sign = 1
    for create, mode in reversed(operators):
        occupied = bool((bits >> mode) & 1)
        if occupied == create:
            return None
        if (bits & ((1 << mode)-1)).bit_count() % 2:
            sign = -sign
        bits ^= 1 << mode
    return bits, sign


def inner(a: Wave, b: Wave) -> complex:
    if len(a) > len(b):
        return inner(b, a).conjugate()
    return sum((value.conjugate()*b.get(key, 0j) for key, value in a.items()), 0j)


def normalized(a: Wave) -> tuple[Wave, float]:
    norm = float(inner(a, a).real)
    if norm <= 0:
        raise ValueError('Actual null source, no dummy normalization')
    scale = 1/np.sqrt(norm)
    return {key: value*scale for key, value in a.items() if value != 0}, norm


class Law:
    def __init__(self, p: dict[str, Any]):
        self.p = p
        self.positions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)/(2*np.sqrt(2))
        self.edges = list(itertools.combinations(range(4), 2))
        self.faces = list(itertools.combinations(range(4), 3))
        self.U, self.left, self.right = link_operators()
        self.transitions = {(a, b, adjoint, state): [(int(out), float(matrix[out, state])) for out in np.flatnonzero(matrix[:, state])]
            for (a, b, adjoint), matrix in self.U.items() for state in range(19)}
        self.local_energies, self.local_vectors = [], []
        n = p['bare_boson_levels']
        b = np.diag(np.sqrt(np.arange(1, n+4)), 1)
        extended_x = (b+b.T)/np.sqrt(2)
        x = extended_x[:n, :n]
        x2 = (extended_x@extended_x)[:n, :n]
        x4 = np.linalg.matrix_power(extended_x, 4)[:n, :n]
        self.bare_X, self.bare_X2 = x, x2
        for load in range(19):
            h = p['omega']*np.diag(np.arange(n))+p['quartic']*x4-p['deformation_g']*load*x
            e, v = eigh(h, subset_by_index=[0, p['local_cloud_states']-1])
            for col in range(v.shape[1]):
                if v[np.argmax(np.abs(v[:, col])), col] < 0:
                    v[:, col] *= -1
            self.local_energies.append(e); self.local_vectors.append(v)
        if p['local_cloud_states'] != 2:
            raise ValueError('This explicit binary vibronic encoding owns exactly2 local states')
        self.overlap = {(a, b): self.local_vectors[a].T@self.local_vectors[b] for a in range(19) for b in range(19)}
        self.X = {(a, b): self.local_vectors[a].T@x@self.local_vectors[b] for a in range(19) for b in range(19)}
        self.X2 = [v.T@x2@v for v in self.local_vectors]
        self.top = [v[-1, :]**2 for v in self.local_vectors]
        self.sigma = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
        self.hops = []
        for edge, (i, j) in enumerate(self.edges):
            direction = self.positions[i]-self.positions[j]
            spin = -p['hopping']*np.eye(2)-1j*p['spin_orbit']*sum((direction[d]*self.sigma[d] for d in range(3)), np.zeros((2, 2), dtype=complex))
            for species in range(3):
                charge = [2/3, -1/3, -2/3][species]
                for a, c, target_spin, source_spin in itertools.product(range(3), range(3), range(2), range(2)):
                    coefficient = spin[target_spin, source_spin]
                    if coefficient == 0:
                        continue
                    target = 18*i+6*species+3*target_spin+a
                    source = 18*j+6*species+3*source_spin+c
                    self.hops.append((target, source, edge, a, c, species == 2, coefficient, i, j, charge))
                    self.hops.append((source, target, edge, a, c, species != 2, coefficient.conjugate(), j, i, charge))
        self.pairs = []
        for site, colour, spin in itertools.product(range(4), range(3), range(2)):
            q = 18*site+3*spin+colour
            anti = 18*site+12+3*(1-spin)+colour
            self.pairs.append((q, anti, site, p['sea_yukawa']*(1 if spin == 0 else -1)/np.sqrt(2)))
        self.count_cache: dict[int, tuple[int, ...]] = {}

    def counts(self, bits: int) -> tuple[int, ...]:
        if bits not in self.count_cache:
            self.count_cache[bits] = tuple(((bits >> (18*s)) & ((1 << 18)-1)).bit_count() for s in range(4))
        return self.count_cache[bits]

    def quantum_channel(self, key: int) -> tuple[int, ...]:
        bits = key & FERMION_MASK
        net_u, down, antiparticles = [], [], 0
        for site in range(4):
            local = bits >> (18*site)
            u, d, anti = (local & 63).bit_count(), ((local >> 6) & 63).bit_count(), ((local >> 12) & 63).bit_count()
            net_u.append(u-anti); down.append(d); antiparticles += anti
        return tuple(net_u+down+[antiparticles])

    def cloud_change(self, key: int, sites: tuple[int, ...], matrices: list[Any]) -> list[tuple[int, complex]]:
        shifts = [CLOUD_SHIFT+s for s in sites]
        initial = [(key >> shift) & 1 for shift in shifts]
        out = []
        for final in itertools.product(range(2), repeat=len(sites)):
            amplitude = complex(np.prod([matrix[b, a] for matrix, b, a in zip(matrices, final, initial)]))
            if amplitude != 0:
                changed = key
                for shift, a, b in zip(shifts, initial, final):
                    changed ^= (a ^ b) << shift
                out.append((changed, amplitude))
        return out

    def action_key(self, key: int, kind: str = 'H', axis: int = 0) -> dict[int, complex]:
        out: dict[int, complex] = defaultdict(complex)
        bits = key & FERMION_MASK; counts = self.counts(bits)
        if kind == 'H':
            diagonal = self.p['constituent_mass']*bits.bit_count()
            diagonal += sum(self.local_energies[n][(key >> (CLOUD_SHIFT+s)) & 1] for s, n in enumerate(counts))
            diagonal += self.p['colour_electric']*4/3*sum(((key >> (FERMIONS+5*e)) & 31) != 0 for e in range(6))
            out[key] += diagonal
            gradient = self.p['deformation_gradient']
            if gradient != 0:
                for i, j in self.edges:
                    for site in [i, j]:
                        for target, value in self.cloud_change(key, (site,), [self.X2[counts[site]]]):
                            out[target] += gradient/2*value
                    for target, value in self.cloud_change(key, (i, j), [self.X[counts[i], counts[i]], self.X[counts[j], counts[j]]]):
                        out[target] -= gradient*value
            if self.p['sea_yukawa'] != 0:
                for q, anti, site, amplitude in self.pairs:
                    for create in [True, False]:
                        ops = [(True, q), (True, anti)] if create else [(False, anti), (False, q)]
                        changed = fermion(bits, ops)
                        if changed is None:
                            continue
                        newbits, sign = changed
                        target_load = counts[site]+(2 if create else -2)
                        base = (key & ~FERMION_MASK) | newbits
                        for target, value in self.cloud_change(base, (site,), [self.X[target_load, counts[site]]]):
                            out[target] += sign*amplitude*value
        if kind in ['H', 'hop', 'J', 'K']:
            for target_mode, source_mode, edge, a, b, adjoint, amplitude, target_site, source_site, charge in self.hops:
                changed = fermion(bits, [(True, target_mode), (False, source_mode)])
                if changed is None:
                    continue
                newbits, sign = changed
                shift = FERMIONS+5*edge; link = (key >> shift) & 31
                transitions = self.transitions[a, b, adjoint, link]
                if not transitions:
                    continue
                factor = amplitude*sign
                dx = self.positions[target_site, axis]-self.positions[source_site, axis]
                if kind == 'J':
                    factor *= -1j*charge*dx
                if kind == 'K':
                    factor *= -(charge*dx)**2
                if factor == 0:
                    continue
                matrices = [self.overlap[counts[target_site]+1, counts[target_site]], self.overlap[counts[source_site]-1, counts[source_site]]]
                for newlink, link_value in transitions:
                    base = ((key & ~FERMION_MASK) | newbits) ^ ((link ^ newlink) << shift)
                    for target, value in self.cloud_change(base, (target_site, source_site), matrices):
                        out[target] += factor*link_value*value
        if kind == 'H' and self.p['colour_magnetic'] != 0:
            for face in self.faces:
                for target, value in self.loop_key(key, face).items():
                    out[target] -= self.p['colour_magnetic']*value/2
                for target, value in self.loop_key(key, tuple(reversed(face))).items():
                    out[target] -= self.p['colour_magnetic']*value/2
        return {target: value for target, value in out.items() if value != 0}

    def action(self, wave: Wave, kind: str = 'H', axis: int = 0, allowed: set[int] | None = None) -> Wave:
        out: dict[int, complex] = defaultdict(complex)
        for key, coefficient in wave.items():
            for target, amplitude in self.action_key(key, kind, axis).items():
                if allowed is None or target in allowed:
                    out[target] += coefficient*amplitude
        return {key: value for key, value in out.items() if value != 0}

    def loop_key(self, key: int, vertices: tuple[int, ...]) -> Wave:
        result: dict[int, complex] = defaultdict(complex)
        for colours in itertools.product(range(3), repeat=3):
            partial = {key: 1.+0j}
            for k, i in enumerate(vertices):
                j = vertices[(k+1) % 3]; a, b = colours[k], colours[(k+1) % 3]
                edge = self.edges.index(tuple(sorted((i, j))))
                adjoint = i > j
                if adjoint:
                    a, b = b, a
                shift = FERMIONS+5*edge
                next_partial: dict[int, complex] = defaultdict(complex)
                for source, coefficient in partial.items():
                    link = (source >> shift) & 31
                    for target, amplitude in self.transitions[a, b, adjoint, link]:
                        next_partial[source ^ ((link ^ target) << shift)] += coefficient*amplitude
                partial = next_partial
            for target, amplitude in partial.items():
                result[target] += amplitude
        return dict(result)

    def loop(self, wave: Wave, vertices: tuple[int, ...] = (0, 1, 2)) -> Wave:
        result: dict[int, complex] = defaultdict(complex)
        for key, amplitude in wave.items():
            for target, coefficient in self.loop_key(key, vertices).items():
                result[target] += amplitude*coefficient
        return dict(result)

    def baryon(self, site: int, down_spin: bool = False) -> tuple[Wave, float]:
        out: dict[int, complex] = defaultdict(complex)
        for colours in itertools.permutations(range(3)):
            sign = epsilon(*colours)
            for spins, amplitude in [((0, 0, 1), np.sqrt(2/3)/np.sqrt(12)), ((0, 1, 0), -1/np.sqrt(18))]:
                if down_spin:
                    spins = tuple(1-s for s in spins)
                modes = [18*site+3*spins[0]+colours[0], 18*site+3*spins[1]+colours[1], 18*site+6+3*spins[2]+colours[2]]
                created = fermion(0, [(True, mode) for mode in modes])
                assert created is not None
                bits, car_sign = created
                out[bits] += sign*car_sign*amplitude
        return normalized(dict(out))

    def gauss(self, wave: Wave, site: int, generator: int) -> Wave:
        t = generators()[generator]
        result: dict[int, complex] = defaultdict(complex)
        for key, coefficient in wave.items():
            bits = key & FERMION_MASK
            for species, spin, a, b in itertools.product(range(3), range(2), range(3), range(3)):
                matrix = -t.T if species == 2 else t
                amplitude = matrix[a, b]
                if amplitude == 0:
                    continue
                moved = fermion(bits, [(True, 18*site+6*species+3*spin+a), (False, 18*site+6*species+3*spin+b)])
                if moved is not None:
                    target, sign = moved
                    result[(key & ~FERMION_MASK) | target] += coefficient*amplitude*sign
            for edge, (i, j) in enumerate(self.edges):
                if site not in [i, j]:
                    continue
                matrix = self.left[generator] if site == i else self.right[generator]
                shift = FERMIONS+5*edge; initial = (key >> shift) & 31
                for target in np.flatnonzero(matrix[:, initial]):
                    result[key ^ ((initial ^ int(target)) << shift)] += coefficient*matrix[target, initial]
        return {key: value for key, value in result.items() if value != 0}
