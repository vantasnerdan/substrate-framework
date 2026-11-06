"""Differentiable ordinary-NC uniform-field operator, distinct from saved H0.

The reference is the original UPF/QE radial Simpson sum with exact spherical
Bessel functions, not QE's forward cubic reciprocal interpolation. For s=p.p,
B_lm(p)=(-i)^l S_lm(p) F_l(s), S_lm=|p|^l Y_lm(p),
F_l(s)=4*pi/sqrt(Omega) sum_r w_r r beta_r(r) r^l j_l(r sqrt(s))/(r sqrt(s))^l.
The ratio is entire in s. No angular division, origin cutoff or fitted term
occurs. F_l^(n)=pref sum w r beta_r (-1/2)^n r^(l+2n)
                     j_(l+n)(r sqrt(s))/(r sqrt(s))^(l+n).
The optional quintic Hermite representation of F(s) is C2, including knots;
it matches F,F',F'' at nodes, and is a declared approximation, not the UPF sum.

field_action(..., anchored=False) returns the quadrature/Hermite Peierls family. The
explicit anchored family is H_A25(0)+H_UPF(a)-H_UPF(0), preserving A25 H0 but
not claiming its constant correction is Peierls dressed. Both have the same
analytic L and C. Measure field_action(0, anchored=False)-H_A25(0) before using this
reference approximation. Shift a=e A/(hbar c) in bohr^-1, positive e; current
in electronic atomic units is -L/2, energy is H0+a L+a^2 C/2 in Ry.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.interpolate import BPoly
from scipy.linalg.blas import zgemm
from scipy.special import hyp0f1, spherical_jn

from .qe_material_current import MaterialOpticalCurrent
from .qe_nc_projectors import place_projectors
from .qe_nc_radial import qe_simpson_weights


def _direction(value):
    value = np.asarray(value, dtype=float)
    if value.shape != (3,) or not np.isfinite(value).all():
        raise ValueError('finite Cartesian direction/shift (3,) required')
    return value


def regular_bessel(l, x):
    """j_l(x)/x^l, with its analytic value at zero, for l=0..4.

    Two stable evaluations of the SAME entire function, not two fitted
    branches: 0F1 for |x|<=1 and scipy spherical_jn otherwise.
    """
    x = np.asarray(x, dtype=float)
    if l not in range(5) or np.any(x < 0) or not np.isfinite(x).all():
        raise ValueError('l=0..4 and finite nonnegative arguments required')
    result = np.empty_like(x)
    small = x <= 1
    result[small] = hyp0f1(l+1.5, -x[small]**2/4) / math.prod(range(1, 2*l+2, 2))
    result[~small] = spherical_jn(l, x[~small])/x[~small]**l
    return result


@dataclass(frozen=True)
class UPFRadial:
    """Exact analytic function of the retained original radial quadrature."""
    tables: object
    row_block: int = 256

    def __post_init__(self):
        if self.row_block < 1 or np.any(~np.isin(self.tables.upf.angular_l, (0, 1, 2))):
            raise ValueError('positive row block and ordinary s/p/d UPF required')

    def evaluate(self, s):
        """Return F,F_s,F_ss, each (npoint,nbeta), not radial q derivatives."""
        s = np.asarray(s, dtype=float)
        if s.ndim != 1 or np.any(s < 0) or not np.isfinite(s).all():
            raise ValueError('finite nonnegative squared moduli required')
        upf = self.tables.upf
        r = upf.r[:upf.beta_mesh]
        weight = (qe_simpson_weights(upf.rab[:upf.beta_mesh])*r)[:, None]*upf.beta_r[:upf.beta_mesh]
        pref = 4*np.pi/np.sqrt(self.tables.omega_bohr3)
        result = tuple(np.empty((s.size, upf.nbeta)) for _ in range(3))
        for start in range(0, s.size, self.row_block):
            stop = min(start+self.row_block, s.size)
            x = np.sqrt(s[start:stop, None])*r
            for ell in np.unique(upf.angular_l):
                ell = int(ell)
                cols = np.flatnonzero(upf.angular_l == ell)
                for order in range(3):
                    result[order][start:stop, cols] = (
                        regular_bessel(ell+order, x)
                        @ (weight[:, cols]*r[:, None]**(ell+2*order))
                        * (pref*(-0.5)**order))
        return result


class HermiteRadial:
    """Explicit C2 reciprocal representation; refinement is numerical, not physical.

    Nodes uniform in q with zero included, interpolating F in s=q^2. Values,
    first and second s derivatives come from UPFRadial, not finite stencils.
    No extrapolation is allowed. Increasing qmax does not alter existing
    intervals; decreasing spacing is the controlled representation study.
    """
    def __init__(self, exact, qmax_bohr_inv, spacing_bohr_inv):
        qmax, spacing = float(qmax_bohr_inv), float(spacing_bohr_inv)
        if not np.isfinite([qmax, spacing]).all() or qmax <= 0 or spacing <= 0:
            raise ValueError('positive finite qmax and reciprocal spacing required')
        self.tables = exact.tables
        self.spacing_bohr_inv = spacing
        q = np.arange(int(np.ceil(qmax/spacing))+1)*spacing
        self.smax = q[-1]**2
        values = exact.evaluate(q*q)
        self.polynomial = BPoly.from_derivatives(q*q, np.stack(values, axis=1), extrapolate=False)
        self.derivatives = (self.polynomial, self.polynomial.derivative(), self.polynomial.derivative(2))

    def evaluate(self, s):
        s = np.asarray(s, dtype=float)
        if s.ndim != 1 or not np.isfinite(s).all() or np.any(s < 0) or np.any(s > self.smax):
            raise ValueError('s outside finite Hermite radial representation')
        return tuple(poly(s) for poly in self.derivatives)


def solid_harmonics(momenta, left, right):
    """QE-ordered regular s/p/d polynomials and exact directional derivatives."""
    p = np.asarray(momenta, dtype=float)
    e, f = _direction(left), _direction(right)
    if p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all():
        raise ValueError('finite momenta (npw,3) required')
    # S = constant + linear.p + p.T Q p; Q is a real polynomial coefficient,
    # not a numerical Hermitization of any physical matrix.
    constant = np.zeros(9)
    constant[0] = 1/np.sqrt(4*np.pi)
    linear = np.zeros((9, 3))
    c = np.sqrt(3/(4*np.pi))
    linear[1, 2], linear[2, 0], linear[3, 1] = c, -c, -c
    quadratic = np.zeros((9, 3, 3))
    c = np.sqrt(5/(4*np.pi))
    quadratic[4] = np.diag([-c/2, -c/2, c])
    quadratic[5, 0, 2] = quadratic[5, 2, 0] = -c*np.sqrt(3)/2
    quadratic[6, 1, 2] = quadratic[6, 2, 1] = -c*np.sqrt(3)/2
    quadratic[7] = np.diag([c*np.sqrt(3)/2, -c*np.sqrt(3)/2, 0])
    quadratic[8, 0, 1] = quadratic[8, 1, 0] = c*np.sqrt(3)/2
    value = constant + p@linear.T + np.einsum('ni,kij,nj->nk', p, quadratic, p)
    de = linear@e + 2*np.einsum('ni,kij,j->nk', p, quadratic, e)
    df = linear@f + 2*np.einsum('ni,kij,j->nk', p, quadratic, f)
    mixed = np.broadcast_to(2*np.einsum('i,kij,j->k', e, quadratic, f), value.shape)
    return value, de, df, mixed


def centered_derivatives(radial, momenta, left, right=None):
    """Return B,B_e,B_f,B_ef; all angular/radial derivative legs included."""
    p = np.asarray(momenta, dtype=float)
    e = _direction(left)
    f = e if right is None else _direction(right)
    solid, se, sf, sef = solid_harmonics(p, e, f)
    upf = radial.tables.upf
    r, rs, rss = radial.evaluate(np.einsum('ni,ni->n', p, p))
    r, rs, rss = (a[:, upf.channel_radial] for a in (r, rs, rss))
    solid, se, sf, sef = (a[:, upf.channel_lm] for a in (solid, se, sf, sef))
    ds_e, ds_f = 2*(p@e)[:, None], 2*(p@f)[:, None]
    phase = ((-1j)**upf.channel_l)[None, :]
    return tuple(a*phase for a in (
        solid*r,
        se*r+solid*rs*ds_e,
        sf*r+solid*rs*ds_f,
        sef*r+se*rs*ds_f+sf*rs*ds_e+solid*(rss*ds_e*ds_f+2*(e@f)*rs)))


def product_actions(beta, left, right, mixed, active_d, vectors):
    """L_e and C_ef actions, with both legs and both mixed cross terms."""
    overlap = lambda b: active_d @ zgemm(1., b, vectors, trans_a=2)
    b, e, f, ef = (overlap(a) for a in (beta, left, right, mixed))
    linear = zgemm(1., left, b)+zgemm(1., beta, e)
    contact = (zgemm(1., mixed, b)+zgemm(1., left, f)
               +zgemm(1., right, e)+zgemm(1., beta, ef))
    return linear, contact


class MaterialFieldDerivatives:
    """Full-row selected-column actions, no dense PW matrix or new DFT state."""
    def __init__(self, hamiltonian, *, spacing_bohr_inv=None, qmax_bohr_inv=None):
        self.hamiltonian = hamiltonian
        self.original = MaterialOpticalCurrent(hamiltonian)
        self.momenta = self.original.momenta
        self.radial = {}
        for species in hamiltonian.factors.species:
            exact = UPFRadial(hamiltonian.background.inputs.tables[species.name])
            if spacing_bohr_inv is None:
                self.radial[species.name] = exact
            else:
                if qmax_bohr_inv is None:
                    raise ValueError('explicit qmax required for Hermite representation')
                self.radial[species.name] = HermiteRadial(exact, qmax_bohr_inv, spacing_bohr_inv)

    def derivative_actions(self, vectors, left, right=None):
        """Return analytic (L_left,C_left,right); no finite difference parameter."""
        vectors, single = self.original._vectors(vectors)
        e = _direction(left)
        f = e if right is None else _direction(right)
        linear = np.asarray(2*(self.momenta@e)[:, None]*vectors, dtype=complex)
        contact = np.asarray(2*(e@f)*vectors, dtype=complex)
        inputs = self.hamiltonian.background.inputs
        for species in self.hamiltonian.factors.species:
            centered = centered_derivatives(self.radial[species.name], self.momenta, e, f)
            for atom in species.atom_indices:
                placed = tuple(place_projectors(b, self.momenta, inputs.reference.positions_bohr[:, atom-1]) for b in centered)
                for start in range(0, vectors.shape[1], self.hamiltonian.block_size):
                    stop = min(start+self.hamiltonian.block_size, vectors.shape[1])
                    l, c = product_actions(*placed, species.active_d, vectors[:, start:stop])
                    linear[:, start:stop] += l
                    contact[:, start:stop] += c
        return (linear[:, 0], contact[:, 0]) if single else (linear, contact)

    def nonlocal_delta_action(self, vectors, shift_bohr_inv, *, anchored):
        """V_rep(a)-V_rep(0) if anchored, else V_rep(a)-V_A25(0).

        Exposes the small difference without subtracting large local-H actions.
        """
        shift = _direction(shift_bohr_inv)
        vectors, single = self.original._vectors(vectors)
        out = np.zeros_like(vectors, dtype=complex)
        if anchored and not np.any(shift):
            return out[:, 0] if single else out
        inputs = self.hamiltonian.background.inputs
        for species in self.hamiltonian.factors.species:
            radial = self.radial[species.name]
            shifted = centered_derivatives(radial, self.momenta+shift, np.zeros(3))[0]
            center = centered_derivatives(radial, self.momenta, np.zeros(3))[0] if anchored else None
            nh = species.active_d.shape[0]
            for slot, atom in enumerate(species.atom_indices):
                tau = inputs.reference.positions_bohr[:, atom-1]
                beta = place_projectors(shifted, self.momenta, tau)
                original = (place_projectors(center, self.momenta, tau) if anchored
                            else species.beta[:, slot*nh:(slot+1)*nh])
                for start in range(0, vectors.shape[1], self.hamiltonian.block_size):
                    stop = min(start+self.hamiltonian.block_size, vectors.shape[1])
                    v = vectors[:, start:stop]
                    out[:, start:stop] += (zgemm(1., beta, species.active_d @ zgemm(1., beta, v, trans_a=2))
                                          -zgemm(1., original, species.active_d @ zgemm(1., original, v, trans_a=2)))
        return out[:, 0] if single else out

    def field_increment_action(self, vectors, shift_bohr_inv):
        """H_rep(a)-H_rep(0), avoiding local-H cancellation."""
        shift = _direction(shift_bohr_inv)
        vectors, single = self.original._vectors(vectors)
        out = self.nonlocal_delta_action(vectors, shift, anchored=True)
        out += (2*(self.momenta@shift)+shift@shift)[:, None]*vectors
        return out[:, 0] if single else out

    def field_action(self, vectors, shift_bohr_inv, *, anchored=False):
        """Physical UPF family by default; the A25 offset is an explicit approximation."""
        shift = _direction(shift_bohr_inv)
        vectors, single = self.original._vectors(vectors)
        out = self.hamiltonian.apply(vectors)
        out += self.nonlocal_delta_action(vectors, shift, anchored=anchored)
        out += (2*(self.momenta@shift)+shift@shift)[:, None]*vectors
        return out[:, 0] if single else out
