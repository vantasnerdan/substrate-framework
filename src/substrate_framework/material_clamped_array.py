"""Periodic ionic/q0-field CONTROL with explicitly frozen prepared electrons.

H = sum_R [p_R-Z E.T a]^2/(2M) - F.sum_R u_R
  + 1/2 sum_RS u_R.C(R-S).u_S + sum_R u_R.T.E.T.a
  + N ell.E.T.a + N/2 a.E.contact.E.T.a + Pi.Pi/(2 Mf).
C includes the fixed onsite deltaC. FFT coordinates are orthonormal; only
Gamma couples to the ONE global field, with sqrt(N), not N, cross blocks.
Raw complex coefficients survive unchanged. This is not joint fermion dynamics.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm


def canonical(size):
    eye = np.eye(size)
    return np.block([[np.zeros_like(eye), eye], [-eye, np.zeros_like(eye)]])


def _array(value, shape, name):
    out = np.asarray(value, dtype=complex)
    if out.shape != shape or not np.all(np.isfinite(out)):
        raise ValueError(f'{name}: expected finite {shape}, got {out.shape}')
    return out.copy()


def relative_defect(value, reference):
    return float(np.linalg.norm(value-reference)/max(1., np.linalg.norm(reference)))


@dataclass(frozen=True)
class ArrayEvolution:
    time_inverse_ry: float
    transitions: tuple
    affine_gamma: np.ndarray


class ClampedArray:
    """All q blocks and the explicit product preparation, never clipped/repaired.

    Block order: (u_q,p_q); Gamma alone appends (a,Pi). Ordered covariance
    W_q(t,s)=S_q(t) W0_q S_-q(s).T; there is deliberately no dagger.
    For imperfect input reciprocity the raw generator is retained and defects
    are exposed, rather than presenting repaired dynamics as the same scalar.
    """
    def __init__(self, force_lattice, masses, charges, delta_curvature, contact,
                 mixed, force, field_gradient, polarizations, field_mass,
                 ionic_ordered, field_ordered, trap_omega):
        lattice = np.asarray(force_lattice, dtype=complex)
        if lattice.ndim == 7:
            lattice = lattice.reshape(*lattice.shape[:3], 3*lattice.shape[3], -1)
        if (lattice.ndim != 5 or lattice.shape[-1] != lattice.shape[-2]
                or min(lattice.shape) < 1 or lattice.shape[-1] % 3
                or not np.all(np.isfinite(lattice))):
            raise ValueError('complete periodic lattice with all Cartesian atom coordinates required')
        self.shape, self.d = lattice.shape[:3], lattice.shape[-1]
        self.n = int(np.prod(self.shape))
        self.indices = np.indices(self.shape).reshape(3, -1).T
        self.minus = np.ravel_multi_index(tuple((-self.indices % self.shape).T), self.shape)
        self.masses = _array(masses, (self.d,), 'masses')
        self.omega = _array(trap_omega, (self.d,), 'trap frequencies')
        if np.any(self.masses.imag) or np.any(self.masses.real <= 0) or np.any(self.omega.imag) or np.any(self.omega.real <= 0):
            raise ValueError('positive real original masses and preparation frequencies required')
        z = _array(charges, (self.d//3,), 'original valence charges')
        if np.any(z.imag) or np.any(z.real <= 0):
            raise ValueError('positive real original valence charges required')
        e = np.asarray(polarizations, dtype=float)
        if e.ndim != 2 or e.shape[1] != 3 or not np.all(np.isfinite(e)) or not np.allclose(e@e.T, np.eye(len(e)), atol=1e-12, rtol=1e-12):
            raise ValueError('orthonormal field polarizations required')
        self.f = len(e)
        if not np.isfinite(field_mass) or field_mass <= 0:
            raise ValueError('positive Maxwell mass required')
        dc = _array(delta_curvature, (self.d,self.d), 'fixed deltaC')
        c = _array(contact, (3,3), 'fixed contact')
        t = _array(mixed, (self.d,3), 'fixed mixed')
        force = _array(force, (self.d,), 'actual initial force')
        ell = _array(field_gradient, (3,), 'actual field gradient')
        self.w = _array(ionic_ordered, (2*self.d,2*self.d), 'saved ionic preparation')
        self.wfield = _array(field_ordered, (2*self.f,2*self.f), 'saved field preparation')
        cq = np.fft.fftn(lattice, axes=(0,1,2)).reshape(self.n,self.d,self.d)+dc
        charge_map = (z[:,None,None]*np.eye(3)[None,:,:]).reshape(self.d,3)@e.T
        self.h, self.j, self.w0 = [], [], []
        for k in range(self.n):
            size = 2*self.d+(2*self.f if k == 0 else 0)
            h = np.zeros((size,size), complex)
            h[:self.d,:self.d] = cq[k]
            h[self.d:2*self.d,self.d:2*self.d] = np.diag(1/self.masses)
            j = np.zeros_like(h)
            j[:2*self.d,:2*self.d] = canonical(self.d)
            w = np.zeros_like(h)
            w[:2*self.d,:2*self.d] = self.w
            if k == 0:
                a = slice(2*self.d,2*self.d+self.f)
                pi = slice(2*self.d+self.f,size)
                h[:self.d,a] = np.sqrt(self.n)*t@e.T
                h[a,:self.d] = h[:self.d,a].T
                h[self.d:2*self.d,a] = -np.sqrt(self.n)*charge_map/self.masses[:,None]
                h[a,self.d:2*self.d] = h[self.d:2*self.d,a].T
                h[a,a] = self.n*(e@c@e.T+charge_map.T@(charge_map/self.masses[:,None]))
                h[pi,pi] = np.eye(self.f)/field_mass
                j[2*self.d:,2*self.d:] = canonical(self.f)
                w[2*self.d:,2*self.d:] = self.wfield
            self.h.append(h)
            self.j.append(j)
            self.w0.append(w)
        self.linear = np.zeros(len(self.h[0]), complex)
        self.linear[:self.d] = -np.sqrt(self.n)*force
        self.linear[2*self.d:2*self.d+self.f] = self.n*e@ell
        self.generators = tuple(j@h for j,h in zip(self.j,self.h))

    def coefficient_diagnostics(self):
        return dict(
            paired_transpose_defect=max(relative_defect(h,self.h[self.minus[k]].T) for k,h in enumerate(self.h)),
            paired_reality_defect=max(relative_defect(h,self.h[self.minus[k]].conj()) for k,h in enumerate(self.h)),
            hermiticity_defect=max(relative_defect(h,h.conj().T) for h in self.h),
            preparation_hermiticity_defect=max(relative_defect(w,w.conj().T) for w in self.w0),
            preparation_ccr_defect=max(relative_defect(w-w.T,1j*j) for w,j in zip(self.w0,self.j)),
            electrons='held at original prepared full-column state',
            repair_applied=False)

    def evolve(self, time_inverse_ry):
        time = float(time_inverse_ry)
        if not np.isfinite(time):
            raise ValueError('finite explicit time in inverse Ry required')
        s = tuple(expm(time*a) for a in self.generators)
        a = self.generators[0]
        augmented = np.zeros((len(a)+1,len(a)+1), complex)
        augmented[:-1,:-1] = a
        augmented[:-1,-1] = self.j[0]@self.linear
        drift = expm(time*augmented)[:-1,-1]
        if not all(np.all(np.isfinite(v)) for v in (*s,drift)):
            raise FloatingPointError('nonfinite exact quadratic evolution; no instability clipping')
        return ArrayEvolution(time,s,drift)

    def kicks(self, cell, coordinate):
        """Coherent momentum displacement with added energy one local trap quantum.

        Local delta p=sqrt(2 M omega). Common delta p/sqrt(N) at every cell
        has the SAME total initial energy, not the same per-cell amplitude.
        """
        cell = self.cell(cell)
        if not 0 <= coordinate < self.d:
            raise ValueError('original zero-based Cartesian coordinate out of range')
        amplitude = np.sqrt(2*self.masses[coordinate]*self.omega[coordinate])
        local, common = [], []
        phase = np.exp(-2j*np.pi*(self.indices/np.asarray(self.shape))@cell)
        for k,h in enumerate(self.h):
            l, c = np.zeros(len(h),complex), np.zeros(len(h),complex)
            l[self.d+coordinate] = amplitude*phase[k]/np.sqrt(self.n)
            if k == 0:
                c[self.d+coordinate] = amplitude
            local.append(l)
            common.append(c)
        return tuple(local), tuple(common)

    def cell(self, cell):
        raw = np.asarray(cell)
        if raw.shape != (3,) or not np.issubdtype(raw.dtype,np.integer):
            raise ValueError('three integer original cell indices required')
        if np.any(raw < 0) or np.any(raw >= self.shape):
            raise ValueError('cell index outside original periodic grid')
        return raw

    def mean(self, evolution, initial=None):
        q = [np.zeros(len(h),complex) for h in self.h] if initial is None else [s@x for s,x in zip(evolution.transitions,initial)]
        q[0] = q[0]+evolution.affine_gamma
        ions = np.fft.ifftn(np.asarray([x[:2*self.d] for x in q]).reshape(*self.shape,2*self.d), axes=(0,1,2), norm='ortho')
        return dict(ions=ions, field=q[0][2*self.d:], q=tuple(q))

    def ordered_pair(self, left, right, cell, other, coordinates):
        """Selected ionic (u,p) entries, connected <x_cell(t) x_other(s)>.

        coordinates index the full cell phase space [u(0:d),p(0:d)].
        All other atoms/cells still participate in every transition.
        """
        delta = self.cell(cell)-self.cell(other)
        ix = np.asarray(coordinates,dtype=int)
        if ix.ndim != 1 or np.any(ix < 0) or np.any(ix >= 2*self.d):
            raise ValueError('phase-space readout indices out of range')
        phase = np.exp(2j*np.pi*(self.indices/np.asarray(self.shape))@delta)/self.n
        out = np.zeros((len(ix),len(ix)),complex)
        for k in range(self.n):
            out += phase[k]*(left.transitions[k][ix]@self.w0[k]@right.transitions[self.minus[k]][ix].T)
        return out

    def field_ordered_pair(self, left, right):
        sl = slice(2*self.d,None)
        return left.transitions[0][sl]@self.w0[0]@right.transitions[0][sl].T

    def ion_field_ordered_pair(self, left, right, coordinates):
        ix = np.asarray(coordinates,dtype=int)
        return left.transitions[0][ix]@self.w0[0]@right.transitions[0][2*self.d:].T/np.sqrt(self.n)

    def diagnostics(self, evolution, mean_q):
        covariance_energy, mean_energy, ccr, symplectic, hermiticity = 0j, 0j, 0., 0., 0.
        for k,(h,j,s) in enumerate(zip(self.h,self.j,evolution.transitions)):
            other = self.minus[k]
            sm = evolution.transitions[other]
            w = s@self.w0[k]@sm.T
            wm = sm@self.w0[other]@s.T
            ccr = max(ccr,relative_defect(w-wm.T,1j*j))
            hermiticity = max(hermiticity,relative_defect(w,w.conj().T))
            symplectic = max(symplectic,relative_defect(s@j@sm.T,j))
            covariance_energy += np.einsum('ij,ij->',h,wm)/2
            mean_energy += mean_q[other]@h@mean_q[k]/2
        mean_energy += self.linear@mean_q[0]
        return dict(mean_energy_ry=mean_energy,covariance_energy_ry=covariance_energy,
                    energy_ry=mean_energy+covariance_energy,
                    energy_zero='constant frozen electronic scalar omitted, never inferred',
                    ordered_ccr_relative_defect=ccr,paired_symplectic_relative_defect=symplectic,
                    covariance_hermiticity_relative_defect=hermiticity,
                    mean_q_reality_relative_defect=max(relative_defect(x,mean_q[self.minus[k]].conj()) for k,x in enumerate(mean_q)))
