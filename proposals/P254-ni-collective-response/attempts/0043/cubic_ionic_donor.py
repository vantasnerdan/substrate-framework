"""Third derivative of the finite ordinary 3D Ewald BvK donor, Ry/bohr^3.

Fixed alpha=1 bohr^-2, real radius=12 bohr, reciprocal radius=12 bohr^-1.
The supercell reciprocal set is spherical in physical K. At primitive-periodic
coordinates its scalar equals N times the primitive finite scalar. Derivatives
hold the finite branch fixed; a real-cutoff contact is rejected. Self/background
constants have zero third derivative. This is NOT the old d2ionq Hessian.

No dense cell/q tensor: each call produces only selected external rows of one
Cartesian three-leg block. Lattice gauge u_R=N^-1/2 sum_q exp(+iq.R) u_q.
"""
from collections import OrderedDict
from itertools import product
from math import pi, sqrt
import numpy as np
from scipy.special import erfc


def radial_third(d, alpha=1.):
    """D^3[erfc(sqrt(alpha)*r)/r], before the Ry charge factor."""
    r = np.linalg.norm(d)
    if r == 0:
        raise ValueError('singular radial derivative')
    e = np.exp(-alpha*r*r); c = 2*sqrt(alpha/pi); f = erfc(sqrt(alpha)*r)
    first = -f/r**2-c*e/r
    second = 2*f/r**3+2*c*e/r**2+2*alpha*c*e
    third = -6*f/r**4-6*c*e/r**3-4*alpha*c*e/r-4*alpha**2*c*r*e
    a = third/r**3-3*second/r**4+3*first/r**5
    b = second/r**2-first/r**3
    eye = np.eye(3)
    return (a*np.einsum('a,b,c->abc',d,d,d)
            +b*(np.einsum('ab,c->abc',eye,d)+np.einsum('ac,b->abc',eye,d)
                +np.einsum('bc,a->abc',eye,d)))


class IonicThird:
    def __init__(self, cell, positions, charges, mesh=(4,4,4), g_block=512):
        self.cell = np.asarray(cell); self.positions = np.asarray(positions)
        self.charges = np.asarray(charges); self.mesh = np.asarray(mesh,int)
        self.n = int(np.prod(mesh)); self.nat = len(charges)
        if self.positions.shape != (self.nat,3) or self.cell.shape != (3,3):
            raise ValueError('cell columns and original atom rows required')
        self.reciprocal = 2*pi*np.linalg.inv(self.cell).T
        self.volume = abs(np.linalg.det(self.cell)); self.g_block = g_block
        self.real_pairs = tuple(self._real_pairs())
        self._reciprocal_cache = OrderedDict()

    def _real_pairs(self):
        inverse = np.linalg.inv(self.cell)
        bounds = np.ceil(12*np.linalg.norm(inverse,axis=1)+.5).astype(int)
        integers = np.array(list(product(*(range(-b,b+1) for b in bounds))))
        guard = 64*np.finfo(float).eps*max(12.,np.linalg.norm(self.cell))
        for i in range(self.nat):
            for j in range(self.nat):
                delta = self.positions[i]-self.positions[j]
                shifts = integers+np.rint(inverse@delta).astype(int)
                translations = shifts@self.cell.T
                vectors = delta-translations
                radii = np.linalg.norm(vectors,axis=1)
                if np.any(np.abs(radii-12.) <= guard):
                    raise ValueError('finite real support is not differentiable at cutoff')
                for lattice,d,r in zip(shifts,vectors,radii):
                    if i == j and not np.any(lattice):
                        continue
                    if r <= guard:
                        raise ValueError('coincident periodic ions')
                    if r < 12.:
                        # 1/2 ordered pairs times QE e2=2.
                        yield i,j,lattice,self.charges[i]*self.charges[j]*radial_third(d)

    def _shifted_reciprocal(self, shift):
        # K=G-shift. Bounds include the entire shifted physical sphere.
        center = np.linalg.solve(self.reciprocal,shift)
        bounds = np.ceil(12*np.linalg.norm(self.cell,axis=0)/(2*pi)+np.abs(center)).astype(int)
        integers = np.array(list(product(*(range(-b,b+1) for b in bounds))))
        vectors = integers@self.reciprocal.T-shift
        squared = np.einsum('ga,ga->g',vectors,vectors)
        return vectors[(squared>1e-24)&(squared<=144.)]

    def _reciprocal_pairs(self,shift_integer,selected):
        """Independent pair orientations; bounded512-shift small intermediate LRU.

        Exact integer subset shifts are retained, not folded/rounded. All signed
        vectors from the original physical sphere participate. No q/Phi bank is
        stored; only selected-external-atom pair tensors survive this call.
        """
        key = (tuple(shift_integer),selected)
        if key in self._reciprocal_cache:
            self._reciprocal_cache.move_to_end(key)
            return self._reciprocal_cache[key]
        shift = self.reciprocal@(np.asarray(shift_integer)/self.mesh)
        k = self._shifted_reciprocal(shift)
        squared = np.einsum('ga,ga->g',k,k)
        radial = np.exp(-squared/4)/squared
        phase = np.exp(1j*k@self.positions.T)
        reverse_phase = phase.conj()
        chosen = phase[:,selected]; chosen_reverse = reverse_phase[:,selected]
        forward = np.empty((len(selected),self.nat,3,3,3),complex)
        reverse = np.empty((self.nat,len(selected),3,3,3),complex)
        for a,b,c in product(range(3),repeat=3):
            weight = -1j*radial*k[:,a]*k[:,b]*k[:,c]
            forward[:,:,a,b,c] = (chosen*weight[:,None]).T@reverse_phase
            reverse[:,:,a,b,c] = phase.T@(chosen_reverse*weight[:,None])
        self._reciprocal_cache[key] = (forward,reverse)
        if len(self._reciprocal_cache)>512:
            self._reciprocal_cache.popitem(last=False)
        return forward,reverse


    def block(self, q_indices, external=None, piece='total'):
        """q_indices[3,3] integer mesh labels, summing to zero MOD mesh.

        external selects atom-major Cartesian rows on FIRST leg. Remaining
        legs retain every Cartesian mode. No atomic-basis exp(iq.tau) gauge.
        """
        q = np.asarray(q_indices,int)
        if q.shape != (3,3) or np.any(np.sum(q,axis=0)%self.mesh):
            raise ValueError('exact three-leg mesh conservation required')
        d = 3*self.nat
        rows = np.arange(d) if external is None else np.asarray(external,int)
        if np.any(rows<0) or np.any(rows>=d) or len(set(rows.tolist())) != len(rows):
            raise ValueError('distinct valid Cartesian external rows required')
        rowmap = {int(row): index for index,row in enumerate(rows)}
        selected = tuple(sorted({int(row)//3 for row in rows}))
        selected_set = set(selected)
        out = np.zeros((len(rows),d,d),complex)
        if piece in ('total','real'):
            for i,j,lattice,tensor in self.real_pairs:
                if i not in selected_set and j not in selected_set:
                    continue
                phases = np.exp(2j*pi*(q/self.mesh)@lattice)
                for mask in range(8):
                    atoms = [i if mask & (1<<leg) else j for leg in range(3)]
                    if atoms[0] not in selected_set:
                        continue
                    coefficient = np.prod([1 if mask & (1<<leg) else -phases[leg] for leg in range(3)])
                    self._accumulate(out,rowmap,atoms,coefficient*tensor/sqrt(self.n))
        if piece in ('total','reciprocal'):
            for mask in range(8):
                subset = [leg for leg in range(3) if mask & (1<<leg)]
                shift_integer = np.sum(q[subset],axis=0) if subset else np.zeros(3,int)
                forward,reverse = self._reciprocal_pairs(shift_integer,selected)
                sign = (-1)**(3-len(subset))
                first = selected if mask & 1 else range(self.nat)
                second = range(self.nat) if mask & 1 else selected
                for ii,i in enumerate(first):
                    for jj,j in enumerate(second):
                        atoms = [i if mask & (1<<leg) else j for leg in range(3)]
                        tensor = forward[ii,j] if mask & 1 else reverse[i,jj]
                        coefficient = sign*4*pi/self.volume/sqrt(self.n)*self.charges[i]*self.charges[j]
                        self._accumulate(out,rowmap,atoms,coefficient*tensor)
        if piece not in ('total','real','reciprocal'):
            raise ValueError('unknown donor contribution')
        return out

    @staticmethod
    def _accumulate(out,rowmap,atoms,tensor):
        for a in range(3):
            row = rowmap.get(3*atoms[0]+a)
            if row is not None:
                out[row,3*atoms[1]:3*atoms[1]+3,3*atoms[2]:3*atoms[2]+3] += tensor[a]
