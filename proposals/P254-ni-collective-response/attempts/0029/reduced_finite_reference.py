"""Small-lattice oracle only: independently assembled dense shifted scalar.

This intentionally does not call dynamics block maps, densities, Hamiltonians,
force maps, energies, or packers. It uses forward physical operators and traces;
reverse and density export coefficients are therefore independently challenged.
"""
import numpy as np


class DenseFiniteReference:
    def __init__(self, c):
        self.c = c
        self.length = c.nk*c.m
        self.h0 = self.diagonal(c.h0)
        self.rho0 = self.diagonal(c.rho0)
        shape = (c.nspin,self.length,self.length)
        self.ionic = np.zeros((c.nq,c.d)+shape,complex)
        self.profiles = [np.zeros((c.feature_count(q),)+shape,complex) for q in range(c.nq)]
        self.partners = []
        centered = np.asarray([np.rint(np.fft.fftfreq(size)*size).astype(int)[c.q_indices[:,axis]]
                               for axis,size in enumerate(c.mesh)]).T
        for q in range(c.nq):
            inverse = int(np.ravel_multi_index(tuple((-c.q_indices[q]) % c.mesh),c.mesh))
            data,other = c.q(q),c.q(inverse)
            fold = (centered[q]+centered[inverse])//np.asarray(c.mesh)
            lookup = {(str(spin),tuple(g)):i*len(other.feature_g)+j
                      for i,spin in enumerate(c.spin_labels) for j,g in enumerate(other.feature_g)}
            self.partners.append(np.asarray([lookup[(str(spin),tuple(-g-fold))]
                for spin in c.spin_labels for g in data.feature_g],dtype=int))
        self.field = np.zeros((c.f,)+shape,complex)
        # Infer q from each dense row/column pair, not from target[q,slot].
        for spin in range(c.nspin):
            for source in range(c.nk):
                slot = spin*c.nk+source
                source_slice = slice(source*c.m,(source+1)*c.m)
                for target in range(c.nk):
                    difference = (c.k_indices[target]-c.k_indices[source]) % c.mesh
                    q = int(np.ravel_multi_index(tuple(difference),c.mesh))
                    data = c.q(q)
                    target_slice = slice(target*c.m,(target+1)*c.m)
                    self.ionic[q,:,spin,target_slice,source_slice] = data.ionic[slot]
                    self.profiles[q][:,spin,target_slice,source_slice] = np.moveaxis(data.potential[slot],-1,0)
                self.field[:,spin,source_slice,source_slice] = np.einsum('fi,iab->fab',c.field_polarizations,c.field[slot])

    def diagonal(self, blocks):
        c = self.c
        result = np.zeros((c.nspin,self.length,self.length),complex)
        for spin in range(c.nspin):
            for k in range(c.nk):
                sl = slice(k*c.m,(k+1)*c.m)
                result[spin,sl,sl] = blocks[spin*c.nk+k]
        return result

    def diagonal_blocks(self, full):
        c = self.c
        return np.asarray([full[spin,k*c.m:(k+1)*c.m,k*c.m:(k+1)*c.m]
                           for spin in range(c.nspin) for k in range(c.nk)])

    @staticmethod
    def trace(operator, density):
        return sum(np.trace(a@b) for a,b in zip(operator,density))

    def densities(self, change):
        c = self.c
        result = []
        for q in range(c.nq):
            inverse = int(np.ravel_multi_index(tuple((-c.q_indices[q]) % c.mesh),c.mesh))
            result.append(np.asarray([self.trace(self.profiles[inverse][partner],change)/np.sqrt(c.n)
                                      for partner in self.partners[q]]))
        return result

    def displacement_potential(self, x):
        c = self.c
        potential = np.zeros_like(self.h0)
        for q in range(c.nq):
            for axis in range(c.d):
                potential += self.ionic[q,axis]*x[q][axis]/np.sqrt(c.n)
        for mode in range(c.f):
            potential += self.field[mode]*x[0][2*c.d+mode]
        return potential

    def hamiltonian(self, rho, x):
        c = self.c
        density = self.densities(rho-self.rho0)
        h = self.h0+self.displacement_potential(x)
        for q in range(c.nq):
            potential = c.q(q).kernel@density[q]
            for r in range(len(potential)):
                h += self.profiles[q][r]*potential[r]/np.sqrt(c.n)
        return h

    def rhs(self, rho, x):
        c = self.c
        h = self.hamiltonian(rho,x)
        rho_dot = -1j*(h@rho-rho@h)
        change = rho-self.rho0
        x_dot = []
        for q in range(c.nq):
            gradient = np.zeros_like(x[q])
            inverse = int(np.ravel_multi_index(tuple((-c.q_indices[q]) % c.mesh),c.mesh))
            for axis in range(c.d):
                gradient[axis] = self.trace(self.ionic[inverse,axis],change)/np.sqrt(c.n)
            if q == 0:
                for mode in range(c.f):
                    gradient[2*c.d+mode] = self.trace(self.field[mode],change)
                gradient += c.mechanical.linear
            x_dot.append(c.mechanical.generators[q]@x[q]+c.mechanical.j[q]@gradient)
        return rho_dot,tuple(x_dot)

    def energy(self, rho, x):
        c = self.c
        change = rho-self.rho0
        mechanical = c.mechanical.linear@x[0]
        for q in range(c.nq):
            inverse = int(np.ravel_multi_index(tuple((-c.q_indices[q]) % c.mesh),c.mesh))
            mechanical += x[inverse]@c.mechanical.h[q]@x[q]/2
        density = self.densities(change)
        screening = 0j
        for q in range(c.nq):
            inverse = int(np.ravel_multi_index(tuple((-c.q_indices[q]) % c.mesh),c.mesh))
            screening += density[inverse][self.partners[q]]@c.q(q).kernel@density[q]/2
        return dict(mechanical=mechanical,band=self.trace(self.h0,change),
                    interaction=self.trace(self.displacement_potential(x),change),screening=screening)

    def pack(self, rho, x):
        return np.concatenate((rho.ravel(),*x))

    def unpack(self, vector):
        c = self.c
        stop = c.nspin*self.length**2
        rho = vector[:stop].reshape(c.nspin,self.length,self.length)
        x = []
        for size in c.mechanical_sizes:
            x.append(vector[stop:stop+size])
            stop += size
        return rho,tuple(x)

    def packed_rhs(self, time, vector):
        del time
        return self.pack(*self.rhs(*self.unpack(vector)))
