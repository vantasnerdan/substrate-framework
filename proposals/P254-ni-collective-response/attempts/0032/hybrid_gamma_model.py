"""Moving source0 finite-bath replacement of the original Gamma Hamiltonian.

Mean and fluctuation mechanical coordinates are the SAME original A29 canonical
coordinates. Ionic insertions/forces carry 1/sqrt(N); the two Hartree insertions
carry 1/N; field insertions/forces retain the original unscaled normalization.
Only the source0 density uses an interaction-picture similarity. No occupation,
Hamiltonian, density or independently supplied adjoint is repaired here.
"""
from __future__ import annotations
import numpy as np


def _transpose(value):
    return value.swapaxes(-1,-2)


def _maximum(value):
    return float(np.max(np.abs(value),initial=0))


def _commutator(h,rho):
    return -1j*(h@rho-rho@h)


def _commutator_duals(h,rho,z):
    return (-1j*(_transpose(h)@z-z@_transpose(h)),
            -1j*(z@_transpose(rho)-_transpose(rho)@z))


class HybridGammaModel:
    """Matrix-free exact tangent and bilinear transpose of the moving mean.

    Flat layout: original mechanical x, all old sectors EXCEPT source_slot,
    full source_slot interaction-picture density. Matrix inputs batch probes in
    columns. All duals pair by sum(z*value), never by vdot or a dagger.
    """
    def __init__(self,coefficients,initial_mechanical,bath_bank,mode_directions,*,source_slot=0):
        self.c=coefficients
        c=coefficients
        if not isinstance(source_slot,(int,np.integer)) or not 0<=source_slot<c.ns:
            raise ValueError('source_slot must index an original explicit source')
        self.source_slot=int(source_slot)
        self.other_slots=np.delete(np.arange(c.ns),self.source_slot)
        self.data=c.q(0)
        self.nx=int(c.mechanical_sizes[0])
        self.root_n=float(np.sqrt(c.n))
        self.q=np.asarray(bath_bank['projected_active_coordinates'])
        self.h0=np.asarray(bath_bank['hamiltonian'])
        self.b=self.h0.shape[0]
        self.frame=np.asarray(bath_bank['projected_frame'])
        self.gram=np.asarray(bath_bank['gram'])
        self.density_coefficients=np.asarray(bath_bank['density_coefficients'])
        self.rho_b0=np.asarray(bath_bank['projected_density'])
        self.directions=np.asarray(mode_directions)
        self._require(self.h0,(self.b,self.b),'literal H')
        self._require(self.q,(self.b,c.m),'saved active embedding')
        self._require(self.rho_b0,(self.b,self.b),'original projected density')
        rank=self.frame.shape[1]
        self._require(self.frame,(self.b,rank),'original projected frame')
        self._require(self.gram,(rank,rank),'original Gram')
        self._require(self.density_coefficients,(rank,rank),'original raw D')
        self._require(self.directions,(c.d,2),'two physical force directions')
        x=np.asarray(initial_mechanical,dtype=complex)
        self._require(x,(self.nx,),'original canonical initial mechanical state')
        # g_i is the physical per-cell derivative, before canonical 1/sqrtN.
        old_force=-np.einsum('im,iab->mab',self.directions,self.data.ionic[self.source_slot])
        self.raw_forces=np.asarray([bath_bank['force_mode0'],bath_bank['force_mode1']])
        self.raw_force_adjoints=np.asarray([bath_bank['adjoint_force_mode0'],bath_bank['adjoint_force_mode1']])
        self._require(self.raw_forces,(2,self.b,self.b),'physical forward force legs')
        self._require(self.raw_force_adjoints,(2,self.b,self.b),'physical adjoint force legs')
        self.delta_force=self.raw_forces-self.q@old_force@self.q.conj().T
        self.delta_force_adjoint=self.raw_force_adjoints-self.q@_transpose(old_force).conj()@self.q.conj().T
        # Diagonalize ONLY literal H. General similarity retains all raw defects.
        self.energies,self.v=np.linalg.eig(self.h0)
        self.vi=np.linalg.inv(self.v)
        self._left0=self.q.conj().T@self.v
        self._right0=self.vi@self.q
        self._delta_eigen=np.concatenate((self.vi@self.delta_force@self.v,
                                          self.vi@self.delta_force_adjoint@self.v))
        self._initial_force_traces=np.einsum('kab,ba->k',
            np.concatenate((self.delta_force,self.delta_force_adjoint)),self.rho_b0)
        self._other_size=len(self.other_slots)*c.m*c.m
        self.size=self.nx+self._other_size+self.b*self.b
        rb=self.vi@self.rho_b0@self.v
        self.initial=self._join(x[:,None],c.rho0[self.other_slots][None],rb[None],1)
        active=c.rho0.copy()
        active[self.source_slot]=self.q.conj().T@self.rho_b0@self.q
        self._initial_active_h=self._active_h(x[:,None],(active-c.rho0)[None])[0,self.source_slot]
        ha=np.asarray(bath_bank['hamiltonian_adjoint'])
        self._require(ha,(self.b,self.b),'independent literal H adjoint')
        self.metadata=dict(schema='P254-moving-hybrid-Gamma-model-v1',source_slot=self.source_slot,
            source_count=c.ns,other_active_sectors=len(self.other_slots),active_dimension=c.m,
            source0_dimension=self.b,original_raw_rank=rank,mechanical_size=self.nx,state_size=self.size,
            state_layout='original canonical mechanical; other source slots in ascending order; source0 full interaction density',
            mean_coordinates='unchanged original A29 canonical mechanical/field coordinates; no per-cell rescaling',
            electronic_measure='ionic forward/reverse 1/sqrtN; Hartree density and potential each 1/sqrtN; field original unscaled',
            extra_hamiltonian='-1/(2 sqrtN) sum_m [(conj(U_m).T x) DeltaF_m + (U_m.T x) DeltaFadj_m]',
            quadrature='explicit physical real-displacement Hermitian-quadrature definition with BOTH independent force legs; complex tangent extension linear in x',
            interaction_picture='general eigensimilarity of raw literal H only; full noncommuting density retained',
            initial_covariance='original mechanical/field W0; other active Wick; source0 original frame/S/D CAR including finite-basis complement; no extra Fourier/cell/spin factor',
            scope='one Gamma/source0 two-mode finite-bath enrichment; not allq/fullsource/fullPW calibration, physical loss or lifetime')
        self.diagnostics=dict(
            hamiltonian_hermiticity_frobenius=float(np.linalg.norm(self.h0-self.h0.conj().T)),
            hamiltonian_supplied_adjoint_defect_frobenius=float(np.linalg.norm(ha-self.h0.conj().T)),
            eigenvector_condition=float(np.linalg.cond(self.v)),
            eigen_reconstruction_frobenius=float(np.linalg.norm((self.v*self.energies)@self.vi-self.h0)),
            eigen_inverse_defect_max_abs=_maximum(self.vi@self.v-np.eye(self.b)),
            eigen_energy_imag_max_abs=_maximum(self.energies.imag),
            force_forward_hermiticity_frobenius=[float(np.linalg.norm(f-f.conj().T)) for f in self.raw_forces],
            force_supplied_adjoint_defect_frobenius=[float(np.linalg.norm(a-f.conj().T)) for a,f in zip(self.raw_force_adjoints,self.raw_forces)],
            force_direction_imag_max_abs=_maximum(self.directions.imag),
            extra_force_adjoint_defect_frobenius=[float(np.linalg.norm(a-f.conj().T)) for a,f in zip(self.delta_force_adjoint,self.delta_force)],
            original_gram_frame_defect_frobenius=float(np.linalg.norm(self.frame.conj().T@self.frame-self.gram)),
            original_D_hermiticity_frobenius=float(np.linalg.norm(self.density_coefficients-self.density_coefficients.conj().T)),
            projected_state_transfer_frobenius=float(np.linalg.norm(self.rho_b0-self.frame@self.density_coefficients@self.frame.conj().T)),
            active_embedding_metric_defect_max_abs=_maximum(self.q.conj().T@self.q-np.eye(c.m)),
            initial_source_compression_difference_frobenius=float(np.linalg.norm(active[self.source_slot]-c.rho0[self.source_slot])),
            initial_interaction_density_offdiagonal_frobenius=float(np.linalg.norm(rb-np.diag(np.diag(rb)))),
            initial_state_hermiticity_frobenius=float(np.linalg.norm(self.rho_b0-self.rho_b0.conj().T)),
            initial_state_trace_real=float(np.trace(self.rho_b0).real),
            initial_state_trace_imag=float(np.trace(self.rho_b0).imag),
            initial_mean_force_extra_reference_subtraction=True,
            density_diagonalizations=0,state_repairs=0,hamiltonian_repairs=0)

    @staticmethod
    def _require(value,shape,name):
        if np.shape(value)!=shape or not np.isfinite(value).all():
            raise ValueError(name+' must be finite with shape '+str(shape))

    def _split(self,value):
        value=np.asarray(value,dtype=complex)
        if value.ndim not in (1,2) or value.shape[0]!=self.size:
            raise ValueError('state/dual must have shape (hybrid state size[, probes])')
        matrix=value[:,None] if value.ndim==1 else value
        p=matrix.shape[1]
        x=matrix[:self.nx]
        other=matrix[self.nx:self.nx+self._other_size].T.reshape(p,len(self.other_slots),self.c.m,self.c.m)
        bath=matrix[self.nx+self._other_size:].T.reshape(p,self.b,self.b)
        return x,other,bath,value.ndim

    @staticmethod
    def _join(x,other,bath,ndim):
        p=x.shape[1]
        result=np.concatenate((x,other.reshape(p,-1).T,bath.reshape(p,-1).T))
        return result[:,0] if ndim==1 else result

    def _moving(self,time):
        phase=np.exp(-1j*self.energies*float(time))
        inverse=np.exp(1j*self.energies*float(time))
        left=self._left0*phase
        right=inverse[:,None]*self._right0
        force=inverse[None,:,None]*self._delta_eigen*phase[None,None,:]
        return left,right,force

    def _compress(self,other,bath,left,right):
        active=np.empty((len(bath),self.c.ns,self.c.m,self.c.m),complex)
        active[:,self.other_slots]=other
        active[:,self.source_slot]=left@bath@right
        return active

    def _active_h(self,x,delta):
        c=self.c; data=self.data
        density=np.einsum('srab,psab->pr',data.density,delta,optimize=True)/self.root_n
        potential=np.einsum('rt,pt->pr',data.kernel,density,optimize=True)
        return (c.h0[None]+np.einsum('siab,ip->psab',data.ionic,x[:c.d],optimize=True)/self.root_n
            +np.einsum('sabr,pr->psab',data.potential,potential,optimize=True)/self.root_n
            +np.einsum('sfab,fp->psab',c.field_modes,x[2*c.d:2*c.d+c.f],optimize=True))

    def _extra_h(self,x,force):
        c=self.c
        amplitudes=np.concatenate((self.directions.conj().T@x[:c.d],self.directions.T@x[:c.d]))
        return -.5*np.einsum('kp,kab->pab',amplitudes,force,optimize=True)/self.root_n

    def _gradient(self,active,bath,force,*,subtract_initial):
        c=self.c
        gradient=np.zeros((self.nx,len(bath)),complex)
        gradient[:c.d]=np.einsum('siab,psab->ip',self.data.ionic_reverse,active,optimize=True)/self.root_n
        gradient[2*c.d:2*c.d+c.f]=np.einsum('sfab,psab->fp',c.field_reverse_modes,active,optimize=True)
        traces=np.einsum('kab,pba->kp',force,bath,optimize=True)
        if subtract_initial:
            traces=traces-self._initial_force_traces[:,None]
        gradient[:c.d]-=.5*(self.directions.conj()@traces[:2]+self.directions@traces[2:])/self.root_n
        return gradient

    def _mean_terms(self,time,mean):
        x,other,bath,ndim=self._split(mean)
        if ndim!=1:
            raise ValueError('mean must be one flat state vector')
        left,right,force=self._moving(time)
        active=self._compress(other,bath,left,right)
        h=self._active_h(x,active-self.c.rho0[None])
        hb=right@(h[:,self.source_slot]-self._initial_active_h)@left+self._extra_h(x,force)
        return x,other,bath,left,right,force,active,h,hb

    def mean_rhs(self,time,mean):
        x,other,bath,left,right,force,active,h,hb=self._mean_terms(time,mean)
        gradient=self._gradient(active-self.c.rho0[None],bath,force,subtract_initial=True)
        dx=self.c.mechanical.generators[0]@x+self.c.mechanical.j[0]@(self.c.mechanical.linear[:,None]+gradient)
        return self._join(dx,_commutator(h[:,self.other_slots],other),_commutator(hb,bath),1)

    def tangent(self,time,mean,value):
        """Exact complex-linear mean derivative, with batched canonical columns."""
        _,rho,rhob,left,right,force,_,h,hb=self._mean_terms(time,mean)
        x,sigma,sigmab,ndim=self._split(value)
        active=self._compress(sigma,sigmab,left,right)
        # Remove only the constant H0 of the active Hamiltonian builder.
        dh=self._active_h(x,active)-self.c.h0[None]
        dhb=right@dh[:,self.source_slot]@left+self._extra_h(x,force)
        dx=self.c.mechanical.generators[0]@x+self.c.mechanical.j[0]@self._gradient(active,sigmab,force,subtract_initial=False)
        ds=_commutator(h[:,self.other_slots],sigma)+_commutator(dh[:,self.other_slots],rho)
        db=_commutator(hb,sigmab)+_commutator(dhb,rhob)
        return self._join(dx,ds,db,ndim)

    def transpose(self,time,mean,dual):
        """True NONCONJUGATING transpose of tangent; no conjugate/time reversal."""
        c=self.c; data=self.data
        _,rho,rhob,left,right,force,_,h,hb=self._mean_terms(time,mean)
        zx,zs,zb,ndim=self._split(dual)
        outs,dh_other=_commutator_duals(h[:,self.other_slots],rho,zs)
        outb,dh_bath=_commutator_duals(hb,rhob,zb)
        p=zx.shape[1]
        dh=np.empty((p,c.ns,c.m,c.m),complex)
        dh[:,self.other_slots]=dh_other
        dh[:,self.source_slot]=right.T@dh_bath@left.T
        outx=c.mechanical.generators[0].T@zx
        outx[:c.d]+=np.einsum('siab,psab->ip',data.ionic,dh,optimize=True)/self.root_n
        outx[2*c.d:2*c.d+c.f]+=np.einsum('sfab,psab->fp',c.field_modes,dh,optimize=True)
        force_dual=-.5*np.einsum('kab,pab->kp',force,dh_bath,optimize=True)/self.root_n
        outx[:c.d]+=self.directions.conj()@force_dual[:2]+self.directions@force_dual[2:]
        potential_dual=np.einsum('sabr,psab->pr',data.potential,dh,optimize=True)/self.root_n
        density_dual=potential_dual@data.kernel
        active_dual=np.einsum('srab,pr->psab',data.density,density_dual,optimize=True)/self.root_n
        gradient_dual=c.mechanical.j[0].T@zx
        active_dual+=np.einsum('siab,ip->psab',data.ionic_reverse,gradient_dual[:c.d],optimize=True)/self.root_n
        active_dual+=np.einsum('sfab,fp->psab',c.field_reverse_modes,gradient_dual[2*c.d:2*c.d+c.f],optimize=True)
        outs+=active_dual[:,self.other_slots]
        outb+=left.T@active_dual[:,self.source_slot]@right.T
        trace_dual=np.concatenate((self.directions.conj().T@gradient_dual[:c.d],self.directions.T@gradient_dual[:c.d]))
        outb-=.5*np.einsum('kp,kba->pab',trace_dual,force,optimize=True)/self.root_n
        return self._join(outx,outs,outb,ndim)

    def probe(self,coordinate):
        if not isinstance(coordinate,(int,np.integer)) or not 0<=coordinate<self.c.d:
            raise ValueError('original ionic Cartesian displacement coordinate required')
        out=np.zeros(self.size,complex)
        out[int(coordinate)]=1
        return out

    def preparation_observables(self,dual):
        """Physical bath observable and its actual physical adjoint at t=0.

        V is not assumed unitary: O=V Z_I.T V^-1; Odagger is formed in
        physical coordinates, NOT by pretending Z_I.conj() is its dual.
        """
        _,_,bath,_=self._split(dual)
        observable=self.v@_transpose(bath)@self.vi
        return observable,_transpose(observable).conj()

    def initial_ordered(self,left,right):
        """Original product preparation ordered covariance, no cell/spin factor."""
        lx,ls,_,_=self._split(left)
        rx,rs,_,_=self._split(right)
        result=lx.T@self.c.mechanical.w0[0]@rx
        rho=self.c.rho0[self.other_slots]
        hole=np.eye(self.c.m)[None]-rho
        result+=np.einsum('psab,rscd,scb,sad->pr',ls,rs,rho,hole,optimize=True)
        ol,_=self.preparation_observables(left)
        ort,_=self.preparation_observables(right)
        # The first CAR term acts through the FULL bath identity, including
        # complement of the original frame; it is not an independent noise bank.
        frame=self.frame; density=self.density_coefficients
        small_l=frame.conj().T@ol@frame
        small_r=frame.conj().T@ort@frame
        for p in range(len(ol)):
            for r in range(len(ort)):
                car=frame.conj().T@ol[p]@ort[r]@frame
                result[p,r]+=np.trace(density@car-density@small_l[p]@density@small_r[r])
        return result

    def state_diagnostics(self,time,mean):
        """Raw physical state defects at a requested checkpoint; no repairs."""
        _,other,bath,_,_,_,_,_,_=self._mean_terms(time,mean)
        phase=np.exp(-1j*self.energies*float(time))
        inverse=np.exp(1j*self.energies*float(time))
        physical=(self.v*phase)@bath[0]@(inverse[:,None]*self.vi)
        return dict(source0_hermiticity_frobenius=float(np.linalg.norm(physical-physical.conj().T)),
            source0_trace_real=float(np.trace(physical).real),source0_trace_imag=float(np.trace(physical).imag),
            source0_trace_change_real=float((np.trace(physical)-np.trace(self.rho_b0)).real),
            source0_trace_change_imag=float((np.trace(physical)-np.trace(self.rho_b0)).imag),
            other_active_hermiticity_frobenius=float(np.linalg.norm(other-_transpose(other).conj())))
