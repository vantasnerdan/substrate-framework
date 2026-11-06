"""Source-complete cubic scalar/mean-force map and true bilinear transpose.

Cartesian lattice gauge, Ry/bohr units; canonical1/sqrt64 is retained. No
all-q111^3 object is constructed. The ordinary-NC electronic term is onsite;
the finite Ewald term is represented by real pairs and compact reciprocal
pair tensors from the already-qualified A43 provider. All111 coordinates are
available, including arbitrary complex non-Hermitian probe covariances.
"""
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent;A43=HERE.parent/'0043'
sys.path.insert(0,str(A43))
from cubic_ionic_donor import IonicThird
from cubic_material_source import MATCHING,BOUNDARY
from substrate_framework.streamed_material_state import authenticate,file_binding

PIECES=('total','ionic_real','ionic_reciprocal','electronic_local','electronic_nonlocal')


class CompactCubic:
    def __init__(self,ionic,local,nonlocal_term):
        self.ionic=ionic;self.nat=ionic.nat;self.d=3*self.nat
        self.mesh=np.asarray(ionic.mesh,int);self.n=int(np.prod(self.mesh))
        self.local=np.asarray(local);self.nonlocal_term=np.asarray(nonlocal_term)
        if self.local.shape!=(self.nat,3,3,3) or self.nonlocal_term.shape!=self.local.shape:
            raise ValueError('complete onsite thirds for every original atom required')
        pairs=ionic.real_pairs
        self.i=np.asarray([pair[0] for pair in pairs],int)
        self.j=np.asarray([pair[1] for pair in pairs],int)
        self.lattice=np.asarray([pair[2] for pair in pairs],int).reshape(-1,3)
        self.real=np.asarray([pair[3] for pair in pairs]).reshape(-1,3,3,3)/np.sqrt(self.n)
        self.selected=tuple(range(self.nat))
        self.charge_pairs=np.outer(ionic.charges,ionic.charges)
        self.reciprocal_prefactor=4*np.pi/ionic.volume/np.sqrt(self.n)
        self.indices=np.indices(tuple(self.mesh)).reshape(3,-1).T
        self.minus=np.ravel_multi_index(tuple((-self.indices%self.mesh).T),tuple(self.mesh))

    @classmethod
    def from_source(cls,path):
        binding=file_binding(path);record=json.loads(Path(path).read_bytes())
        if (record['status']!='complete' or record['schema']!='p254-a43-full-state-third-v2'
                or record['electronic_atom_indices']!=list(range(37))
                or record['matching']!=MATCHING or record['original_boundary']!=BOUNDARY
                or record['atom_names']!=['Ni']*31+['H']*6
                or record['coefficient_manifest']['sha256']!='f851d2df2b2f1136d22b606f059cf6cd9bf6df7c4681340153e12c8569dcd479'
                or len(record['full_sector_lineage'])!=128
                or len(record['completed_full37_sector_archives'])!=128):
            raise ValueError('actual source-complete128sector/all37electronic-atom source required')
        for source in record['source_bindings']:
            authenticate(source)
        with np.load(authenticate(record['arrays']),allow_pickle=False) as saved:
            if not np.array_equal(saved['electronic_atom_indices'],np.arange(37)):
                raise ValueError('saved full37 atom coverage differs')
            ionic=IonicThird(saved['cell'],saved['positions'],saved['charges'],(4,4,4))
            result=cls(ionic,saved['local'],saved['nonlocal_term'])
        result.source_binding=binding;result.source_record=record
        return result

    def qindex(self,q):
        return int(np.ravel_multi_index(tuple(np.asarray(q,int)%self.mesh),tuple(self.mesh)))

    def _qs(self,qs,piece):
        q=np.asarray(qs,int)
        if q.shape!=(3,3) or np.any(q.sum(axis=0)%self.mesh) or piece not in PIECES:
            raise ValueError('three exact conserving original mesh labels and explicit donor piece required')
        return q

    def _vector(self,value):
        value=np.asarray(value)
        if value.shape!=(self.d,) or not np.isfinite(value).all():
            raise ValueError('finite complete atom-major Cartesian coordinate vector required')
        return value.reshape(self.nat,3)

    def _covariance(self,value):
        value=np.asarray(value)
        if value.shape!=(self.d,self.d) or not np.isfinite(value).all():
            raise ValueError('finite complete two-leg covariance; no Hermitian repair allowed')
        return value.reshape(self.nat,3,self.nat,3).transpose(0,2,1,3)

    def _reciprocal(self,q):
        for mask in range(8):
            subset=[leg for leg in range(3) if mask&(1<<leg)]
            shift=np.sum(q[subset],axis=0) if subset else np.zeros(3,int)
            pair,_=self.ionic._reciprocal_pairs(shift,self.selected)
            yield mask,(-1)**(3-len(subset))*self.reciprocal_prefactor*self.charge_pairs,pair

    def contract(self,qs,y=None,z=None,covariance=None,piece='total'):
        """Phi[first,second,third]:(y second*z third+covariance), all first rows.

        Inputs remain independent complex objects. No conjugation, implicit
        symmetry, mass factors or covariance ordering replacement occurs.
        """
        q=self._qs(qs,piece)
        if (y is None)!=(z is None) or (y is None and covariance is None):
            raise ValueError('provide both coordinate legs and/or an explicit covariance')
        yv=None if y is None else self._vector(y)
        zv=None if z is None else self._vector(z)
        w=None if covariance is None else self._covariance(covariance)
        out=np.zeros((self.nat,3),complex)
        if piece in ('total','electronic_local','electronic_nonlocal'):
            onsite=(self.local+self.nonlocal_term if piece=='total' else
                    self.local if piece=='electronic_local' else self.nonlocal_term)/np.sqrt(self.n)
            moment=np.zeros((self.nat,3,3),complex)
            if w is not None:
                moment+=w[np.arange(self.nat),np.arange(self.nat)]
            if yv is not None:
                moment+=yv[:,:,None]*zv[:,None,:]
            out+=np.einsum('nabc,nbc->na',onsite,moment,optimize=True)
        if piece in ('total','ionic_real'):
            phases=np.exp(2j*np.pi*(self.lattice@(q/self.mesh).T))
            moment=np.zeros((len(self.i),3,3),complex)
            if w is not None:
                moment+=(w[self.i,self.i]-phases[:,2,None,None]*w[self.i,self.j]
                    -phases[:,1,None,None]*w[self.j,self.i]
                    +(phases[:,1]*phases[:,2])[:,None,None]*w[self.j,self.j])
            if yv is not None:
                left=yv[self.i]-phases[:,1,None]*yv[self.j]
                right=zv[self.i]-phases[:,2,None]*zv[self.j]
                moment+=left[:,:,None]*right[:,None,:]
            value=np.einsum('nabc,nbc->na',self.real,moment,optimize=True)
            np.add.at(out,self.i,value)
            np.add.at(out,self.j,-phases[:,0,None]*value)
        if piece in ('total','ionic_reciprocal'):
            diagonal=None if w is None else w[np.arange(self.nat),np.arange(self.nat)]
            for mask,coefficient,pair in self._reciprocal(q):
                moment=np.zeros((self.nat,self.nat,3,3),complex)
                if w is not None:
                    if mask&2 and mask&4:moment+=diagonal[:,None]
                    elif not mask&2 and not mask&4:moment+=diagonal[None,:]
                    elif mask&2:moment+=w
                    else:moment+=w.transpose(1,0,2,3)
                if yv is not None:
                    left=yv[:,None,:] if mask&2 else yv[None,:,:]
                    right=zv[:,None,:] if mask&4 else zv[None,:,:]
                    moment+=left[:,:,:,None]*right[:,:,None,:]
                value=np.einsum('ijabc,ijbc->ija',pair,moment,optimize=True)*coefficient[:,:,None]
                out+=value.sum(axis=1 if mask&1 else 0)
        return out.reshape(self.d)

    def covariance_transpose(self,qs,cotangent,piece='total'):
        """Bilinear transpose of W -> Phi:W, returning the complete W dual."""
        q=self._qs(qs,piece);cot=self._vector(cotangent)
        out=np.zeros((self.nat,self.nat,3,3),complex)
        if piece in ('total','electronic_local','electronic_nonlocal'):
            onsite=(self.local+self.nonlocal_term if piece=='total' else
                    self.local if piece=='electronic_local' else self.nonlocal_term)/np.sqrt(self.n)
            out[np.arange(self.nat),np.arange(self.nat)]+=np.einsum('nabc,na->nbc',onsite,cot,optimize=True)
        if piece in ('total','ionic_real'):
            phases=np.exp(2j*np.pi*(self.lattice@(q/self.mesh).T))
            effective=cot[self.i]-phases[:,0,None]*cot[self.j]
            value=np.einsum('nabc,na->nbc',self.real,effective,optimize=True)
            np.add.at(out,(self.i,self.i),value)
            np.add.at(out,(self.i,self.j),-phases[:,2,None,None]*value)
            np.add.at(out,(self.j,self.i),-phases[:,1,None,None]*value)
            np.add.at(out,(self.j,self.j),(phases[:,1]*phases[:,2])[:,None,None]*value)
        if piece in ('total','ionic_reciprocal'):
            for mask,coefficient,pair in self._reciprocal(q):
                effective=cot[:,None,:] if mask&1 else cot[None,:,:]
                effective=np.broadcast_to(effective,(self.nat,self.nat,3))
                value=np.einsum('ijabc,ija->ijbc',pair,effective,optimize=True)*coefficient[:,:,None,None]
                if mask&2 and mask&4:
                    out[np.arange(self.nat),np.arange(self.nat)]+=value.sum(axis=1)
                elif not mask&2 and not mask&4:
                    out[np.arange(self.nat),np.arange(self.nat)]+=value.sum(axis=0)
                elif mask&2:out+=value
                else:out+=value.transpose(1,0,2,3)
        return out.transpose(0,2,1,3).reshape(self.d,self.d)

    def ordered_scalar(self,qs,x,y,z):
        """One ordered cubic term Phi[x,y,z]/6; full field scalar sums all q."""
        return np.dot(self._vector(x).reshape(self.d),self.contract(qs,y,z))/6

    def scalar_jvp(self,qs,x,y,z,dx,dy,dz):
        return (self.ordered_scalar(qs,dx,y,z)+self.ordered_scalar(qs,x,dy,z)
                +self.ordered_scalar(qs,x,y,dz))

    def scalar_transpose(self,qs,x,y,z,cotangent):
        q=self._qs(qs,'total')
        return (cotangent*self.contract(q,y,z)/6,
                cotangent*self.contract(q[[1,0,2]],x,z)/6,
                cotangent*self.contract(q[[2,0,1]],x,y)/6)

    def mean_force(self,qs,y,z,covariance):
        return -.5*self.contract(qs,y,z,covariance)

    def mean_force_jvp(self,qs,y,z,dy,dz,dcovariance):
        return -.5*(self.contract(qs,dy,z)+self.contract(qs,y,dz)+self.contract(qs,covariance=dcovariance))

    def mean_force_transpose(self,qs,y,z,cotangent):
        q=self._qs(qs,'total');cot=self._vector(cotangent).reshape(self.d)
        return (-.5*self.contract(q[[1,0,2]],cot,z),
                -.5*self.contract(q[[2,0,1]],cot,y),
                -.5*self.covariance_transpose(q,cot))

    def routed_mean_force(self,external_q,means,covariance_provider):
        """Full64 internal-q sum for any externalQ and explicit covariance(k,l).

        The provider must supply its actual111x111 ordered covariance for every
        routed pair; a missing/unknown covariance is never fabricated as zero.
        """
        means=np.asarray(means)
        if means.shape!=(self.n,self.d) or not np.isfinite(means).all():
            raise ValueError('complete64q/all111-coordinate means required')
        result=np.zeros(self.d,complex)
        for k in range(self.n):
            l=self.qindex(self.indices[external_q]-self.indices[k])
            qs=self.indices[[int(self.minus[external_q]),k,l]]
            result+=self.mean_force(qs,means[k],means[l],covariance_provider(k,l))
        return result

    def routed_mean_jvp(self,external_q,means,dmeans,dcovariance_provider):
        means=np.asarray(means);dmeans=np.asarray(dmeans)
        if means.shape!=(self.n,self.d) or dmeans.shape!=means.shape:
            raise ValueError('complete original mean and tangent axes required')
        result=np.zeros(self.d,complex)
        for k in range(self.n):
            l=self.qindex(self.indices[external_q]-self.indices[k])
            qs=self.indices[[int(self.minus[external_q]),k,l]]
            result+=self.mean_force_jvp(qs,means[k],means[l],dmeans[k],dmeans[l],dcovariance_provider(k,l))
        return result

    def routed_mean_transpose(self,external_q,means,cotangent):
        """Return mean dual and one covariance dual for each of64 routed pairs."""
        means=np.asarray(means)
        if means.shape!=(self.n,self.d):raise ValueError('complete original mean axes required')
        mean_dual=np.zeros_like(means,dtype=complex)
        covariance_dual=np.empty((self.n,self.d,self.d),complex)
        for k in range(self.n):
            l=self.qindex(self.indices[external_q]-self.indices[k])
            qs=self.indices[[int(self.minus[external_q]),k,l]]
            left,right,matrix=self.mean_force_transpose(qs,means[k],means[l],cotangent)
            mean_dual[k]+=left;mean_dual[l]+=right;covariance_dual[k]=matrix
        return mean_dual,covariance_dual

    def field_scalar(self,means):
        """Complete coherent cubic scalar; skip only exactly zero input terms."""
        means=np.asarray(means)
        if means.shape!=(self.n,self.d) or not np.isfinite(means).all():raise ValueError('complete finite64q/all111-coordinate means required')
        active=np.flatnonzero(np.any(means!=0,axis=1));result=0j
        for q1 in active:
            for q2 in active:
                q3=self.qindex(-self.indices[q1]-self.indices[q2])
                if np.any(means[q3]!=0):
                    result+=self.ordered_scalar(self.indices[[q1,q2,q3]],means[q1],means[q2],means[q3])
        return result

    def field_force(self,means):
        means=np.asarray(means)
        if means.shape!=(self.n,self.d) or not np.isfinite(means).all():raise ValueError('complete finite64q/all111-coordinate means required')
        active=np.flatnonzero(np.any(means!=0,axis=1));result=np.zeros_like(means,dtype=complex)
        for k in active:
            for l in active:
                external=self.qindex(self.indices[k]+self.indices[l])
                result[external]-=.5*self.contract(self.indices[[int(self.minus[external]),k,l]],means[k],means[l])
        return result

    def field_force_jvp(self,means,dmeans):
        means=np.asarray(means);dmeans=np.asarray(dmeans)
        if means.shape!=(self.n,self.d) or dmeans.shape!=means.shape:
            raise ValueError('complete64q/all111-coordinate mean/tangent required')
        active=np.flatnonzero(np.any(means!=0,axis=1))
        tangent=np.flatnonzero(np.any(dmeans!=0,axis=1))
        result=np.zeros_like(means,dtype=complex)
        for k in tangent:
            for l in active:
                external=self.qindex(self.indices[k]+self.indices[l])
                result[external]-=.5*self.contract(self.indices[[int(self.minus[external]),k,l]],dmeans[k],means[l])
        for k in active:
            for l in tangent:
                external=self.qindex(self.indices[k]+self.indices[l])
                result[external]-=.5*self.contract(self.indices[[int(self.minus[external]),k,l]],means[k],dmeans[l])
        return result

    def field_force_transpose(self,means,cotangent):
        means=np.asarray(means);cotangent=np.asarray(cotangent)
        if means.shape!=(self.n,self.d) or cotangent.shape!=means.shape:
            raise ValueError('complete64q/all111-coordinate mean/cotangent required')
        active=np.flatnonzero(np.any(means!=0,axis=1))
        result=np.zeros_like(means,dtype=complex)
        for k in range(self.n):
            for l in active:
                external=self.qindex(self.indices[k]+self.indices[l])
                if np.any(cotangent[external]!=0):
                    result[k]-=.5*self.contract(self.indices[[k,int(self.minus[external]),l]],cotangent[external],means[l])
        for k in active:
            for l in range(self.n):
                external=self.qindex(self.indices[k]+self.indices[l])
                if np.any(cotangent[external]!=0):
                    result[l]-=.5*self.contract(self.indices[[l,int(self.minus[external]),k]],cotangent[external],means[k])
        return result

    def field_mean_force(self,means,covariances):
        """All external coherent forces plus a declared translation-invariant
        connected Gaussian covariance W[k,-k], with the full64 internal sum.

        This input class is explicit, not an inferred FDT/stationarity closure.
        General non-diagonal momentum covariances use routed_mean_force.
        """
        covariances=np.asarray(covariances)
        if covariances.shape!=(self.n,self.d,self.d):
            raise ValueError('all64 raw ordered111x111 Gaussian covariances required')
        result=self.field_force(means)
        for k in range(self.n):
            qs=self.indices[[0,k,int(self.minus[k])]]
            result[0]-=.5*self.contract(qs,covariance=covariances[k])
        return result

    def field_mean_jvp(self,means,dmeans,dcovariances):
        dcovariances=np.asarray(dcovariances)
        if dcovariances.shape!=(self.n,self.d,self.d):
            raise ValueError('all64 raw covariance tangents required')
        result=self.field_force_jvp(means,dmeans)
        for k in range(self.n):
            qs=self.indices[[0,k,int(self.minus[k])]]
            result[0]-=.5*self.contract(qs,covariance=dcovariances[k])
        return result

    def field_mean_transpose(self,means,cotangent):
        """True complex-bilinear transpose; no cotangent is conjugated."""
        cotangent=np.asarray(cotangent)
        if cotangent.shape!=(self.n,self.d):
            raise ValueError('all64 external111-coordinate cotangents required')
        mean_dual=self.field_force_transpose(means,cotangent)
        covariance_dual=np.empty((self.n,self.d,self.d),complex)
        for k in range(self.n):
            qs=self.indices[[0,k,int(self.minus[k])]]
            covariance_dual[k]=-.5*self.covariance_transpose(qs,cotangent[0])
        return mean_dual,covariance_dual

    def field_scalar_jvp(self,means,dmeans):
        dmeans=np.asarray(dmeans)
        if dmeans.shape!=(self.n,self.d):raise ValueError('complete64q coordinate tangent required')
        return -np.sum(self.field_force(means)*dmeans[self.minus])

    def field_scalar_transpose(self,means,cotangent):
        return -cotangent*self.field_force(means)[self.minus]

    def field_mean_scalar(self,means,covariances):
        """Actual Gaussian expectation E3(mu,W), using the same complete Phi.

        E3=<Phi u^3/6>=E3(mu)+1/2 sum_k Phi[0,k,-k]:mu0 Wk.
        Connected third Gaussian moment is zero; no stationarity is assumed.
        Its coordinate derivative is exactly the mean force above.
        """
        means=np.asarray(means);covariances=np.asarray(covariances)
        if covariances.shape!=(self.n,self.d,self.d):raise ValueError('all64 ordered covariances required')
        scalar=self.field_scalar(means)
        for k in range(self.n):
            qs=self.indices[[0,k,int(self.minus[k])]]
            scalar+=.5*np.dot(means[0],self.contract(qs,covariance=covariances[k]))
        return scalar

    def field_mean_scalar_jvp(self,means,covariances,dmeans,dcovariances):
        means=np.asarray(means);dmeans=np.asarray(dmeans)
        covariances=np.asarray(covariances);dcovariances=np.asarray(dcovariances)
        if dmeans.shape!=(self.n,self.d) or covariances.shape!=(self.n,self.d,self.d) or dcovariances.shape!=covariances.shape:
            raise ValueError('all original mean/covariance tangent axes required')
        value=self.field_scalar_jvp(means,dmeans)
        for k in range(self.n):
            qs=self.indices[[0,k,int(self.minus[k])]]
            value+=.5*(np.dot(dmeans[0],self.contract(qs,covariance=covariances[k]))
                +np.dot(means[0],self.contract(qs,covariance=dcovariances[k])))
        return value

    def field_mean_scalar_transpose(self,means,covariances,cotangent):
        means=np.asarray(means)
        mean_dual=-cotangent*self.field_mean_force(means,covariances)[self.minus]
        covariance_dual=np.empty((self.n,self.d,self.d),complex)
        for k in range(self.n):
            qs=self.indices[[0,k,int(self.minus[k])]]
            covariance_dual[k]=.5*cotangent*self.covariance_transpose(qs,means[0])
        return mean_dual,covariance_dual
