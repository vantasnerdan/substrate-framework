"""A29 common scalar derivatives plus the A47 Gaussian E3 at homogeneous mean.

Only Gamma has a coherent mean. All64 ordered covariances contribute to E3.
The fixed-rho phase Hessian is not an electronic bubble: that response belongs
in the common Pi, not a second contact insertion. Canonical coordinates retain
original masses, field momentum coupling, contacts, and passive complement.
"""
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parents[3]/'src')]+[str(HERE.parent/n) for n in ('0055','0052','0051','0047','0045','0042','0041','0030','0029')]
from refinement_common import binding,save,authenticate
from reduced_electronic_dynamics import ActiveBandDynamics,BaselineState
from cubic_feedback import CompactCubic
from cubic_contour import vertex

SOURCE=Path('/var/tmp/p254-a47-full37-source-001/source.json')

class CommonMeanFunctional:
    def __init__(self,c,source,phi=None):
        self.c=c;self.source=source;self.base=ActiveBandDynamics(c)
        self.data=c.q(0)
        self.phi=vertex(source,source.indices[[0,0,0]]) if phi is None else phi
        self.nx=c.mechanical_sizes[0]

    def electronic_gradient(self,rho):
        return self.base._gradient(self.data,rho-self.c.rho0,True)

    def hamiltonian(self,rho,x):
        return self.base.baseline_hamiltonian(BaselineState(rho,x),self.data)

    def hamiltonian_jvp(self,drho,dx):
        c=self.c;d=self.data
        density=self.base.density(d,drho)
        return (np.einsum('siab,i->sab',d.ionic,dx[:c.d])/c.sqrt_n
            +np.einsum('sabr,r->sab',d.potential,d.kernel@density)/c.sqrt_n
            +np.einsum('sfab,f->sab',c.field_modes,dx[2*c.d:2*c.d+c.f]))

    def cubic_scalar(self,x,covariance_force):
        u=x[:self.c.d]
        return np.einsum('ijk,i,j,k->',self.phi,u,u,u,optimize=True)/6-u@covariance_force

    def cubic_gradient(self,x,covariance_force):
        return .5*np.einsum('ijk,j,k->i',self.phi,x[:self.c.d],x[:self.c.d],optimize=True)-covariance_force

    def cubic_covariance_derivative(self,x,internal_q):
        s=self.source;k=int(internal_q)
        return .5*s.covariance_transpose(s.indices[[0,k,int(s.minus[k])]],x[:self.c.d])

    def curvature(self,x):
        h=self.c.mechanical.h[0].astype(complex,copy=True)
        h[:self.c.d,:self.c.d]+=np.einsum('ijk,k->ij',self.phi,x[:self.c.d],optimize=True)
        return h

    def mean_rhs(self,x,electronic_gradient,covariance_force):
        c=self.c;g=c.mechanical.linear+electronic_gradient
        # No sign guessing from momentum indices: use the original full J.
        extra=np.zeros(self.nx,complex);extra[:c.d]=self.cubic_gradient(x,covariance_force)
        return c.mechanical.generators[0]@x+c.mechanical.j[0]@(g+extra)

    def mean_jvp(self,x,dx,drho,dcovariance_force):
        c=self.c;g=self.base._gradient(self.data,drho,True)
        g[:c.d]-=dcovariance_force
        return c.mechanical.j[0]@(self.curvature(x)@dx+g)

    def scalar(self,rho,x,covariance_force):
        """Common scalar at fixed passive preparation (constant omitted).

        Gamma density bilinear uses explicit feature partner, not conjugation.
        Its rho derivative includes h0, ionic/field potential and active K.
        """
        c=self.c;d=self.data;density=self.base.density(d,rho-c.rho0)
        return (.5*x@c.mechanical.h[0]@x+c.mechanical.linear@x
            +np.einsum('sab,sba->',c.h0,rho)
            +x@self.electronic_gradient(rho)
            +.5*density[d.feature_partner]@(d.kernel@density)
            +self.cubic_scalar(x,covariance_force))
