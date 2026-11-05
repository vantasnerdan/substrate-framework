"""One new composition oracle; admitted A29/A47 proofs are not replayed."""
import numpy as np
from reduced_electronic_dynamics import BaselineState

def boundary_oracle(model,rho,w,reference):
    c=model.c;d=c.d;nx=model.nx
    rng=np.random.default_rng(5701)
    x=(rng.normal(size=nx)+.2j*rng.normal(size=nx))*.003
    dx=(rng.normal(size=nx)+.3j*rng.normal(size=nx))*.01
    dr=(rng.normal(size=rho.shape)+.7j*rng.normal(size=rho.shape))*.001
    dw=(rng.normal(size=w.shape)+.6j*rng.normal(size=w.shape))*.001
    qs=model.source.indices[[0,0,0]]
    cf=-.5*model.source.contract(qs,covariance=w)
    dcf=-.5*model.source.contract(qs,covariance=dw)
    eps=2e-5
    hfd=(model.hamiltonian(rho+eps*dr,x+eps*dx)-model.hamiltonian(rho-eps*dr,x-eps*dx))/(2*eps)
    hj=model.hamiltonian_jvp(dr,dx)
    g=model.electronic_gradient(rho)
    fd=(model.mean_rhs(x+eps*dx,model.electronic_gradient(rho+eps*dr),cf+eps*dcf)-model.mean_rhs(x-eps*dx,model.electronic_gradient(rho-eps*dr),cf-eps*dcf))/(2*eps)
    tangent=model.mean_jvp(x,dx,dr,dcf)
    # Scalar differentiated with independent complex coordinate/density probes.
    sfd=(model.scalar(rho+eps*dr,x+eps*dx,cf+eps*dcf)-model.scalar(rho-eps*dr,x-eps*dx,cf-eps*dcf))/(2*eps)
    phase_gradient=c.mechanical.h[0]@x+c.mechanical.linear+g
    phase_gradient[:d]+=model.cubic_gradient(x,cf)
    sj=phase_gradient@dx+np.einsum('sab,sba->',model.hamiltonian(rho,x),dr)-x[:d]@dcf
    cubic_direct=model.source.mean_force(qs,x[:d],x[:d],w)
    tensor_force=-model.cubic_gradient(x,cf)
    cubic_j=model.source.mean_force_jvp(qs,x[:d],x[:d],dx[:d],dx[:d],dw)
    tensor_j=-np.einsum('ijk,k,j->i',model.phi,x[:d],dx[:d],optimize=True)+dcf
    covdual=model.cubic_covariance_derivative(x,0)
    e3fd=(model.cubic_scalar(x,cf+eps*dcf)-model.cubic_scalar(x,cf-eps*dcf))/(2*eps)
    expected_cov=np.sum(covdual*dw)
    # Actual newly used diagonal conversion against genuine saved old provider.
    dg,dl=reference.ordered(0,1,1,'positions')
    reconstructed=dl+(dg-dl)
    ix=np.r_[np.arange(d),np.arange(2*d,2*d+c.f)]
    expected_w=reference.S(0,reference.times[1])[ix]@c.mechanical.w0[0]@reference.S(0,reference.times[1])[ix].T
    rec=dict(scope='new homogeneous A29+A47 raw density/Gamma Gaussian force and same-scalar tangent composition',
        hamiltonian_jvp_error=float(np.max(np.abs(hfd-hj))),
        phase_mean_jvp_error=float(np.max(np.abs(fd-tangent))),
        common_scalar_jvp_error=float(abs(sfd-sj)),
        source_tensor_force_error=float(np.max(np.abs(cubic_direct-tensor_force))),
        source_tensor_jvp_error=float(np.max(np.abs(cubic_j-tensor_j))),
        covariance_scalar_derivative_error=float(abs(e3fd-expected_cov)),
        genuine_bare_ordered_covariance_error=float(np.max(np.abs(1j*reconstructed-expected_w))),
        wrong_cubic_sign_sensitivity=float(np.linalg.norm(cubic_direct+tensor_force)),
        missing_active_K_sensitivity=float(np.linalg.norm(hj-model.hamiltonian_jvp(np.zeros_like(dr),dx))),
        wrong_covariance_conjugation_sensitivity=float(abs(np.sum(covdual*(dw-dw.conj())))),
        frozen_contact_is_bubble=False,inherited_proofs_replayed=False)
    errors=[v for k,v in rec.items() if k.endswith('_error')]
    if not np.isfinite(errors).all() or max(errors)>2e-7 or min(rec[k] for k in ('wrong_cubic_sign_sensitivity','missing_active_K_sensitivity','wrong_covariance_conjugation_sensitivity'))<=1e-18:
        raise AssertionError(rec)
    return rec
