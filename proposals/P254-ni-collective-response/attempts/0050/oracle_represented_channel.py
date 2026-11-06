"""New electron-only boundary against original coupled tangent and Fock noise."""
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from represented_electron_channel import RepresentedElectronChannel,polarization_difference
from reduced_electronic_correlations import ActiveBandCorrelations
from oracle_fixture import fock_ensemble,bilinear


def check(c,moving_baseline):
    baseline=lambda t:moving_baseline(0.)
    channel=RepresentedElectronChannel(c,baseline,rtol=2e-10,atol=2e-12)
    original=ActiveBandCorrelations(c,baseline)
    nx,ne=4,4
    full=original.tangent(0,0.,np.eye(nx+ne,dtype=complex))
    E=channel.tangent(0.,np.eye(ne).reshape(1,2,2,ne)).reshape(ne,ne)
    B=channel.forcing(0.,np.eye(2)).reshape(ne,2)
    C=channel.reverse.transpose(1,0,2,3).reshape(2,ne)
    ids=np.array([0,2]);P=np.eye(4)[:,ids]
    E_error=float(np.linalg.norm(E-full[nx:,nx:]))
    B_error=float(np.linalg.norm(B-full[nx:,ids]))
    C_error=float(np.linalg.norm(c.mechanical.j[0]@P@C-full[:nx,nx:]))
    assert max(E_error,B_error,C_error)<1e-13
    t=.73
    duals,receipt=channel.pullback(t,[0.,t])
    response=channel.response_from_dual(duals[0.],0.)
    expected=C@expm(t*E)@B
    response_error=float(np.linalg.norm(response-expected))
    assert response_error<1e-9
    # Schur/Volterra realization driven by original mechanical trajectories.
    initial=np.array([.13,-.07,.21,.04,.03+.02j,.01j,-.02j,-.03])
    def eliminated(t,y):
        x,s=y[:4],y[4:].reshape(1,2,2,1)
        dx=c.mechanical.generators[0]@x+c.mechanical.j[0]@P@channel.gradient(s)[:,0]
        ds=channel.tangent(t,s,x[ids,None])
        return np.r_[dx,ds.ravel()]
    solved=solve_ivp(eliminated,(0,t),initial,rtol=2e-11,atol=2e-13,method='DOP853')
    explicit=expm(t*full)@initial
    elimination_error=float(np.linalg.norm(solved.y[:,-1]-explicit))
    assert elimination_error<1e-9
    density,ann=fock_ensemble(c.rho0[0])
    left=duals[0.];right=channel.reverse.transpose(0,2,3,1)
    ordered=channel.ordered_noise(left,right)
    fock=np.zeros((2,2),complex)
    for i in range(2):
        a=bilinear(left[0,:,:,i],ann)
        for j in range(2):
            b=bilinear(right[0,:,:,j],ann)
            fock[i,j]=np.trace(density@a@b)-np.trace(density@a)*np.trace(density@b)
    fock_error=float(np.linalg.norm(ordered-fock));assert fock_error<1e-10
    # Original full Gaussian tangent: mechanical noise due to precisely this
    # same electronic initial sector, no additional marginal force bath.
    W=np.empty((ne,ne),complex)
    probes=np.eye(ne).reshape(1,2,2,ne)
    W[:]=channel.ordered_noise(probes,probes)
    F=expm(t*full)
    exact_noise=F[:nx,nx:]@W@F[:nx,nx:].T
    # Independent Duhamel convolution C Ue(t,s) with bare mechanical evolution.
    def transfer_rhs(t,z):
        X=z[:nx*ne].reshape(nx,ne);U=z[nx*ne:].reshape(ne,ne)
        return np.r_[(c.mechanical.generators[0]@X+c.mechanical.j[0]@P@C@U).ravel(),
                     (E@U+B@P.T@X).ravel()]
    transfer=solve_ivp(transfer_rhs,(0,t),np.r_[np.zeros(nx*ne),np.eye(ne).ravel()].astype(complex),
                      rtol=2e-11,atol=2e-13,method='DOP853')
    T=transfer.y[:nx*ne,-1].reshape(nx,ne)
    noise_error=float(np.linalg.norm(T@W@T.T-exact_noise));assert noise_error<1e-9
    kernel=channel.block.kernel.copy();channel.block.kernel[:]=0
    E_missing=channel.tangent(0.,probes).reshape(ne,ne)
    channel.block.kernel[:]=kernel
    missing_k=float(np.linalg.norm(C@expm(t*E_missing)@B-expected))
    stacked_response=float(np.linalg.norm(2*expected-expected))
    stacked_noise=float(np.linalg.norm(2*exact_noise-exact_noise))
    wrong_hole=np.einsum('sabp,scdr,scb,sad->pr',left,right,c.rho0,c.rho0,optimize=True)
    pauli=float(np.linalg.norm(wrong_hole-ordered))
    assert min(missing_k,stacked_response,stacked_noise,pauli)>1e-6
    packed={name:np.ones((2,2)) for name in ('pi_retarded','pi_greater','pi_lesser')}
    assert all(np.max(abs(v))==0 for v in polarization_difference(packed,packed).values())
    # Time-dependent original mean actually traverses this public transpose API.
    moving=RepresentedElectronChannel(c,moving_baseline,rtol=2e-10,atol=2e-12)
    dm,_=moving.pullback(t,[0.,t])
    moving_difference=float(np.linalg.norm(dm[0.]-duals[0.]))
    assert moving_difference>1e-6
    batched,_=moving.pullback(t,[0.,t],[0])
    batch_error=float(np.linalg.norm(batched[0.]-dm[0.][:,:,:,0:1]))
    assert batch_error<1e-10
    return dict(block_errors=dict(E=E_error,B=B_error,C=C_error),response_error=response_error,
        elimination_error=elimination_error,fock_error=fock_error,explicit_coupled_electronic_noise_error=noise_error,batch_error=batch_error,
        sensitivities=dict(missing_K=missing_k,stacked_response=stacked_response,stacked_noise=stacked_noise,
                           incorrect_Pauli_hole=pauli,moving_mean=moving_difference),pullback=receipt)
