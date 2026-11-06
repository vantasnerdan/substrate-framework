"""Small explicit differentiated scalar, canonical mechanics, and CAR ensemble."""
from types import SimpleNamespace
import numpy as np
from scipy.linalg import expm


def fixture():
    h = np.array([[.3,.1+.07j],[.1-.07j,-.2]])
    rho = np.array([[.73,.12+.06j],[.12-.06j,.24]])
    g = np.array([[.4,.2+.11j],[.2-.11j,-.3]])
    f = np.array([[-.1,.31-.08j],[.31+.08j,.27]])
    v = np.array([[.2,.4+.13j],[.4-.13j,-.17]])
    block = SimpleNamespace(ionic=g[None,None], ionic_reverse=g.T[None,None],
        density=v.T[None,None], potential=v[None,:,:,None], kernel=np.array([[1.7]]))
    # x=(u,p,a,Pi), including nonzero contact, mixed, minimal coupling blocks.
    mh = np.array([[.8,0,.13,0],[0,.6,-.19,0],[.13,-.19,.9,0],[0,0,0,.7]])
    j = np.array([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0.]])
    w0 = np.diag([.8,.6,.7,.9]).astype(complex)+.5j*j
    mech = SimpleNamespace(generators=[j@mh], j=[j], h=[mh], w0=[w0])
    c = SimpleNamespace(n=4,sqrt_n=2.,ns=1,m=2,d=1,f=1,target=np.array([[0]]),
        minus=np.array([0]),rho0=rho[None],h0=h[None],field_modes=f[None,None],
        field_reverse_modes=f.T[None,None],mechanical=mech,q=lambda q:block)
    def baseline(t):
        u=expm(-1j*t*h)
        return (u@rho@u.conj().T)[None],h[None]
    return c,baseline


def annihilators(n):
    result=[]
    for i in range(n):
        a=np.zeros((2**n,2**n),complex)
        for state in range(2**n):
            if state&(1<<i):
                a[state^(1<<i),state]=(-1)**((state&((1<<i)-1)).bit_count())
        result.append(a)
    return result


def fock_ensemble(rho):
    # Product probabilities in the actual natural modes; no rho/hole inverse.
    occ,u=np.linalg.eigh(rho)
    a=annihilators(len(occ))
    natural=[sum(u[k,i].conjugate()*a[k] for k in range(len(a))) for i in range(len(a))]
    eye=np.eye(2**len(a)); density=eye.copy()
    for p,mode in zip(occ,natural):
        number=mode.conj().T@mode
        density=density@(p*number+(1-p)*(eye-number))
    return density,a


def bilinear(coefficients,a):
    # sigma_ab is delta<c_b^dagger c_a>; coeff_ab contracts it bilinearly.
    return sum(coefficients[i,j]*(a[j].conj().T@a[i])
               for i in range(len(a)) for j in range(len(a)))
