# Actual one-wave physical mean: complete fixed-time second spatial jet

## 1. Same Euler field, actual preparation and whole-field law

Fix a unit vector n, orthonormal e1,e2 perpendicular to n, q>0 and v>0.
The field

    u(x)=v[e1 cos(q n·x)+e2 sin(q n·x)]

is smooth, bounded, divergence free and stationary Euler, with constant
pressure: u·∇u=0. It is a constant-curl Beltrami field, with the curl sign
fixed by the orientation of (e1,e2,n). Its actual cell covariance is

    C=〈u⊗u〉=c Pn,   c=v²/2,   Pn=I−n⊗n.

Translate the whole field by a uniform phase and rotate the whole field
by one Haar SO(3) variable. Whole-field time reversal can also be paired;
it changes none of the second-order coefficients below. No independent
cell rotations or projected/rescaled SH velocities are used.

Put K=k κ, |κ|=1, and prescribe D,V perpendicular to κ. The two actual
linear Euler initial data, combined by superposition, are

    w_D(0)=P_K[D exp(iK·x) × curl u],
    w_V(0)=V exp(iK·x).

The first is a solenoidal Kelvin displacement preparation. The second
is the independently allowed full Euler/Lin velocity preparation; its
kinetic mass is exactly rho. Each nonlinear trajectory conserves its
own circulation. Their initial circulation classes need not coincide.
Let m be the full physical slow Fourier mean of w and define
X(t)=D+∫_0^t m(r)dr. In particular m(0)=V exactly.

## 2. The cell operator actually reduces

Use the exact material equations and current of 0144, not a supplied
oscillator. Write η=U+i k χ+O(k²), 〈χ〉=0, div χ=0. All first cell
functions depend only on q n·x. Thus n·χ=0, u·∇χ=0, and Hess p=0.
The actual zero-mean solenoidal cell equation is exactly

    χ_tt=−2 a Pn U_t,               a=κ·u.                 (1)

The full pressure is still present: Pn is the microscopic pressure
projector, and Pκ is the separate slow projector. The initial rows are

    χ(0)=0,
    χ_t(0)=S_D=−Pn[a D+κ(u·D)].                            (2)

Equation (2) follows by expanding the actual material Kelvin operator
−P_K[(u·∇)η+(Dη)^T u]; it is not chosen to obtain the desired sign.
The independent common-V phase has zero initial cell position and rate.
At the retained order U=D+tV and hence

    χ(t)=t S_D−t² a Pn V.                                 (3)

The physical current and stress are

    m=U_t−k² Pκ〈aχ−u(κ·χ)〉+O_T(k³),
    m_t=k²{B U+Jκ χ_t}+O_T(k³),
    B=〈a²〉=c(1−s²), s=κ·n,
    Jκ f=Pκ〈a f+u(κ·f)〉.                                (4)

These retain the complete ambient-fluid momentum. The microscopic
mean material rate is not silently substituted for m.

## 3. Both initial phases give the same whole second jet

Direct covariance contraction gives, on the transverse macro plane,

    B D+Jκ S_D = −C_n D,
    B V−2 Jκ(a Pn V) = −C_n V,

    C_n = c(1−s²) I_T
          −2c(1−2s²)(Pκ n)⊗(Pκ n).                       (5)

The unobserved cell vectors S_D and −2a Pn D are generally different.
Their difference is annihilated by Jκ for this circular covariance.
It is this actual observation identity that removes the first-cell
memory from the physical second spatial jet. It does not identify the
cell vectors or discard the physical current.

Consequently, for every fixed T,

    X(t)=D+tV−k² C_n(t²D/2+t³V/6)+O_T(k³),
    X_tt+k² C_n X=O_T(k³).                                (6)

For SH polarization parallel to n×κ and SV polarization parallel to
Pκ n, the two coefficients are respectively

    c_SH²=c(1−s²),
    c_SV²=c(1−s²)(4s²−1).                                (7)

At |s|=1 both vanish and the polarization limit is unambiguous. The SV
coefficient is negative for |s|<1/2. At s=0 it is −v²/2, reproducing
0129's genuine unstable Euler acoustic branch rather than contradicting
it. No orientation has been selected by a favorable sign.

For a fixed transverse unit vector e, Haar averaging gives

    E[(n·e)²]=1/3,   E[s²(n·e)²]=1/15,
    E[C_n]=c[2/3−2(1/3−2/15)] I_T
           = (2v²/15) I_T.                               (8)

All realizations have the same leading D+tV; averaging (6) therefore
gives the actual full physical mean, on both transverse polarizations,

    Xbar_tt+k² c_eff² Xbar=O_T(k³),
    c_eff²=2v²/15>0.                                     (9)

The density is rho, not the 2rho obtained by reconstructing two copies
of an SH-only population. A formal local quadratic action reproducing
(9) is rho(|Xbar_t|²−k² c_eff²|Xbar|²)/2. This observation-level
representation is not being promoted as the complete microscopic
action pullback for arbitrary histories; the actual action is retained
in section 5.

0132's equal-amplitude ABC axes have isotropic second covariance but
not the whole Haar fourth-rank orientation law. For κ along the third
axis its two orthogonal helical components give opposite even second
jets on a given transverse polarization; their cancellation is exactly
compatible with (5). The chiral cubic response of that source is not
changed here. No magnetic energy is imported into (9).

## 4. Uniform fixed-window error from the full, untruncated Euler operator

The one-coordinate structure permits a direct remainder estimate,
without discretization or an assumed cell spectral gap. On periodic
L² vector functions of θ=q n·x, the exact Bloch Euler operator is

    L_k w=−P_k[i k(κ·u)w+q(n·w)u_θ].                     (10)

Here the projector on harmonic ell is P_(ell q n+kκ); on the slow
harmonic it is the fixed Pκ. Equation (10) includes all Fourier modes
and is bounded: there is no missing u·∂θ transport, since u·n=0.
For complex δ=k/q with |δ|≤1/4, every nonzero-harmonic projector is
holomorphic and uniformly bounded. Indeed the denominator divided by
ell²q² differs from 1 by at most 2|δ|/|ell|+|δ|²/ell²≤9/16.
The numerator is uniformly bounded as well. Thus universal constants
C0,C1 bound ||L_k||≤C0 qv, uniformly in orientation and phase.

The two initial columns are holomorphic with respective L² bounds
C1 qv|D| and |V|. The exact solution exp(t L_k) and its mean are
therefore holomorphic, with norm bounded by C1 exp(C0 qv t).
The Cauchy estimate for the Taylor remainder on any smaller disk,
for example |δ|≤1/8, proves (6) and (9) uniformly over the whole law.
Writing τ=qvT and A=|D|+|V|/(qv), one convenient residual bound is

    sup_[0,T] |Xbar_tt+k² c_eff² Xbar|
        ≤ C2 q²v² |δ|³ (1+τ) exp(C3τ) A.              (11)

The harmless k² times the second-order displacement correction is
absorbed in (11). Dividing by the actual positive natural scale
k² c_eff² yields

    residual/(k² c_eff²)
        ≤ C4 |δ| (1+τ) exp(C3τ) A.                     (12)

Thus the optical/microscopic time window can be selected first and the
relative acoustic second-jet error then made arbitrarily small by k/q.
The constants do not conceal a discretization, box inversion floor or
an orientation-dependent small denominator.

This estimate does not control T proportional to 1/|k|. The retained
negative SV sector and 0129's actual growing branches are concrete
reasons not to substitute an orientation-averaged frequency at that
scale. Even an ideal collection of already closed orientation waves
would generally retain higher moments E[C_n^r], not powers of E[C_n].
The scalar diagnostic on a fixed transverse e is

    E|C_n e|²=136c²/315,
    E|C_n e|²−|E C_n e|²=568c²/1575>0;

its strict inequality is checked directly below, without claiming that
this diagnostic alone computes the actual acoustic-time fourth jet.
An acoustic-time mean normal form remains a different construction.

## 5. Actual action and material current retained for the coupled join

For this field the full pressure-constrained Euler/Lin action is

    L=rho/2〈|(∂t+i k a)η|²〉,   div_K η=0.

Expanding it through second spatial order, and retaining the exact
time boundary rho k² d〈a U·χ〉/dt, gives

    L=rho|U_t|²/2
       +rho k²/2[〈|χ_t|²〉+B|U|²−4〈χ·a U_t〉]
       +rho k² d〈a U·χ〉/dt
       +O(k³).                                         (13)

The mean X is related to U by the current in (4). The microscopic cell
state (3), its endpoint data, and this observation map are part of the
action join. An on-shell positive force coefficient does not authorize
varying a prepared on-shell cell history as though it were an arbitrary
off-shell history with identical endpoint data. This is precisely why
the formal action after (9) is labeled observation-level.

Actual circulation preservation, mass rho, the full action (13), and
the positive closed physical second jet (9) are established together.
The intrinsic material angle, spin normalization and same-field
translation/rotation action cross terms remain the parent's active
construction; they are not supplied by renaming a wave phase.
