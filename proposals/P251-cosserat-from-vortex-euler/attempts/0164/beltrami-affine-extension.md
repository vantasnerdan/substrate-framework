# Exact constant-lambda extension of the acoustic affine response

This continuation uses the new exact identity communicated by the parent;
the signs, domain and full Lin transport are derived here. It replaces
the missing low-frequency pressure estimate in section5 of the earlier
affine proof for the actual constant-lambda EPS/periodic fields. It does
not replace the exact 2D3C identity or change a frozen earlier attempt.

## 1. The Leray term is exactly local on the nonzero Helmholtz shell

Let curl u=lambda u, lambda!=0, div u=0, and A x=n cross x. For the
actual periodic constant-lambda field, or the bounded Fourier/Herglotz
realization on its nonzero Helmholtz shell,

    curl(Au)=-partial_n u,
    curl[-partial_n u/lambda]=-partial_n u.

The difference is curl free. There is no zero Fourier harmonic: the
field lies on |g|=|lambda|. Thus its full Leray decomposition gives

    P(Au)=-partial_n u/lambda.                           (1)

This uses the actual pressure projector, not a local truncation of it.
There is no unmeasured large-domain Riesz tail. The multiplier identity
also applies to the affine derivatives of shell-supported tempered
fields: the projector is smooth near that shell. A nonzero imposed
harmonic velocity would need its separate harmonic term and is not
silently included in(1).

## 2. Exact affine velocity history and local material correction

Let L w=-P[(u dot grad)w+(w dot grad)u] be the actual linear Euler
operator. The static rotational symmetry tangent is

    r=A u-Du A x.

Both r and partial_n u are stationary Euler symmetry tangents, hence
Lr=L(partial_n u)=0. Equation(1) yields

    L(Ax)=r+2 partial_n u/lambda,
    w_V(t)=Ax+t[r+2 partial_n u/lambda].                 (2)

As before this is the covering-space AFFINE DERIVATIVE of the actual
Bloch family, not a finite-energy periodic rigid velocity. The demodulation
can be written explicitly. The zero-K velocity histories are
w_D=-Du D and w_V=V-t Du V. In the rotational affine sector their
polynomial modulation pieces are -Du Ax and Ax-t Du Ax. Subtracting
them from the full histories leaves the PERIODIC first velocity-cell
rows Au and t(Au+2 partial_n u/lambda), respectively. The corresponding
material demodulation subtracts Ax and t Ax, leaving0 and zeta below.
Thus a derivative of a moving lattice is not substituted for a fixed-cell
Bloch corrector.

Put B_L=u dot grad-Du. The material equation is eta_t+B_L eta=w.
Since B_L(Ax)=r and B_L n=-partial_n u, the exact displacement is

    eta_V=t A x+zeta(t),
    zeta=(2/lambda)[-t n+integral_0^t exp(-s B_L)n ds]
         =(2/lambda)integral_0^t(t-s)exp(-s B_L)partial_n u ds. (3)

Direct differentiation of(3) gives zeta_t+B_L zeta=
2t partial_n u/lambda, including its sign and initial zero displacement/
rate. Thus(2)-(3) are actual Euler AND actual material histories.

The propagator in(3) is ordinary vector-field transport by the local
flow: exp(-s B_L)f(x)=Dg_s(g_(-s)x) f(g_(-s)x).
On the material tube swept out by the tag and its backward trajectories,
the usual differentiated flow equation therefore gives

    ||zeta||_(C^r,tag;[0,T])
       <= C_(r,T,local) T² ||partial_n u||_(C^r,tube)/|lambda|. (4)

The constant uses the needed LOCAL flow derivatives, not a global norm
growing with the torus radius. Apply the actual central angle, spin,
displacement and shape rows to this correction. Their finite operator
norms include the small reference quadrupole and the packet's own scale.
The corresponding time derivatives follow from the same local transport
equation, rather than assuming that a small position error controls rates.

## 3. The actual EPS geometry can make this correction small

For the selected large-radius constant-lambda toroidal field, n is the
fixed tangent at the packet center. On a fixed finite packet arc and its
finite-time material neighborhood,

    ||partial_n u||_(C^r,tube)
        <= C_local (1+L_arc+U_local T)/R

with the dimensionful local coefficient scales included in C_local.
This follows by differentiating the rotating cylindrical basis and the
fixed-tangent/actual-tangent difference along the arc. All packet length,
tag width, observation conditioning and time scales are fixed before R
is chosen. A subsequent sufficiently accurate actual periodic/Herglotz
transfer preserves this C^r bound on the same swept tube. The full field
remains constant-lambda, so(1) continues to hold exactly.

Consequently the Euclidean physical angle has the affine rows

    theta_D=alpha,
    theta_V=t alpha+E_n(t)alpha,
    ||Omega E_n||_(time derivatives scaled by Omega)
                  <= epsilon_R,

where epsilon_R can be selected after the packet's own normalization.
Whole-field O(2)/mirror averaging removes the other affine tensor sectors
exactly. The physical mean/angle bracket is therefore

    {X,theta}=curl E(t)/(2rho)+O_T(K³),                  (5)

with its actual retained rate-rate phase form and time connections.
The Euclidean observer has exact static unit response; for a calibrated
observer also retain its known R_E(t) factor and measured correction.

## 4. What this extension earns

The symmetry route now supplies the ACTUAL same-EPS first-jet material
response with arbitrarily small controlled realization error. The error
is compared at the physical mixed scale and, where needed, at the optical
curvature scale of the companion current proof. It is not a demand that
the coefficient vanish at every nonzero K. The full joint action and
all current/shape connections remain the ones derived in0158/0162.

The previous local-field/low-frequency-pressure gap is closed for this
constant-lambda route by(1), not by asserting that every Leray projector
is local. This does not prove a general non-Beltrami transfer theorem,
the same-field acoustic scalar closure for every geometry, or a uniform
large-cell homogenization limit. Those separate parent obligations keep
their original scope.
