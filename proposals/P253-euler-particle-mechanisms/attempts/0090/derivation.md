# Comoving Maxwell shell-to-flux derivation

## 1. Typed current, Fourier signs, and the outgoing field

Use

    f_hat(k)=(2 pi)^(-3/2) integral exp(-i k dot x) f(x) dx

and write the real source as

    (rho,J)(t,x)=Re{exp(-i omega t)(rho_omega,J_omega)(y)},
    y=x-c_g t e_z.                                            (1)

Its laboratory temporal frequency at spatial wave vector `k` is

    Omega(k)=omega+c_g k_z.                                   (2)

The continuity equation becomes

    -i Omega rho_hat+i k dot J_hat=0,
    Omega rho_hat=k dot J_hat.                                (3)

Assume either explicitly `integral rho_omega dx=0`, or that `rho_omega` and
the profile flux are integrable so their Fourier transforms are continuous at
zero and (3) may be evaluated there.  For `omega!=0`, (3) at `k=0` then gives
`rho_hat(0)=0`.  The oscillatory mode is net neutral.  The nonzero charge of a
charged carrier is an `omega=0` affine background and is not part of this
shell calculation.  The later assumption `<x>J_omega in L2` supplies the shell
trace; by itself it does not supply this pointwise zero-frequency row.

In Coulomb gauge, or equivalently after transverse projection in Lorenz
gauge, the vector potential satisfies

    [partial_t^2-c_EM^2 Delta]A_T=mu_EM c_EM^2 J_T.            (4)

For the convention (1), the retarded/Sommerfeld boundary value is

    A_T_hat=mu_EM c_EM^2
      [D_omega-i0 sign Omega]^-1 P_T J_hat,
    D_omega=c_EM^2|k|^2-Omega^2.                              (5)

The subluminal condition implies `sign Omega=sign omega` everywhere on a
real shell component: on the shell `|Omega|=c_EM|k|`, while
`|c_g k_z|<c_EM|k|`.  Because `E_T=-partial_t A_T`,

    E_T_hat=i Omega A_T_hat.                                  (6)

The longitudinal field is fixed by Gauss and (3).  It contributes reversible
Coulomb/near-field energy, but its numerator does not occur in the transverse
outgoing flux.

## 2. Exact source-work sign and mean outgoing power

Poynting's theorem is

    partial_t u_EM+div S=-J dot E,
    u_EM=epsilon_EM|E|^2/2+|B|^2/(2mu_EM),
    S=E cross B/mu_EM.                                        (7)

Thus positive outward power supplied by a stationary harmonic source is

    P_out=-(1/2)Re integral J_hat dot conjugate(E_T_hat) dk.   (8)

Now

    conjugate(E_T_hat)
      =-i Omega mu_EM c_EM^2
        [D_omega+i0 sign Omega]^-1 conjugate(P_T J_hat),

and

    [D+i0 s]^-1=PV(1/D)-i pi s delta(D).                      (9)

The principal-value contribution in (8) is reactive.  The real delta term in
`J dot conjugate(E)` is

    -pi mu_EM c_EM^2 |Omega| delta(D_omega)|P_T J_hat|^2.

Therefore

    P_out=(pi mu_EM c_EM^2/2)
      integral |Omega| delta(D_omega)|P_T J_hat|^2 dk
      =pi/(2 epsilon_EM)
      integral |Omega| delta(D_omega)|P_T J_hat|^2 dk.        (10)

The minus sign in (8), complex conjugation, and outgoing `-i0 sign Omega`
are all load bearing.  Reversing only one of them would give the wrong sign.

The phasor factor `1/2` in (8) already accounts for the conjugate real-field
component.  When the complex `+omega` profile is integrated over all `k`, no
second identical shell is added.  Equivalently, a two-frequency calculation
assigns half-amplitude to each conjugate component and their sum reproduces
(10).

## 3. Shell geometry and coarea constants

Put `k=r n`, `|n|=1`, `s=sign omega`.  The positive-radius shell satisfies

    r(n)=|omega|/[c_EM-s c_g n_z],
    Omega(rn)=s c_EM r.                                       (11)

Direct differentiation gives

    grad_k D_omega=2[c_EM^2 k-c_g Omega e_z],                 (12)

and on (11)

    |grad D|^2=4 c_EM^2 r^2
       [c_EM^2+c_g^2-2s c_EM c_g n_z].                        (13)

This is bounded below by
`4 c_EM^2 r^2(c_EM-|c_g|)^2`, so the shell is smooth for every fixed
subluminal margin.  The radial derivative is especially simple:

    |partial_r D|=2|c_EM^2 r-c_g n_z Omega|
                 =2 c_EM |omega|.                            (14)

Coarea in star-shaped coordinates turns (10) into

    P_out=pi mu_EM c_EM^2/(4|omega|)
      integral_S2 r(n)^3 |P_T J_hat(r(n)n)|^2 dn              (15)

or

    P_out=pi |omega|^2/(4 epsilon_EM)
      integral_S2
        |P_T J_hat(r(n)n)|^2
        /[c_EM-s c_g n_z]^3 dn.                              (16)

Assume explicitly that `<x> J_omega` lies in `L2`, equivalently that
`J_omega` and every `x_j J_omega` lie in `L2`.  Then `J_hat_omega` is `H1`
near the compact shell, and the trace theorem restricts it continuously to
that shell.  Smooth compact currents satisfy this hypothesis.  Since
`omega!=0`, the shell stays away from `k=0`, where `P_T(k)` is singular.  These
facts also give the standard weighted outgoing resolvent boundary value
locally away from threshold.  More explicitly, with

    alpha=c_EM^2-c_g^2>0,

completion of the square gives

    D_omega=c_EM^2|k_perp|^2
      +alpha[k_z-omega c_g/alpha]^2
      -omega^2 c_EM^2/alpha.                                (16a)

The translation in `k_z` is a spatial modulation and the positive anisotropic
rescaling `p_z=sqrt(alpha)[k_z-omega c_g/alpha]/c_EM` reduces (16a) to the
standard Helmholtz symbol.  On a fixed subluminal margin these maps and their
inverses are bounded on weighted spaces, so the standard outgoing
`L2_1 -> L2_-1` Helmholtz limiting-absorption estimate supplies the
anisotropic bridge used here.

For an exposing exact source set `c_g=0`, `rho_omega=0`, and

    J_hat(k)=i j_0 (k cross a) exp(-sigma^2|k|^2/2).           (17)

It is a smooth divergence-free Gaussian curl current.  At
`r_0=|omega|/c_EM`,

    integral_S2 |k cross a|^2 dn=(8 pi/3)r_0^2|a|^2.

Equation (16) gives exactly

    P_G=2 pi^2 j_0^2 |a|^2 |omega|^4
          exp[-sigma^2 omega^2/c_EM^2]
          /(3 epsilon_EM c_EM^5).                            (18)

This source checks the sign and every normalization constant without a
remaining angular quadrature.

## 4. Conserved switching and an operational finite-time observable

Multiplying both profiles in (1) by a real envelope `a_T(t)` produces the
continuity defect `a_T' rho_omega`.  Since the oscillatory charge has zero
mean, a compact-domain Bogovskii construction, or the corresponding weighted
Hodge construction, gives `K_omega` with

    div K_omega=-rho_omega.                                   (19)

Here `rho_omega` is required to belong to the declared zero-mean
Bogovskii/Hodge source class, and `P_T K_omega` belongs to the weighted
energy-source class used in (22)--(24).  Smooth compact sources satisfy both
requirements.

The complex source

    rho_T=a_T exp(-i omega t)rho_omega(y),
    J_T=exp(-i omega t)[a_T J_omega(y)+a_T' K_omega(y)]       (20)

is exactly conserved.  In Fourier space, with
`A_T(xi)=integral exp(i xi t)a_T(t)dt`, integration by parts gives the
transverse switching numerator

    A_T(xi) P_T[J_hat_omega-i xi K_hat_omega],
    xi=c_EM|k|-Omega(k).                                     (21)

In particular, the completion vanishes at the carrier shell `xi=0`; it does
not alter the coefficient (10), but it contributes switching radiation away
from that shell.

For a compactly supported smooth envelope and zero incoming radiation, the
late-time outgoing transverse free-field energy, after the source switches
off and after subtracting the static charged background, is exactly the
oscillator spectral norm

    F_T(k)=integral_R exp(i c_EM|k| t)
                    P_T(k) J_T^real(t,k) dt,
    E_late=1/(2 epsilon_EM) integral |F_T(k)|^2 dk,            (22)

where the spatial Fourier transform inside `J_T^real(t,k)` is the unitary
convention of Section 1.  Thus `F_T` is defined once from the full real current.
Writing that real current as its two conjugate phasors merely decomposes this
same `F_T` into one half the sum of their expressions (21); it does not define
a second shell or add another copy of the power.  Formula (22) retains their
finite-window interference.

Here is the precise plateau theorem.  First isolate the rectangular resonant
phasor and define its one-dimensional pushforward density

    H(xi)=1/(8 epsilon_EM) integral delta(xi-[c_EM|k|-Omega(k)])
                         |P_T J_hat_omega(k)|^2 dk.           (23)

The factor `1/8` is the late-energy factor `1/(2 epsilon_EM)` times the
real-phasor amplitude squared `1/4`.  With

    I_T(xi)=integral_0^T exp(i xi t)dt,
    C(tau)=integral_R exp(i xi tau)H(xi)dxi,                  (23a)

the resonant rectangle energy and shell power are exactly

    E_rect(T)=integral |I_T(xi)|^2 H(xi)dxi,
    P_out=2 pi H(0).                                         (23b)

The autocorrelation identity

    |I_T(xi)|^2=integral_(-T)^T
                  (T-|tau|)exp(i xi tau)d tau                (23c)

therefore gives

    E_rect-T P_out
      =-integral_(|tau|<T)|tau|C(tau)d tau
       -T integral_(|tau|>=T)C(tau)d tau.                    (23d)

Thus the explicit first-moment hypothesis

    M_rect=integral_R |tau| |C(tau)|d tau < infinity         (23e)

implies `|E_rect-T P_out|<=2 M_rect`, with the exact coefficient in (10).

For a fixed smooth ramp class, its Fourier source differs from the rectangle
by fixed endpoint profiles

    R_-(xi,k)+exp(i xi T)R_+(xi,k).                           (23f)

These profiles include the entire conserved `a_T'K` completion: the factor
`-i xi I_T(xi)K=(1-exp(i xi T))K` is itself an endpoint pair.  Define the
rectangle--endpoint pushforward correlations by replacing one factor in
(23) with `R_-` or `R_+`.  Assume their inverse Fourier transforms are `L1`,
and assume the endpoint profiles have finite late-field energy.  Then

    |integral I_T(xi) H_(0,-)(xi)dxi|
      <=||C_(0,-)||_L1,

with the same bound after translating the right endpoint.  Endpoint--endpoint,
`J-K`, and `K-K` terms are bounded by their fixed energy norms.  This types
the phrase "fixed ramps" rather than hiding their contribution in an
approximate delta.

For `omega>0`, the conjugate temporal phasor has positive-free-frequency
detuning

    xi_-(k)=c_EM|k|+omega-c_g k_z
            >=omega+(c_EM-|c_g|)|k|>=omega.                 (23g)

Hence `|I_T(xi_-)|<=2/omega`; its self energy is uniformly bounded by the
corresponding fixed weighted source norm.  The cross term with the resonant
branch has the form

    integral_0^T integral_0^T
      C_+-(t-s)exp(-2i omega s)dt ds.                         (23h)

If `integral (1+|tau|)|C_+-(tau)|d tau` is finite, its constant interior
part is bounded by `1/|omega|` and both endpoint differences are bounded by
the first moment.  Negative `omega` exchanges the two conjugate roles.

Collecting (23e)--(23h), define `M_1` to be the sum of the rectangle first
moment, the endpoint cross-correlation `L1` norms, endpoint energies,
opposite-branch weighted source norm, and conjugate cross first moment.  For
this full-real conserved source,

    E_late(T)=T P_out+E_switch(T),
    |E_switch(T)|<=C(omega,ramp) M_1.                         (24)

The bound is independent of the plateau length.  Without these explicitly
listed correlation and endpoint norms, (22) remains exact and the supported
verdict is the spectral formula rather than an asserted `O(1)` remainder.

A sharp characteristic envelope has distributional endpoint currents
`a_T' K`.  If `P_T K` lies in the energy source space, the impulses may create
a finite field-velocity jump; infinite endpoint energy is not automatic.
The sharp source belongs to a weak/distributional class and its endpoint
energy must be computed as a limit of declared smooth ramps before comparison
with (24).

The same operational observable follows from a moving control volume
`V(t)=V_0+c_g t e_z`.  Reynolds transport and (7) give

    d/dt integral_V u_EM
      +integral_boundary(V) (S-c_g u_EM e_z) dot n
      =-integral_V J dot E.                                  (25)

Thus, for every finite moving control volume, relative boundary flux plus the
change of stored energy equals source work.  A fixed bounded laboratory sphere
does not contain a source translating distance `O(T)` and is not the
observable used in (24).  Equation (22) defines the late-free-field `E_rad`;
(25) is an exact finite-volume balance.  Equality between its infinite-world-
tube limit and (22), including endpoint stored energy, requires a separate
decay/limit argument and is not claimed here.

## 5. Reciprocal one-block conversion and its boundary

Suppose a normalized positive-frequency carrier mode has physical action
`A_mode |z|^2`, energy `nu_phys A_mode |z|^2`, and current

    J=z j+conjugate(z j).                                     (26)

Only when (26) and the mode equation arise from the same joint action does
eliminating the retarded Maxwell field produce a legitimate self-energy.
Its boundary imaginary part is the quadratic shell functional (10), with
the sign fixed by (8).  If `P_out[j]` is the power at `|z|=1`, energy balance
for an amplitude law `dot z=(-i nu_phys-gamma)z+...` gives

    2 gamma nu_phys A_mode=P_out[j],
    gamma=P_out[j]/(2 nu_phys A_mode).                        (27)

Hence the leading energy/action loss fraction over a gate is controlled by

    T_gate P_out[j]/(nu_phys A_mode),                         (28)

not by power divided by action.  Equation (27) is an identity conditional on
a reciprocal closed amplitude equation.  The prescribed-current Maxwell
calculation alone does not construct the pole, establish analytic Fredholm
continuation, determine its real shift, or prove a nonzero source-specific
current trace.  Exact vanishing of the transverse trace is a BIC candidate;
an exponentially small but nonzero trace is a quasimode/radiation-leakage
input rather than an `L2` eigenmode.

## 6. Route verdicts and next achievements

- Route A is established for prescribed conserved currents with the stated
  weighted trace regularity: (10), (15), and (16) are the exact outgoing
  power identities, and (18) is an exact source-derived oracle.
- Route B is established at the exact conserved-source and late-energy levels,
  (19)--(22), and yields (24) under the explicit full-real pushforward,
  correlation, endpoint, and fixed-ramp hypotheses (23)--(23h).  Sharp
  switching is a separately typed weak limit with endpoint energy retained;
  (25) is a finite moving-volume identity rather than an asserted equivalence
  with the late-field limit.
- Route C establishes the conditional reciprocal conversion (27)--(28).  A
  Cao damping rate or resonance remains open on the actual current map,
  physical KKS normalization, and one-block analytic continuation.

These results supply classical radiation conversion only.  They do not select
a universal action, Born weights, reset, exchange/fermionic character, an
electron, or a neutrino.  In the admitted U(1) extension,
`q^2/(4 pi epsilon_EM c_EM)` has action units and the same sector supplies a
causal speed.  Turning that dimensional combination into a selected quantum
requires two new achievements: an integral compact gauge character/topological
tag fixing `q/g`, and same-carrier dynamics selecting `A_mode` in proportion
to it.  Classical radiation alone favors zero excitation and supplies neither
selection theorem.  The bare-Euler routes therefore remain active in parallel.
