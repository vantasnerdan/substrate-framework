# Exact reduction, thin-core bridge, and route verdicts

## 1. Full Euler relative-equilibrium equation

Let `u_0=u_epsilon` be the fixed Cao velocity, translating in the positive
`z` direction at speed `c_0`.  Put `R(x)=e_z cross x`.  For the physical
push-forward ansatz

    u(t,x)=Q_{Omega t} V(Q_{-Omega t}(x-c t e_z)),           (1)

direct differentiation gives

    partial_t u|_{t=0}=-c partial_z V-Omega[R,V].            (2)

Thus the exact stationary map is

    F(V,c,Omega)=P_L[(V dot grad)V]-c partial_z V
                   -Omega[R,V].                             (3)

At `(u_0,c_0,Omega)` the base is a trivial solution for every `Omega`, since
`[R,u_0]=0`.  With

    Gv=-P_L[((u_0-c_0 e_z) dot grad)v+(v dot grad)u_0],      (4)

the velocity derivative is

    D_V F(u_0,c_0,Omega)=-G-Omega ad_R.                     (5)

This derives the sign in the README from the laboratory solution, without
using Cao's unrelated stationary elliptic linearization.

Suppose, as a hypothesis to be earned, that in the reduced dynamically
accessible `l` sector

    Gq=i nu_l q,   [R,q]=i l q.                             (6)

Then (5) has a kernel at `Omega_0=-nu_l/l`.  On the real reflection slice,
the cosine/sine degeneracy becomes the one real kernel required by
Crandall--Rabinowitz.  If `q^*` spans the adjoint kernel, then

    partial_Omega D_V F q=-[R,q]=-i l q,                    (7)
    <q^*,partial_Omega D_V F q>=-i l <q^*,q>.               (8)

For an algebraically simple eigenvalue `<q^*,q>` is nonzero.  Therefore the
rotation-parameter transversality is automatic once the simple isolated mode
and Fredholm range have been constructed.  This removes transversality as an
independent unknown; it does not supply isolation or a range theorem.

## 2. KKS coefficient follows from the same eigenmode

Use the orbit convention `i_{X_J} Omega_KKS=dJ_z`.  Let real eigenvectors be
phased so that

    [R,q_c]=l q_s,
    sigma_l=Omega_KKS(q_c,q_s) != 0.                        (9)

For a prospective branch at orientation zero,

    omega_a=omega_0+a q_c+O(a^2).                           (10)

Its physical rotation tangent and amplitude tangent are

    X_J(omega_a)=-a l q_s+O(a^2),
    partial_a omega_a=q_c+O(a).                             (11)

Here the minus sign is fixed by the same physical push-forward used in (1):
`d/dtheta|(0) (R_theta)_* omega=-[R,omega]`.  It is not a freely chosen
Hamiltonian sign.

The moment-map identity gives, without needing the unknown second-order
field explicitly,

    dJ_z/da=Omega_KKS(X_J,partial_a omega_a)
            =l sigma_l a+O(a^2),                           (12)
    J_z(a)=J_z(0)+l sigma_l a^2/2+O(a^3).                  (13)

The centered no-swirl base has `J_z(0)=0`.  Consequently a nonzero
KKS/Krein signature of the actual mode is exactly the nonzero classical
angular-impulse coefficient required by the README.  The full group-cycle
action is `2*pi*J_z(a)`; an honest `C_l` quotient gives `2*pi*J_z(a)/l`.
This is a classical action identity and makes no quantization inference.

## 3. Route A: what CR proves and what the Cao source does not

The source orbital-stability theorem controls nonnegative **axisymmetric**
potential-vorticity rearrangements.  It neither evaluates the Hessian in
`l>=2` sectors nor constructs the full three-dimensional coadjoint chart.
Moreover, the steady Euler map in Sobolev velocity variables loses one
derivative, and its transport linearization is not automatically Fredholm
between an arbitrarily chosen `H^s` domain and range.

Equations (5)--(8) establish the exact kernel and transversality algebra.
Equations (9)--(13) establish the exact action consequence.  The remaining
Route-A inputs are specifically:

1. a simple isolated eigenvalue of the full Leray/Hodge operator (4) on the
   dynamically accessible, impulse-fixed, translation-quotiented sector;
2. a closed codimension-one range and complement inverse in a scale on which
   the nonlinear relative-equilibrium map is differentiable; and
3. an integrable volume-preserving displacement chart whose reconstruction
   retains compact vorticity and finite energy.

Cao's Appendix-C inverse concerns `-Delta-pU_+^(p-1)` in the axisymmetric
elliptic uniqueness proof.  It proves none of these three statements.  Route
A is therefore **blocked by the exact spectral/Fredholm construction**, not
refuted and not blocked by CR transversality.

## 4. Route B1: explicit filament Hamiltonian and positive Kelvin pair

The opened Fischer--Schopohl cutoff model writes a perturbed circular filament
as

    X(phi,t)=(R+r(phi,t)) e_r(phi)+z(phi,t)e_z.              (14)

For Fourier number `l>=2`, their equations (9)--(14) are

    dot z_l=B_l r_l,       dot r_l=-A_l z_l,                (15)

where, with `g=Gamma/(4*pi*R^2)`, `L=log(4R/xi_c)`, and
`S_l=sum_{j=1}^l 1/(2j-1)`,

    A_l/g=l^2(L-2S_l+1/2)+(3/2)S_l,
    B_l/g=(l^2-1)(L-2S_l+1/2)-(3/2)(S_l-1).                (16)

For fixed `l>=2` and sufficiently thin core both are positive, so

    nu_l=sqrt(A_l B_l)
        =|Gamma| l sqrt(l^2-1)L/(4*pi*R^2)
           *(1+O(1/L)).                                    (17)

This is an actual source-derived fixed-mode result, not a fitted dispersion
curve.  It also explains why the mode is a slow bending mode compared with
the `O(Gamma/xi_c^2)` internal core rotation.

The filament vorticity measure is

    omega=Gamma integral delta(x-X(phi)) X_phi dphi.        (18)

The finite vorticity-side rotational moment map gives

    J_z=-rho_0 Gamma/2 integral |X|^2 z_phi dphi
       =-rho_0 Gamma R integral r z_phi dphi+O(||(r,z)||^3).
                                                                    (19)

For the traveling mode

    z=a cos(l phi-nu_l t),
    r=(nu_l/B_l)a sin(l phi-nu_l t),                        (20)

equation (19) evaluates exactly to

    J_z=pi rho_0 Gamma R l (nu_l/B_l)a^2+O(a^3),            (21)

which is nonzero.  More explicitly, take the amplitude tangent `q_c` from
(20) at zero phase and set `q_s=(1/l) partial_phi q_c`, so that
`[R,q_c]=l q_s`.  Fischer--Schopohl's canonical variables `q=z`,
`p=rho_0 Gamma R r` then give

    sigma_l=Omega_KKS(q_c,q_s)
             =2*pi*rho_0*Gamma*R*nu_l/B_l,                 (21a)
    j_2=l*sigma_l/2
       =pi*rho_0*Gamma*R*l*nu_l/B_l.                       (21b)

This agrees in sign and magnitude with (21) and gives the positive quadratic
Hamiltonian when `A_l,B_l>0`.  Thus the thin-ring model contains the requested
period and nonzero classical action.

## 5. Route B2: map to the actual Cao core

The Cao blow-up is the compact radial column `W_p=U_+^p` from
`source-audit.md`.  Let `a_epsilon` be its physical core radius and set

    delta=a_epsilon/R,
    k_delta=l delta.                                        (22)

In the local Frenet frame, a displacement of the vortex center is the
cross-sectional column mode `m=1`, while `e^(i l phi)` around the ring becomes
the small longitudinal column wavenumber `k_delta`.  Gallay--Smets identify
precisely this mechanism: at `k=0` the `m=1` mode is the translation zero
mode, and for small nonzero `k` a simple imaginary Kelvin eigenvalue peels
away from the endpoint of the column's imaginary essential interval.  Its
nonanalytic small-`k` scale is the `k^2 log(1/k)` scale which, after restoring
the `Gamma/a_epsilon^2` core time, agrees with (17).

This produces a materially better bridge than transferring a hollow-core
mode directly: the model operator is now the exact full Euler column
linearization about the **source-defined Cao limiting profile**.  It also
exposes the threshold which a proof has to control.  Write `P_delta` for the
two real translation modes of the local core and `Q_delta=1-P_delta`.  The
required exact ring reduction is the Feshbach operator

    K_delta(z)=P(A_delta-z)P
      -P A_delta Q [Q(A_delta-z)Q]^-1 Q A_delta P.           (23)

To turn (17) into a Cao eigenvalue one must prove, on a contour of radius
comparable to the bending frequency,

    K_delta(z)=K_filament(z)+E_delta(z),
    sup_Gamma ||E_delta(z)||
       < inf_Gamma ||K_filament(z)^-1||^-1,                 (24)

and the corresponding complement resolvent and Riesz projection estimates.

A naive operator-norm comparison on the entire perturbation space is false:
the exact core contains `O(Gamma/a_epsilon^2)` internal modes and an essential
transport interval which the two-coordinate filament model does not contain.
Nor is a uniform zero-frequency inverse available: the translation mode
emerges from an endpoint of that essential interval, so the complement
resolvent degenerates as `k_delta->0`.  The estimate must be a **scaled
threshold/critical-layer estimate**, with the translation projection and the
`k^2 log(1/k)` singularity retained.  Cao's scalar elliptic inverse cannot be
substituted for (23)--(24).

The exterior part of the threshold singularity can nevertheless be evaluated
exactly.  Outside the compact column core the perturbation is irrotational,
so an `m=1`, axial-wave-number `k` potential has the decaying form

    phi(s,theta,z)=C K_1(|k|s) exp(i theta+i k z).           (24a)

Writing `x=|k|a_epsilon`, its boundary Dirichlet-to-Neumann multiplier is

    partial_s phi/phi
      =a_epsilon^-1 x K_1'(x)/K_1(x)
      =-a_epsilon^-1
        +k^2 a_epsilon[log(x/2)+gamma_E]
        +O(k^4 a_epsilon^3 log(x)^2).                       (24b)

This follows directly by substituting the decaying Frobenius jet

    K_1(x)=x^-1+(x/2)[log(x/2)+gamma_E-1/2]
             +O(x^3 log x)                                 (24c)

in the modified-Bessel equation.  Therefore the coefficient of the exterior
`k^2 log(1/k)` term is universal and already exact for the source-defined
compact Cao core.  What remains in (23)--(24) is not an unidentified exterior
tail: it is the finite core response, toroidal-curvature coupling of the two
translation polarizations, and a uniform complement estimate at the embedded
threshold.  In particular (24b) explains the logarithm in (17) while leaving
its order-one core/curvature matrix to the Feshbach calculation.

The profile applicability check is also explicit.  For

    F_p(s)=integral_0^s W_p(t)t dt=-s U'(s),
    Omega_p=F_p/s^2,
    D_p=2F_p-s^2 W_p,                                      (25)

Gallay--Smets' Richardson function becomes

    J_p(s)=2 F_p(s) W_p(s) s^4/D_p(s)^2.                   (26)

Their H2 hypothesis is exactly `J_p'(s)<0` in the core, together with an
admissible compact-edge approximation.  Equation (26), derived from the Cao
Lane--Emden equation, is the concrete scalar profile test.  It can be reduced
one step further without solving `U`.  Define

    y=-s U'/U,       a=s^2 U^(p-1),       t=a/y.              (27)

The Lane--Emden ODE gives

    s y'=y^2+a,
    s a'=a(2-(p-1)y),
    s t'=t(2-t-py),                                  (28)

and `D_p>0` is `0<t<2`.  Direct logarithmic differentiation yields the exact
sign criterion

    s J_p'/J_p
      =[8-t^2-2t-py(t+2)]/(2-t).                           (29)

The inequality in (29) can in fact be closed analytically.  Set

    h=p y-2(2-t).                                             (30)

The regular Lane--Emden expansion at the center, with
`lambda=U(0)^(p-1)`, gives

    y=lambda*s^2/2+(2-p)*lambda^2*s^4/16+O(s^6),
    t=2-p*lambda*s^2/4+p*(3p-4)*lambda^2*s^4/96+O(s^6),
    h=p*lambda^2*s^4/24+O(s^6)>0.                            (31)

Using (28),

    s h'=p y(y-t)+2t(2-t).                                  (32)

At any hypothetical zero of `h`, `p y=2(2-t)`, and the right-hand
side of (32) is `p y^2=4(2-t)^2/p>0`.  Since `h` starts positive, it cannot
cross from positive to nonpositive.  Hence `p y>2(2-t)` throughout the
core.  The numerator in (29) then satisfies

    (2-t)(t+4)-p y(t+2)
      <(2-t)[t+4-2(t+2)]=-t(2-t)<0.                          (33)

Thus `J_p'(s)<0` for every `0<s<1`.  At the edge,
`W_p=O((1-s)^p)` and `J_p` tends monotonically to zero; for `p>=2` the
zero extension is `C1`.  This proves the Cao limiting core's H2 sign and
puts it in the compact-profile closure contemplated by Gallay--Smets
(an explicit strict positive-tail approximation remains part of a fully
formal source-hypothesis map).  The hard unresolved input is now solely the
curved-ring threshold estimate (23)--(24), not the profile sign.

Route B therefore has two verdicts: the fixed-mode filament/action subroute is
**established within its declared model**, while transfer to the exact smooth
Cao ring is **blocked by the named threshold Feshbach estimate**.  The exact
Lane--Emden comparison above closes the limiting core's profile-H2 sign; the
remaining positive-tail source-hypothesis map is a regularization detail, not
the spectral construction.  This is a route gap, not a physical impossibility.

## 6. Route C: embedded-spectrum continuations

### C1. Reversible/KAM or Nash--Moser

If the exact ring mode is embedded rather than separated, ordinary CR fails
at its range hypothesis.  A reversible Lyapunov/KAM construction would need
a reducible normal operator, tame inverse estimates with explicit small
divisors, and a twist coefficient on the dynamically accessible leaf.  The
column endpoint described above is not a discrete normal frequency, and no
source theorem supplies those estimates for the toroidal Leray operator.
The exact formulation is viable only after the same threshold normal form in
(23), now with Cantor exclusions rather than a Riesz contour.  It is blocked
by that construction, not ruled out by continuous spectrum alone.

### C2. Exact helical contour branch and why it does not transfer directly

Cao--Fan--Li--Qin give an exact positive CR construction for 3D helical Euler.
Their contour multiplier is diagonal in Fourier modes; Proposition 3.5 proves
a one-dimensional kernel, a closed codimension-one range, and the explicit
transversality `-m sin(m theta)` outside that range.  This validates the
contour/CR mechanism as a real Euler construction rather than an abstract API.

It cannot itself meet 0040.  Their symmetry satisfies

    H_{2*pi}(x)=x+2*pi*h e_z.                               (34)

Any nonzero helical vorticity therefore repeats under every integer multiple
of `2*pi*h e_z`.  A compactly supported helical field on `R3` must be zero,
and a nonzero one has infinite total kinetic energy even when its energy per
pitch is finite.  This directly refutes only the **unmodified helical-to-Cao
transfer**.

A toroidal contour replacement encounters a second exact difference.  The
helical no-swirl symmetry reduces vorticity to one advected scalar and closes
the boundary equation.  A nonaxisymmetric vortex ring has stretched vector
vorticity; its boundary shape alone does not determine the interior Cauchy
invariants.  An honest toroidal route must evolve a volume-preserving solid-
torus map (boundary plus internal labels), derive its full Biot--Savart/Leray
functional, and prove a Fourier-block range theorem.  The patch construction
is therefore a useful blueprint, but not the missing finite-energy branch.

Route C's direct helical transfer is **refuted by (34)**; the reversible and
solid-torus representation changes are **blocked by their named normal-form
and vector-vorticity closure constructions**.  No global nonexistence claim is
made.

## 7. Strongest exact result and next dependency

The activated attempt does not construct the registered exact Cao branch.
It does establish three pieces that materially narrow the construction:

1. once a simple isolated signed mode exists, the physical rotation parameter
   automatically satisfies CR transversality, equations (7)--(8);
2. the same nonzero KKS signature automatically gives the nonzero quadratic
   angular-impulse/action coefficient, equations (12)--(13); and
3. the source-derived thin-ring Kelvin pair has the explicit positive
   frequency (17) and moment (21), and maps to the `m=1`,
   `k=l a_epsilon/R` threshold branch of the actual Cao limiting column; and
4. the Cao Lane--Emden profile satisfies the exact Richardson sign, while the
   compact-core exterior contributes the universal Dirichlet-to-Neumann term
   `k^2 a_epsilon[log(|k|a_epsilon/2)+gamma_E]`.

The next executable dependency is a **Cao-ring threshold spectral theorem**:
construct the translation-bundle Feshbach map (23) and establish the scaled
critical-layer resolvent/Riesz estimate (24) for the curved full Leray
operator.  The source-defined core-profile sign is already closed by
(30)--(33), and (24a)--(24c) close the exterior logarithmic coefficient.  The
unresolved matrix entries are the finite-core response and toroidal-curvature
coupling, together with the uniform complement bound.  With that theorem,
Route A's transversality and action gates are also closed; the following step
is a tame volume-preserving solid-torus CR/Nash--Moser reconstruction.  Until
those two constructions exist, LP2/P4 and the particle objective remain
active.
