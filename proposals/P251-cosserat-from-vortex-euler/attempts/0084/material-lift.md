# Actual material reconstruction and a local collar/action repair

## 1. Reconstruct the transported set, not its instantaneous orbit generator

Let `u=u0` be the selected smooth stationary Euler field, `D` its bounded
invariant solid torus, `Sigma=boundary D`, and `chi0=1_D`. Thus `u.n=0` and
`u.grad chi0=0` distributionally. For a smooth divergence-free Eulerian
velocity perturbation `v`, actual material displacement satisfies the Lin
identity

    xi_t+[u,xi]=v,    [u,xi]=(u.grad)xi-(xi.grad)u.

Writing `Phi_t` for the background flow, its explicit solution is

    xi(t)=(Phi_t)_*xi_initial
            + integral_0^t (Phi_(t-s))_*v(s) ds.

Smoothness and uniqueness hold on any finite time interval where the
background flow and coefficients are smooth on the transported support.
Pushforward by this volume-preserving flow preserves divergence-free fields.
This is a genuine construction, including the initial material tag. Define

    delta chi=-xi.grad chi0.

Differentiating gives exactly
`(partial_t+u.grad)delta chi=-v.grad chi0`. Thus a compact instantaneous
orbit generator with zero boundary trace is not a substitute for this `xi`.
Nonzero normal Leray reaction drives the reconstructed boundary even when
the instantaneous generator vanishes there. Keeping this transported tag
state is already an exact local-in-time extended description; eliminating
it generally introduces history, as analyzed separately in 0082.

## 2. The full moving-boundary spin becomes local in material displacement

Put `M=rho|D|`, `X=M^-1 integral_D rho x`, and `r=x-X`. Invariance gives
`integral_D u=integral_Sigma x(u.n)=0`. The centroid perturbation is
`delta X=M^-1 integral_D rho xi`, with its actual time derivative.
Differentiate the material definition of intrinsic angular momentum:

    delta S = rho integral_D [xi cross u+r cross Dt xi],
    Dt xi=xi_t+(u.grad)xi.

Terms involving `delta X` and `delta Xdot` vanish by the two centered
integrals, not by assuming a rigid parcel. Integration by parts uses the
actual invariant boundary condition `u.n=0` and gives the useful identity

    delta S = rho integral_D [r cross xi_t+2 xi cross u].       (1)

Equivalently, in Eulerian variables,

    delta S = rho integral_D r cross v
                  +rho integral_Sigma (r cross u)(xi.n).       (2)

Equation (2) includes the full moving-boundary term missing from an
instantaneous velocity moment. Equations (1) and (2) are equal because
`v=Dt xi-(xi.grad)u`. Background invariant domains have stationary centroids;
no pointwise ensemble mean has replaced them.

For a material-first slow-affine restriction `xi=sum_a Xi_a z_a(t)`, the
physical spin is the explicit local map

    delta S=A zdot+C z,
    A_a=rho integral_D r cross Xi_a,
    C_a=2rho integral_D Xi_a cross u.                           (3)

It includes all collar, tube and selected boundary motions. Matched lifts
on a time-reversed pair cancel `C` but retain `A`; this is an identity of
material observables, not a license to add geometric inertia to an already
Kelvin-reduced action.

## 3. Constructing the collar rather than imposing zero normal reaction

For the connected torus surface, any smooth desired normal displacement
`b` with `integral_Sigma b=0` has a unique mean-zero solution
`Delta_Sigma psi=b`. In a tubular collar, extend its tangential potential
`n cross grad_Sigma psi` smoothly along signed normals and multiply by a
smooth cutoff equal to one near the surface and zero at the collar edges.
Then

    E[b]=curl(cutoff(s) n cross grad_Sigma psi)

is smooth, divergence free, supported in the collar, and has normal trace
`b`. The last assertion follows from the intrinsic Stokes formula:
`n.curl(n cross grad_Sigma psi)=Delta_Sigma psi`. The collar can be chosen
disjoint from all physical core-angle observations. Its extension into the
ambient is part of the same field, not an independently slipping interface.
The zero-flux condition is volume preservation, automatically satisfied by
the trace of a global divergence-free velocity/displacement.

This explicitly completes a chosen core material lift `Xi_core` to
`Xi=Xi_core+E[b]`. Derive its velocity rather than retaining an incompatible
old one:

    v=Xi zdot+[u,Xi]z.                                       (4)

The surface normal displacement is `b=(Xi.n)z`. Its exact equation is
`b_t+L_Sigma b=v.n`, where `L_Sigma b=div_Sigma(b u_tangent)`. The density
divergence is important: invariance of one torus need not supply a stationary
first integral throughout a collar. For signed distance `F`, one has
`u.grad F=aF+O(F²)` with `a=-div_Sigma u_tangent`; hence
`[u,Xi].n=L_Sigma(Xi.n)` on the surface. Equation (4) therefore proves the
tag transport equation exactly. If an existing velocity sector is
`v=V1 zdot+V0 z`, its tag is local with this zero-order-in-time material form
precisely when its normal coefficients satisfy

    V0.n=L_Sigma(V1.n).                                     (5)

with compatible prepared initial tag. Choose `b=V1.n` to construct its
boundary lift when (5) holds. This is the same compatibility derived by
0082's transport representation. If it does not hold, the explicit repair
is to retain the transported boundary state, or replace the velocity family
by (4) and recompute the action. The added commutator is not negligible by
definition, and the old orbit Hessian is not preserved by that replacement.

## 4. Recompute the whole Jacobi action and preserve Kelvin data explicitly

For a global smooth lift with its ambient extension, use 0037's full action

    L2=1/2 integral [rho|Dt xi|²-xi.Hess(p0).xi].

The pressure integration-by-parts surface term cancels between the tube and
its ambient. A separated tube calculation instead retains that term as the
pressure reaction. Substitution of the actual lift gives

    L2=zdot^T M zdot/2+zdot^T G z-z^T K z/2,
    M_ab=rho integral Xi_a.Xi_b,
    G_ab=rho integral Xi_a.(u.grad Xi_b),
    K_ab=integral Xi_a.Hess(p0).Xi_b
               -rho integral (u.grad Xi_a).(u.grad Xi_b).      (6)

The integrals include the complete chosen ambient lift. `M` is a genuine
material Gram matrix; `G` and both terms in `K` remain. Eliminating further
shape coordinates requires their full operator, not an isolated tube inverse.

The selected material variations do not automatically preserve a fixed
Kelvin leaf. Their first variation of the pulled-back velocity one-form is
represented in spatial coordinates by

    C_K=xi_t+(u.grad)xi+(Dxi)^T u.                            (7)

Fixed Kelvin momentum requires this one-form to be exact, including zero
periods. Equivalently it gives `v=P(xi cross omega)` on the global
simply-connected decaying/compatible domain, with the appropriate harmonic
conditions retained on other domains. For `curl u=lambda u`, (7) becomes

    xi_t=(lambda P-curl)(xi cross u),                         (8)

because `[u,xi]=curl(xi cross u)`. This checks the reconstruction convention
directly. A material-first restriction preserves the full material transport
and volume constraint; it preserves fixed Kelvin data only after (7) is
imposed or the actual relabeling momentum reduction has been performed.
Equation (8) makes clear why simply declaring a finite generator span
invariant is a new substantive assertion. No such assertion is used here.

For a genuine linearized Euler solution and compatible initial Kelvin data,
the reconstruction in section 1 preserves (7), by the linearized Kelvin
theorem. For a conditionally restricted action, the retained circulation
constraints and reaction terms must be stated; a kinematic material metric
alone does not supply a fixed-leaf positive spin inertia.

## 5. Exposing exact Jacobi sign: a smooth Beltrami core/collar example

The general static identity can first be derived without choosing a profile.
Set `a=(u.grad)Xi`, `b=(Xi.grad)u`, `w=a-b=curl(Xi cross u)`. Stationarity
and one integration by parts along `u` give

    integral Xi.Hess(p0).Xi/rho
        = integral Xi.(w.grad u)+a.b.

Subtract `integral |a|²` and use
`(Du)^T Xi-(Du)Xi=Xi cross omega`. Therefore

    K/rho=integral curl(Xi cross u).(Xi cross omega)
                       -|curl(Xi cross u)|².

For `omega=lambda u` and `F=Xi cross u`, this is precisely

    K/rho=lambda integral F.curl F-integral |curl F|².        (9a)

The integration-by-parts flux is
`-rho div[u (Xi.(Xi.grad u))]`; it vanishes for the invariant-domain
background or the full compatible global action. The separate pressure
surface term in the material action is still retained as specified above.
Equation (9a) supplies a genuine positive-candidate continuation: the
same-helicity curl band `0<k<lambda` is positive, unlike the old opposite-
helicity high-carrier orbit cage. Realizing this band with `F=Xi cross u`
and the actual core/ambient constraints is a substantive construction, not
inferred from the spectral sign alone. Main attempt 0088 owns that search.

Take the smooth stationary Beltrami field

    u=(sin(lambda z),cos(lambda z),0), lambda>0,
    curl u=lambda u, p0=constant.

Every compact smooth material lift has

    K(Xi,Xi)=-rho integral |u.grad Xi|² <= 0.                 (9)

For an explicit local core/collar curl, choose
`Xi=curl(f(x)g(y)h(z)e_z)`, where `f,g,h` are nonzero smooth compact bumps
and `g` has a nonzero plateau. On that plateau the second velocity component
has directional derivative `-sin(lambda z) f''(x)g h(z)`.
Choose `h` supported where `sin(lambda z)!=0`; a nonzero compact `f` cannot
have `f''` identically zero. Thus (9) is strictly negative. These supports
may be placed in a finite flow-box neighborhood; the example uses an actual
smooth Beltrami background and full Jacobi form, not a singular core model.
It is a sign-exposing local example, not a knotted invariant torus.

For one compact lift, `G=rho integral u.grad(|Xi|²/2)=0`. Consequently its
unrestricted one-coordinate material Jacobi projection has negative static
stiffness. This does not refute fixed-Kelvin orbit positivity: the arbitrary
one-coordinate path has not obeyed (7). It does show why a high-carrier
material replacement cannot borrow the positive orbit-Hessian cage sign.

The constructive continuation is therefore explicit: use sections 1–3 to
transport or lift the actual tag, and impose/reduce (7) in the complete
same-action shape sector before using (3) as the spin observation map. The
remaining EPS-sector calculation is that constrained operator reduction;
neither a new phenomenological inertia nor the sign of `M` can replace it.

## Result and effect on the continuing parent construction

`route_verdict: established` for the exact material reconstruction, collar
extension, full local spin map and compatibility/action formulas.
`evidence_scope: REPRESENTATION_SCOPED` — a constructive material transport
license and an exact sign-exposing Beltrami example. This produces an actual
repair representation, not only a rank identity or an obstruction label.
It does not establish that the old 0080 instantaneous sector already obeys
(5), (7), or has positive reduced `K`. The parent must apply the computed
lift or retained tag state and its full Kelvin/reaction reduction to its
chosen sector; that is the next executable construction, not exhaustion.
