# Joint charged-Cao Hessian and radiation-channel derivation

## 1. Instantaneous tangent and the exact comoving field block

Write `epsilon=epsilon_EM`, `mu=mu_EM`,

    c_EM^2=1/(epsilon mu),   beta^2=epsilon mu c_g^2<1.       (3)

At a fixed time the Maxwell tangent is constrained by

    div delta B=0,   epsilon div delta E=g delta chi.         (4)

Ampere and the material current are evolution equations, not extra
instantaneous constraints.  Gauge fixing is needed when potentials are used;
the field calculation below uses only the gauge-invariant `delta E,delta B`.

With field momentum density `epsilon E cross B`, direct expansion gives

    Q_EM,c[E,B]
      = integral {epsilon |E|^2/2+|B|^2/(2mu)
                  -c_g epsilon (E cross B)_z} dx
      = (epsilon/2)||E+c_g e_z cross B||_2^2
        +(1/(2mu))(1-beta^2)||B_perp||_2^2
        +(1/(2mu))||B_z||_2^2.                              (5)

There is no field-equation substitution in (5).  In scaled perpendicular
variables `e=sqrt(epsilon)E_perp`, `b=B_perp/sqrt(mu)`, the two polarization
blocks have eigenvalues `1-beta` and `1+beta`.  Consequently

    Q_EM,c >= (1-|beta|)/2
                  [epsilon||E_perp||_2^2+mu^-1||B_perp||_2^2]
               +epsilon||E_z||_2^2/2+mu^-1||B_z||_2^2/2.    (6)

Thus the field block is positive and norm-equivalent to the finite-energy
field norm for every fixed subluminal margin.  Its lower constant degenerates
as `|c_g|` approaches `c_EM`.

## 2. Exact Gauss-fiber Schur complement

Set `q=g delta chi` and use a unitary Fourier convention, so harmless common
Plancherel factors do not appear.  For each nonzero `k`, impose

    i k dot Ehat=qhat/epsilon,   k dot Bhat=0.               (7)

Let `t=k_z/|k|`, `s=|k_perp|/|k|`, and choose the orthonormal frame
`(n,e_1,e_2)` with `n=k/|k|`, `e_1` in the `(k,e_z)` plane, and
`n cross e_1=e_2`.  The Gauss row fixes
`|E_n|=|qhat|/(epsilon |k|)`.  The only sourced transverse pair in (5) is
`(E_1,B_2)`.  Before minimization its real quadratic is

    epsilon(E_n^2+E_1^2)/2+B_2^2/(2mu)
      -c_g epsilon(t E_1+s E_n)B_2.                        (8)

The stationary values are

    E_1=c_g t B_2,
    B_2=c_g epsilon mu s E_n/(1-beta^2 t^2).               (9)

The other transverse pair minimizes at zero.  Substitution yields

    Q_min(k)
      = |qhat(k)|^2/(2epsilon |k|^2)
          (1-beta^2)/(1-beta^2 k_z^2/|k|^2)
      = (1-beta^2)/(2epsilon)
          |qhat(k)|^2/(|k_perp|^2+(1-beta^2)k_z^2).        (10)

For fixed total charge, `qhat(0)=0`; the affine Coulomb monopole of the base
does not enter the tangent integral.  Since

    (1-beta^2)|k|^2
      <= |k_perp|^2+(1-beta^2)k_z^2 <= |k|^2,              (11)

the exact integrated Schur form obeys

    (1-beta^2)/(2epsilon)||q||_{Hdot^-1}^2
      <= Q_min[q] <= 1/(2epsilon)||q||_{Hdot^-1}^2.        (12)

It is positive, `O(g^2)` on `delta chi`, and controls precisely an anisotropic
homogeneous `H^-1` charge norm.  It supplies no `L2` or positive-order
Sobolev tag estimate.  Equations (7)--(12), rather than Ampere constraints or
a guessed Coulomb sign, are the field Schur complement needed by the
instantaneous Hessian.

## 3. The full joint relative Hessian is a strong saddle

The exact Cao steady relation on the positive core is

    Psi_0=G_1 zeta_0-c_0 r^2/2-mu_0,
    zeta_0=f(Psi_0)=epsilon_core^-2(Psi_0)_+^p.             (13)

In Cao's meridional measure `dnu=r dr dz`, normalize the rearrangement
functional as

    E_c[zeta]=1/2 integral zeta G_1 zeta dnu
               -c_0/2 integral r^2 zeta dnu.              (14)

On a regular positive-core band introduce the local Casimir with
`C'(zeta_0)=Psi_0` (the constant chemical-potential row is absorbed by fixed
circulation).  Then

    C''(zeta_0)=1/f'(Psi_0)
      =epsilon_core^2/[p Psi_0^(p-1)]>0,                  (15)

and, up to the common positive physical normalization, the constrained fluid
second variation is exactly

    K_fl[q,q]=integral q G_1 q dnu
               -integral q^2/f'(Psi_0) dnu.               (16)

This sign also matches the source theorem: the uncharged Cao ring maximizes
`E-c_0 P` on its axisymmetric rearrangement class.  It does not imply that the
joined fluid--Maxwell functional is a maximizer.

Choose a compact regular meridional cell outside the compact tag support but
inside the positive core, with both `f'(Psi_0)>0` and `zeta_0'(I)!=0`.  In its
enclosed-area action-angle coordinates, choose `a(I)` compactly supported in
that cell and the area-preserving Hamiltonian generator

    S_N(I,beta)=-a(I) sin(N beta)/[N zeta_0'(I)].          (17)

Then `xi_N^I=partial_beta S_N`,
`xi_N^beta=-partial_I S_N`, so the displacement preserves `dI d beta` and
its dynamically accessible vorticity is exactly

    q_N=-xi_N dot grad zeta_0=a(I) cos(N beta).            (17a)

The packets have zero angle mean.  The finitely many circulation, center, and translation rows are
removed by graph-continuous corrections of lower order.  On the fixed cell,
`1/f'(Psi_0)` has a positive lower bound, while the compact order-minus-two
Green term satisfies

    integral q_N G_1 q_N dnu=O(N^-2)||q_N||_2^2.          (18)

Therefore `K_fl[q_N,q_N]<0` for every sufficiently large `N`, and disjoint
frequency packets expose the operator mechanism.  More precisely, on the
infinite-dimensional DA zero-angle-mean subspace of this cell,

    K_fl=-M_(1/f'(Psi_0))+G_1,                            (18a)

where `M_(1/f')>=m I` for some `m>0` and the compressed Green operator is
compact self-adjoint.  Compactness gives a finite-codimension subspace on
which the Green norm is less than `m/2`; (18a) is at most `-m/2` there.
Therefore the negative index is infinite, rather than merely containing
infinitely many individually negative vectors whose spans were uncontrolled.
Because the cell avoids the tag support, `delta chi=0`, so the Gauss row
permits `delta E=delta B=0`.

Independently, take any smooth compact axisymmetric toroidal electric field
`delta E=e_theta e(r,z)`, with `delta B=0` and zero fluid/tag variation.  It
is divergence-free, satisfies (4), and (5) gives

    Q_EM,c[delta E,0]=epsilon||delta E||_2^2/2>0.          (19)

There are infinitely many orthogonal such radiation directions.  Hence the
uncharged joint relative Hessian has infinite positive and negative index,
already in the axisymmetric joint phase space.  Multiplying the entire
relative functional by `-1` merely exchanges the signs.

Using the independently reviewed P253/0080/0084 charged branch, its fluid,
speed, and tag change is `O(g^2)` and its base field is
`O(g)`.  On one smaller common cell separated from the tag support, the
charged coefficients converge uniformly in `C^1`, `1/f'_g` keeps one uniform
positive lower bound, and the action-chart compressions of the Green operator
form a uniformly compact family.  After one common finite-dimensional sector
is removed, the charged fluid form is uniformly negative on its
infinite-dimensional complement.  Equivalently, there are common numbers
`g_0>0` and `N_0` and a sparse almost-orthogonal packet basis whose entire
closed tail span, for `|g|<g_0`, is strict negative.
The total momentum, center, circulation, and phase conditions add only
finitely many continuous rows to this cell construction; take their common
kernel or remove them with uniformly bounded finite-rank corrections.  At a
charged base, a pure toroidal `delta E` can likewise have an `O(g)` linear
total-momentum row through the nonzero base field.  The infinite-dimensional
kernel of those finitely many field rows retains the positive lower bound
(6), or equivalently uniform finite-rank projection restores the slice.
Thus the pure field directions in (19) remain infinite-dimensional and strict
positive on the actual constrained tangent.  Hence small charge
cannot turn the full joint relative energy into a one-sign coercive form or a
finite-Morse-index form.  This is a Hessian statement, not a proof of spectral
instability: Hamiltonian systems can be spectrally stable with an indefinite
energy.

## 4. What the material lock removes, and what it does not

Use regular action coordinates `(I,beta,theta)` with physical vorticity

    Omega=zeta(I) partial_theta,   chi=F(I).               (20)

For `delta Omega=-[xi,Omega]` and a toroidal Fourier harmonic `n`,

    delta Omega^I=i n zeta xi^I,
    delta Omega^theta=-zeta' xi^I+i n zeta xi^theta,
    delta chi=-F' xi^I.                                   (21)

On the tagged band `zeta zeta'!=0`, this gives

    n=0:   delta chi=(F'/zeta') delta Omega^theta,
    n!=0:  delta chi=-(F'/(i n zeta)) delta Omega^I.       (22)

In particular a same-orbit stabilizer `delta Omega=0` also has
`delta chi=0` there.  The carefully supported tag therefore removes the
otherwise real defect of changing charge at fixed vorticity within the
selected joint orbit.  Different tag profiles over the same fluid state are
still distinct semidirect-product leaves, and the Coulomb form (12) still
degenerates like `g^2` and controls only `H^-1`.  Equations (21)--(22) do not
create a physical tag restoring energy.

## 5. Maxwell radiation embeds every imaginary internal frequency

Let the Maxwell tangent Hilbert space be the finite-energy source-free fiber

    H_EM,0={(delta E,delta B) in L2 cross L2:
             div delta E=div delta B=0},

with its epsilon/mu energy norm, and let

    D_EM,0=H_EM,0 intersect [H(curl) cross H(curl)].      (22a)

The full joint Hilbert space appends the declared fluid and transported-tag
spaces and imposes the finite circulation, center, impulse, total-momentum,
phase, and gauge rows.  Its closed generator `A_joint,g` has the corresponding
product graph domain, with the finite rows taken as a closed finite-codimension
subspace.  Define `sigma_ess` as the Fredholm essential spectrum: `lambda`
belongs to it when `A_joint,g-lambda` is not Fredholm.  The construction below
in fact shows failure of upper semi-Fredholmness.

At `g=0` the full generator is a direct sum of the fluid, passive-tag, and
free Maxwell generators.  In the frame translating with speed `c_g`, a free
Maxwell plane wave has

    lambda_+(k)=-i c_g k_z+i c_EM|k|,
    lambda_-(k)=-i c_g k_z-i c_EM|k|.                     (23)

For `|c_g|<c_EM`, each branch ranges over an imaginary half-line and their
union is `i R`.  Here is the joint graph-domain construction.  Fix
`lambda=i omega` with `omega!=0`, choose a nonzero `k` and transverse free
Maxwell polarization `U_0=(E_0,B_0)` satisfying (23), and choose constant
potentials `a_E,a_B` with

    i k cross a_E=E_0,
    i k cross a_B=B_0.                                    (23a)

For a fixed smooth compact cutoff `chi`, put

    chi_L(x)=L^(-3/2)chi((x-x_L)/L),

    E_L=curl[chi_L a_E exp(i k dot x)],
    B_L=curl[chi_L a_B exp(i k dot x)].                   (23b)

Both fields are smooth, compact, exactly divergence free, Gauss free, and in
`D_EM,0`.  After energy normalization their main terms are
`chi_L(E_0,B_0)`; every cutoff correction and free-Maxwell residual is
`O(L^-1)`.  Choose the centers so

    dist(x_L,supp carrier)/L -> infinity                 (23c)

and, for the actual `r^-2` base Maxwell tail, fast enough that its coefficient
pairings with the packet and one derivative are `o(1)`.  The fluid vorticity,
tag, and source current are compact, so their direct couplings vanish once the
supports separate.  Translation to infinity gives weak convergence to zero.

For `omega=0`, subluminality leaves no nonzero fixed-`k` root.  Choose instead
`|k_L|=L^(-1/2)` on either Maxwell branch and potentials satisfying (23a).
The symbol eigenvalue and cutoff-potential correction are then both
`O(L^(-1/2))`, so the normalized sequence again has vanishing residual and is
weak null.

The finitely many joint momentum, gauge, center, and constraint functionals
are graph continuous and tend to zero on these sequences.  Subtract their
values against fixed smooth source-free profiles dual to the independent
rows.  The correction is `o(1)` in graph norm and preserves the residual.
Thus for every `lambda in i R` there is a normalized, constrained, weak-null
sequence `U_L in D(A_joint,g)` with

    ||(A_joint,g-lambda)U_L|| ->0.                        (23d)

If `A_joint,g-lambda` admitted a left Fredholm regularizer
`S(A_joint,g-lambda)=I-K` with `K` compact, (23d), weak nullness, and
compactness would force `U_L->0` in energy norm, contradicting normalization.
Therefore

    i R subset sigma_ess(A_joint,g)                       (24)

for the localized charged carrier, independently of whether equality is
proved.  The passive-tag transport block can add essential spectrum and
cannot restore a gap.

Consequently the positive-Krein Cao modes reviewed in P253/0074/0078 are
embedded in the full joint continuum.  Ordinary Riesz/Kato persistence is
available for a genuinely isolated off-axis eigenvalue `Re lambda!=0` under
the usual common-domain relative-bounded perturbation hypotheses.  It is not
available for an imaginary internal mode solely because that mode was
isolated in a fluid fiber.

For a comoving harmonic perturbation
`exp(i k dot(x-c_g t e_z)-i omega t)`, the Maxwell wave denominator is

    D_omega(k)=|k|^2-(omega+c_g k_z)^2/c_EM^2.             (25)

At `omega=0`, subluminality makes (25) elliptic away from `k=0`, explaining
the steady affine Green solve.  For every `omega!=0`, its zero set is a
nonempty smooth real radiation shell in every angular direction: writing
`k=r n`, its positive radial root is

    r(n)=|omega|/[c_EM-sign(omega)c_g n_z]>0.             (25a)

The radial-root denominator in (25a) stays uniformly positive because
`|c_g|<c_EM`; `D_omega` itself vanishes on the shell.  The
linearized material current is

    delta J=g(chi_g delta u+u_g delta chi).                (26)

The shell is star-shaped.  Differentiating (25) radially at (25a) gives

    |partial_r D_omega(r(n)n)|=2|omega|/c_EM,             (26a)

so its exact coarea weight is

    r(n)^2/|partial_r D_omega|
      = |omega| c_EM/
          [2(c_EM-sign(omega)c_g n_z)^2].                (26b)

Accordingly the nonnegative shell functional, up to the frozen Fourier and
Maxwell constants, is

    Gamma_omega[delta J]
      = integral delta(D_omega(k))
          |P_T(k) delta Jhat(k)|^2 dk
      = |omega|c_EM/2 integral_(S2)
          |P_T delta Jhat(r(n)n)|^2/
          (c_EM-sign(omega)c_g n_z)^2 dOmega.            (26c)

It vanishes exactly when the transverse trace vanishes almost everywhere on
the shell, provided that trace exists.

A localized finite-energy time-periodic eigenfield necessarily needs the
transverse on-shell cancellation, with enough trace regularity and vanishing
order to cancel the simple resolvent pole,

    P_T(k) delta Jhat(k)=0
      whenever D_omega(k)=0,                              (27)

together with the longitudinal continuity/Gauss row on every shell component.
This is necessary, not sufficient: sufficiency still requires the weighted
outgoing Maxwell resolvent and coupled fluid feedback.  If (27) is nonzero,
the prescribed real-frequency source has no localized `L^2` Maxwell field.
This activates a weighted limiting-absorption plus analytic-Fredholm/Feshbach
construction of a self-consistent outgoing resonance; neither existence nor
the width/sign is inferred from the trace alone.  Compact support makes `delta Jhat` analytic, but vanishing on a
codimension-one shell does not force the current to vanish.  If (27) vanishes
identically, the mode is a genuine nonradiating/BIC candidate, not yet a
proved eigenmode.

SO(2) character alone does not force (27): free Maxwell waves contain every
integer azimuthal character.  A useful P253/0079 analyzer mode therefore
needs a source-specific dark-current identity or a quantified radiative
leakage bound over its full gate time.

## 6. Route verdicts and continuation

**Route A — full one-sign joint coercivity:** refuted for the declared
whole-space joint phase space.  The exact pure-field positive directions and
high-frequency Cao-orbit negative directions give infinite positive and
negative index.  The positive Gauss-fiber Schur form and tag-locking identities
are retained as useful exact subclaims.

**Route B — ordinary spectral persistence of neutral Cao modes:** refuted.
Whole-space Maxwell radiation puts `i R` in the joint essential spectrum.
Direct persistence remains a valid route for an off-axis unstable mode; no
such mode is asserted here.  The exact shell and dark-current condition
(25)--(27) replace the invalid neutral Riesz route.

**Route C — reduced or radiative persistence:** active with two concrete
branches.  First, restrict only the *steady* Maxwell variables by the exact
subluminal Green solve and test whether the axisymmetric rearrangement
maximizer persists for the reduced functional.  On a compact regular tagged
cell, (22) and the order-minus-two steady Maxwell inverse make its `O(g^2)`
charge/current Hessian compact relative to the negative multiplication term
in (18a).  Thus the high-frequency reduced form remains uniformly negative
and only a finite-dimensional low sector can change sign.  The exact sign of
that low block is still open, and even a favorable answer cannot be promoted
to full-dynamics stability because independent Maxwell radiation remains.
Second, compute (27) and its limiting-absorption dual for the two actual
0079 crossing modes.  A zero gives a symmetry/source-specific dark candidate;
a nonzero value gives a quantitative `O(g^2)` radiation-resonance channel and
the gate-time leakage test.  Modulated local decay with the radiation field
retained is the appropriate positive theorem after the latter outcome.

The bare-Euler Cao and weighted-tail routes remain active.  None of these
results selects a universal action or charge, constructs measurement or
exchange statistics, supplies a full Lorentz cone, or identifies an electron
or neutrino.
