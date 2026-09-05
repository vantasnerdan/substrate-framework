# Positive compact angle/reaction pair on the actual EPS field

This continues `compact-spin-theorem.md`. It recomputes the full
coadjoint Hessian and KKS on the compact-velocity sector; it does not
import a material Jacobi coefficient or the old circular-cage sign.

## 1. Exact identities and the necessary subprincipal return

Let lambda!=0, omega=lambda u, and let S be any nonzero real analytic
finite-order vector differential operator constructed above with

    div(S psi)=div[(S psi) cross omega]=0

for every compact smooth psi. On its support, write xi=S psi and
v=xi cross omega. Exactly,

    H(xi,eta)=rho integral [v_xi dot v_eta-v_xi dot curl v_eta/lambda],
    Omega(xi,eta)=rho integral omega dot(xi cross eta).

The second term of H is the complete isovortical second variation. Its
local appearance here does not assert that the second-order induced
velocity of an arbitrary finite-amplitude orbit is compact. Integration
by parts also gives

    H(xi,xi)=-2rho integral v dot Sym(grad xi) u.

No pressure-gradient contribution has been discarded: div v=0 makes
its integral vanish, and the complete orbit Hessian remains as written.

Fix a constant real unit vector n, a compact smooth envelope, and a
phase convention cancelling the leading i^N of the order-N operator.
After normalization by k^-N, its cosine/sine outputs have expansions

    xi_c=a cos(k n.x)+k^-1 b sin(k n.x)+O(k^-2),
    xi_s=a sin(k n.x)-k^-1 b cos(k n.x)+O(k^-2).

The coefficients include the full derivatives of the envelope and
analytic operator coefficients. They are not independent trial
polarizations. The two exact divergence identities force

    n.a=0, n.(a cross omega)=0,
    n.b=-div a,
    (omega cross n).b=-div(a cross omega).

For generic n with n cross omega!=0 and n.omega!=0 this gives

    a=A P_n omega, f=a cross omega=c(n cross omega), c=-A(n.omega),
    g=b cross omega,
    g.(n cross f)=-c(n.omega) div f.

Thus the leading helicity term is NOT just f.curl f. Direct expansion
of the real oscillations gives

    average(v_c.curl v_c)
       =1/2 [f.curl f-2g.(n cross f)]+O(k^-1).

The return g contributes at the same order as f.curl f and is retained.

## 2. Corrected energy symbol and nonzero symplectic pair

Put t=n cross omega. Since div t=-lambda(n.omega), integration by parts
of the compact envelope coefficient c gives

    H(xi_c,xi_c)=H(xi_s,xi_s)
       =rho/2 integral c^2 h_n+O(k^-1),
    H(xi_c,xi_s)=O(k^-1),

    h_n=2|n cross omega|^2+(n.omega)^2
         +2(n cross omega).(n.grad)omega/lambda.

For clarity, the intermediate helicity expression is

    f.curl f-2g.(n cross f)
      =c^2[t.curl t-2lambda(n.omega)^2]
          +(n.omega)t.grad(c^2).

Integrating its last term and using
grad(n.omega)-(n.grad)omega=lambda(n cross omega) gives h_n above.
For an orthonormal triad n_i,

    sum_i (n_i cross omega).(n_i.grad)omega
        =-omega.curl omega=-lambda|omega|^2,
    sum_i h_(n_i)=3|omega|^2>0.

Therefore some direction has strictly positive h_n. Positivity is open
in n; choose it also away from the parallel/perpendicular exceptional
directions and from the zero set of the nonzero principal symbol of S.
On a sufficiently small support ball, h_n has a positive lower bound,
n.omega has one nonzero sign, and the coefficient c is not identically
zero. This choice depends only on actual local geometry and the
operator, not on a desired numerical modulus.

The SAME pair gives, including the subprincipal return,

    Omega(xi_c,xi_s)
       =-rho lambda/(2k) integral c^2(n.omega)+O(k^-2).

Indeed its leading cross product is -k^-1 a cross b. The constraints
above give omega.(a cross b)=-c div f, and compact integration yields
integral c div f=-lambda/2 integral c^2(n.omega). Hence the leading
coefficient B0 is nonzero on the chosen ball.

All remainders have finite explicit derivative bounds. S has finite
order, so expanding its outputs leaves finitely many powers of k^-1;
oscillatory integrals with phase2k n.x are bounded by one or more exact
integrations by parts in n. Record the resulting finite constants
CH,CB and h0=rho integral c^2 h_n/2>0. Then for a finite sufficiently
large k,

    H_pair >= (h0/2) I,
    |B_c-B0/k|<=CB/k^2,  B_c!=0,
    |H(xi_c,xi_s)|<=CH/k.

Unlike an arbitrary numerical truncation, increasing k here does not
approach a floating-point sign threshold: the strict leading bounds
and the finite analytic remainders specify the finite choice.

## 3. Independent exact exposing oracle

`wkb_pair_verify.py` uses the fully explicit plane Beltrami field
omega=(cos z,sin z,0), lambda=-1, and a differential syzygy whose complete
xi and xi cross omega are divergence-free exactly. For the envelope
1+delta cos z it verifies the integrated leading coefficients

    H=(35delta^2+52)/2048>0,
    Omega=7delta/256.

The degrees in its unnormalized order-five carrier are k^10 and k^9.
These are exact trigonometric Laurent constant terms, not quadrature.
Its nine checks include both divergence identities, the complete
quadrature contribution, the KKS sign, and the universal triad identity.
The example is an exposing algebraic oracle, not the actual EPS
background: the general analytic argument in section2 supplies that
transfer directly to its local field and compact operator.

## 4. Attach an actual core angle and match its physical spin exactly

Work inside ONE actual invariant EPS solid torus with the open
finite-jet rank property. Use mutually disjoint support balls for:
one core observation, one positive compact cage pair Q_c,S_c, and
three compact-spin responses eta_i with mechanical spins e_i. The last
fields are the normalized determinant/bump construction of the compact
spin theorem. All are contained strictly inside D.

A nonzero principal operator symbol gives a genuine core-angle
direction. For generic n its induced vorticity principal part is

    curl v proportional to n cross(n cross omega),

which has a nonzero component perpendicular to omega when both
n.omega and n cross omega are nonzero. A compact oscillatory bump equal
to one near the observation point therefore gives a nonzero transverse
vorticity-jet variation at a finite carrier. Define the physical local
angle from the direction of the actual vorticity jet, and normalize
that direction Q0 to have angle derivative one. This is an Eulerian,
relabeling-invariant core observation, not a label-template frame.
The off-core cage and spin-response supports have zero derivative of
that observation. No desired energy coefficient sets the normalization.

Let L denote the actual three-component mechanical spin functional.
Choose a fixed positive cage amplitude A and set

    Q=Q0+A Q_c-sum_i eta_i L_i(Q0+A Q_c),
    S=S_c+sum_i eta_i [A B_c e_i^(axis)-L_i(S_c)].

Here e^(axis) is the physical angular-momentum axis of the selected
core angle. These formulas give exactly

    core-angle(Q)=1, core-angle(S)=0,
    L(Q)=0, L(S)=A B_c e^(axis),
    Omega(Q,S)=A B_c=:B!=0.

Disjoint supports imply every raw-to-response and response-to-response
KKS cross is zero. Since the induced velocities too are compact in
these supports, the COMPLETE H cross terms between different balls
are zero. This is an exact locality result, not an independent-cell
inverse approximation.

The normalized cage spins satisfy |L(Q_c)|+|L(S_c)|<=CL/k^2 (and in
fact decay faster than any prescribed power) by compact oscillatory
integration of the angular adjoint row. Let E_eta bound the finite
spin-response energy matrix and let l0=L(Q0). The fixed part
Q0-eta l0 has finite energy, of either sign. First choose A so that
A^2 h0 dominates its negative part with a strict margin. Then choose
the finite k large enough that:

    H_QQ >= A^2 h0/4,
    H_SS >= h0/2,
    |H_QS| <= Ccross(A)/k,
    B!=0.

Here the reaction spin correction is O(A/k); its energy is O(A^2/k^2).
The Q spin correction beyond the fixed l0 is O(A/k^2).
For example Ccross(A) follows directly from CH,CL,E_eta,l0 and
|B0|+CB by the triangle inequality. Taking also
Ccross(A)^2/k^2 < A^2 h0^2/8 makes the FULL2 by2 Hessian positive
definite. This hierarchy uses a finite fixed core/response geometry,
then a fixed A, then a finite analytic carrier choice. No sign from an
unmeasured solver or borrowed cage energy is used.

## 5. The actual quadratic same-Euler rotor and physical current

Both Q and S now satisfy div xi=div(xi cross omega)=0, have compact
induced velocity inside D, and have zero actual normal boundary trace.
The reference tag transport has no omitted Q/S normal source. Their
centroid and STF velocity moments vanish automatically. Their actual
tube spin equals their global spin; the ambient induced velocity is
zero. Shared-interface pressure work remains the Euler one, with no
new wall, altered Leray projector or assigned rotor mass.

Write H_SS=P>0, H_QS=N and H_QQ=Hq. The selected SAME Euler action is

    L=B s Phidot-[P s^2+2N s q+Hq q^2]/2,
    q=Phi-beta, beta=curl U_centres/2,

once the existing affine material lift Gamma_K=K-Q is used. The exact
spin and affine moment identities give Omega(Gamma_K,S)=0 on this
compact pair, and common rigid rotation Gamma_K+Q=K preserves the
full stationary Euler energy. The parent assembly supplies its shared
macro action and retained gradient terms; no independent global rotor
is added here.

Reaction elimination gives

    j=B^2/P>0, K=Hq-N^2/P>0,
    L_reduced=j Phidot^2/2-K q^2/2-(BN/P)q Phidot.

On a single background the actual spin is
L_physical=B s=j Phidot-(BN/P)q. Thus N is NOT silently removed by
naming a new momentum. Under the prescribed equal time-reversal
background pair u->-u, H is unchanged while B and physical spin rows
change sign. After full reaction elimination, the mixed term and its
static angle-spin current cancel between the two backgrounds; j and K
stay positive. The paired actual spin is j Phidot. Reflection pairing
and spatial-gradient terms are likewise applied to the full action,
not to a prematurely averaged isolated inverse.

This proves a positive compact physical angle/reaction sector and its
boundary-compatible first-order material tag/current data on the
actual stationary EPS field. It does not declare the parent continuum
gradient/translation/registry obligations complete, nor assert an
unrestricted finite-amplitude invariant manifold. The complete orbit
second variation and all momenta have been retained at the stated
quadratic-action scope.
