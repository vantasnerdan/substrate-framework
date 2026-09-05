# Hybrid tube-centroid / continuous-ambient Euler action

## 1. The partition is material and the ambient remains continuous

Partition reference labels into bounded disjoint tube parcels D_a0 and their
ambient complement A0. The maps g_a and g_A are restrictions of ONE smooth
volume-preserving map g. Their shared traces agree; the pressure multipliers
and their shared-face work are retained. At the stationary reference Euler
state each selected tube D_a is invariant, u0.n=0, and
`integral_Da u0=0`, so its base mass-centroid velocity is zero.

For each tube define mass m_a, actual mass centroid X_a, r_a=g_a-X_a and
centered momentum pi_a. The exact full kinetic/cotangent decomposition is

    T=sum_a [|P_a|²/(2m_a)+integral_Da0 |pi_a|²/(2rho)]
          +integral_A0 |p_A|²/(2rho),
    Theta=sum_a [P_a.dX_a+integral_Da0 pi_a.dr_a]
          +integral_A0 p_A.dg_A.

Here `integral_Da0 r_a rho=integral_Da0 pi_a=0`. The ambient term is a
continuous material integral, not a sum of independently resetting cells.
All ambient kinetic energy and tube/ambient relative motion are present.
Normalizing a stationary ensemble gives total mass density rho because the
tube masses plus the continuous ambient count every particle exactly once.

## 2. The physical hybrid momentum and its exact observable map

As a spatial distribution define

    p_H(x)=sum_a P_a delta(x-X_a)+1_A(x) rho u(x).

It differs from whole-fluid point momentum only inside the tubes. Let
I_a,ij=integral_Da rho r_i r_j,
S_a=integral_Da rho r cross (u-V_a), V_a=Xdot_a. Then

    Q_a,ij=integral_Da rho (u_i-V_a,i) r_j
           =Idot_a,ij/2-epsilon_ijm S_a,m/2.

The exact first-multipole relation is therefore

    p_E=p_H-div(Idot_tubes)/2+curl(S_tubes)/2+R_2.

Only TUBE moments enter this correction. The ambient point momentum is the
same on both sides. R_2 is the distributional second-and-higher Taylor
remainder, not a discarded term in a claimed second-gradient coefficient.
For a test field f its leading coefficient includes
`integral_D rho u_i r_j r_l`, including both translation and internal parts;
the integral Taylor remainder controls the subsequent terms.

This physical observable map is independent of a constitutive spin law.
If a complete material reduction establishes S_tubes=j_tube Phidot and the
shape response, its leading spin part gives
`rho U_Edot=rho U_Hdot+curl(j_tube Phidot)/2`. It does not identify U_H with
the point mean corrected in 0072. It also does not permit dropping Idot,
the second moment or phase-relative velocities when computing the full jet.

## 3. Actual tube spin versus global angular impulse

At a reference tangent whose generator xi vanishes near each tube boundary,
the actual tube spin response is

    A_D(xi)=rho integral_D r cross v_xi,
    v_xi=P_R3(xi cross omega).

This is NOT generally the free-space angular impulse of xi. Vorticity
rearrangements outside the tube can induce velocity inside it. For a general
moving-interface tangent, add the Reynolds term
`rho integral_boundary r cross u0 (xi.n)`; the centroid correction vanishes
at the invariant base because integral_D u0=0. All such interface terms
remain in the full material chart below.

Let Q be a physical local angle direction, S its conjugate reaction,
B=Omega(Q,S). Suppose actual tube centroid and symmetric tracefree velocity
moments of S vanish, while its ACTUAL spin satisfies A_D(S)=B n. These are
the rows independently constructed in 0080, not an identity inferred from
the global impulse. The isotropic trace of the velocity first moment drops
out of incompressible affine probes.

For a locally affine divergence-free U and beta=curl U/2, the exact Euler
pairing and the tube multipole identity give

    Omega(U,S)=integral U.v_S rho,
    <U,p_H(S)>=<U,p_E(S)>-beta.A_D(S)
                [with the retained shape term if its row is nonzero].

The hybrid mean response therefore differs from the point-mean response by
the ACTUAL tube spin, not by a chosen global-angular-impulse label.

## 4. The physical fixed-angle affine lift

At the affine reference tangent use

    Gamma_U=U-sum_a (n_a.beta_a) Q_a,
    xi=Gamma_U+sum_a (n_a.Phi_a) Q_a
       =U+sum_a n_a.(Phi_a-beta_a) Q_a.

Q has zero material centroid displacement and vanishes near the tube
interface. Thus Gamma_U gives the actual centroid displacement U(X_a),
the correct affine interface displacement, and holds the physical core
angle fixed when U varies at fixed Phi. The physical local angle read from
the field is Phi, while the energy deformation relative to the transported
frame is Phi-beta. This is an actual tagged displacement lift, not a change
of probability-law labels.

Using B=Omega(Q,S)=n.A_D(S),

    Omega(Gamma_U,S)=Omega(U,S)-beta.B
                    =<U,p_H(S)>

for the affine probes just specified. Nonaffine gradients carry the full
Taylor/return terms, rather than replacing this affine equality by an
unproved all-k identity. The constant translation, strain and interface
terms are retained separately if their prescribed rows are nonzero.

This is the crucial difference from 0072. There the retained mean coordinate
was the whole-fluid POINT mean and its induced momentum equaled the full
Omega(U,S). Here the physical coordinate is the hybrid centroid/ambient
coordinate, and its conjugate mixed row is Omega(Gamma_U,S). The difference
is the actual tube moment measured in section 3.

## 5. Exact hybrid centering of the complete action

The common-translation ensemble is specified in the exact material chart:
tube centroid velocities and continuous ambient velocities have a coherent
slow component U_Hdot, plus independent phase-relative components whose
weighted HYBRID momentum has zero mean. All relative components, within-tube
shape velocities and pressure/shared-interface constraints are retained.
At zero gradient the coherent mass metric is rho. At gradient order let G
denote its FULL computed metric, including within-tube affine motion and
the fixed-angle lift corrections. Its gradient coefficients are not copied
from an unrelated point-mean normalization.

Write s for the retained internal fluid reaction momenta and z for all
other retained physical shape/interface coordinates. Let R s be their
induced hybrid momentum and D the actual tube-angle KKS coupling. Additional
shape momenta and their induced means can be included by enlarging s. The
uncentered affine-order cotangent/kinetic blocks have the form

    Theta=<P0+R s,dU_H>+<s,D dPhi+D_z dz>+Theta_remaining,
    H=<P0,G^-1 P0>/2+<P0,G^-1 R s>+H_internal.

The SAME R appears in both equations. In particular the velocity used in
the kinetic cross is the hybrid response of Gamma_U, not the point-mean
response of U substituted from 0072. A mismatch here would repeat the exact
double-counting error already diagnosed there.

Changing to physical hybrid momentum p_H=P0+R s gives exactly

    Theta=<p_H,dU_H>+<s,D dPhi+D_z dz>+Theta_remaining,
    H=<p_H,G^-1 p_H>/2+H_internal-<s,R*G^-1 R s>/2.

Eliminating p_H therefore cancels the residual hybrid induced-mean cross,
but leaves D Phidot. It does NOT subtract another B betadot: that difference
has already been accounted for by the physical hybrid moment and Gamma_U.
The reaction operator is P_H=P-R*G^-1R, with the complete analogous block
Gram subtraction for coordinate-dependent induced means. Its inverse is
the FULL fluid operator, never an isolated-cell inverse.

With independently varied time-reversed reactions, the even reduced kinetic
block on (Phi,z) is

    J=[D,D_z]* P_H^-1 [D,D_z].

All mixed shape terms remain. When the physical coordinates form an independent
admissible chart and P_H is coercive on its reaction space, this is positive
on their independent source image. An additional legitimate Routh/shape
reduction uses its full Schur complement. A nonzero D alone does not justify
discarding z or treating a duplicated shape/gauge direction as a new spin.
The physical tube angular-momentum response is the same D* s, so its reduced
law follows from this full block, not from an assigned rigid mass.

Consequently exact hybrid centering does NOT force 0072's point-mean
cancellation. Under a completed admissible two-field reduction with
J_PhiPhi=j>0 and no unresolved shape-rate cross, the physical leading action
has rho|U_Hdot|²/2+j|Phidot|²/2 and the objective locking in Phi-beta. Its
centroid transfer is nonzero. This is a conditional implication from the
actual material lift and full reduction, not a claim that moment rows alone
have already supplied those inputs.

## 6. The material reconstruction block that remains essential

The following concrete issue prevents promoting the tangent algebra alone
to a closed two-field Euler parcel action. A compact isovortical generator
can vanish near the tube interface, giving delta chi_D=-xi.grad chi_D=0,
while v_xi=P(xi cross omega) has nonzero normal trace there. A trajectory
containing ONLY those compact Q/S configuration directions would freeze the
material interface while its proposed physical Euler velocity crosses it.

Therefore z includes the interface embedding and its pressure response,
together with all needed shape/circulation reconstruction directions. The
admissibility identity is the actual material one

    (partial_t delta g).n = v'.n + background transport terms

on the moving interface, equivalently the linearized material tag equation.
It is not replaced by a finite list of spin and centroid moment constraints.
The full material action in section 1 enforces this identity; any further
two-field reduction must either provide a compatible lift or retain/solve
this boundary block. Independent ambient phase-relative velocities are part
of that block, not fluid that can be silently removed from the mass count.

This is not a request for an invariant unrestricted Euler finite-dimensional
ansatz. It is the kinematic requirement that the claimed fluid velocity and
the claimed material boundary belong to the SAME constrained material path.
The original affine conditional ensemble can prescribe the compatible shape
map, but cannot simultaneously prescribe a fixed tag and nonzero normal
crossing velocity. The current rank construction supplies the moment rows;
it does not by itself close that boundary reconstruction.

## 7. Executed tag-transport repair and its exact locality criterion

Treat the tube indicator as an actual advected observable, not a frozen
Eulerian mask. If zeta is its normal displacement at the stationary invariant
surface Sigma, the linearized material tag equation is

    partial_t zeta+div_Sigma(u0_t zeta)=v'.n.

It follows from incompressibility and the moving-surface transport theorem.
In particular the normal strain term is retained inside the surface-density
divergence. Integrating gives the correct volume constraint
`d/dt integral_Sigma zeta dA=integral_Sigma v'.n dA=0`.

The material boundary spin variation is exactly

    delta S_boundary=rho integral_Sigma r cross u0 zeta dA,

in addition to the bulk moment. This term need not vanish when its bulk
spin row has been prescribed. Its phase and amplitude are fixed by tag
transport, not by the scalar identity A_D(S)=B.

Let T_t be the transport semigroup of the surface-density operator. For
initial displacement zeta_0 the exact finite-time solution is

    zeta(t)=T_t zeta_0+integral_0^t T_(t-s) [v'(s).n] ds.

This is a constructive compatible material tag with retained initial data;
no frequency-space inverse is needed to define it. On an invariant torus
whose stationary tangential flow is smoothly conjugate to constant frequency
vector Omega, write zeta dA=b(theta)dtheta. Then

    bdot+L b=f,   L=Omega.grad_theta,   f=J(theta) v'.n.

The factor J is the actual surface area Jacobian. Using the density b avoids
assuming an invariant Euclidean surface area element or a full integrable
foliation around the torus. The semigroup statement is valid even when that
particular torus conjugacy has not been supplied.

Now retain the full two-field tangent velocity, schematically
`v'=V1 Phidot+V0 Phi` for one angle channel; the same calculation applies
columnwise to displacement and other shape channels. Put
`f1=J V1.n`, `f0=J V0.n`. The exact decomposition is

    b=f1 Phi+w,
    wdot+L w=(f0-L f1)Phi,
    w(0)=b(0)-f1 Phi(0).

Consequently a LOCAL instantaneous collar tag `b=f1 Phi` exists if and only
if `f0=L f1` and the initial tag is compatible, for arbitrary Phi histories.
This is a derived local-action compatibility condition, not an added demand
for unrestricted Euler finite-dimensional invariance. It directly prevents
the previously exhibited mismatch between a stationary mask and crossing
fluid velocity.

The same calculation executes the repair when that condition fails: retain
w as the actual transported boundary coordinate with its pressure conjugate.
Its contribution to the physical spin is the finite-time weighted history
`rho integral (r cross u0) w dtheta`. It is not an unspecified missing term.
For a harmonic history exp(-i omega t), the m-th torus Fourier coefficient
in the steady particular solution is

    b_m=[-i omega f1_m+f0_m]/[i(Omega.m-omega)] Phi
       =f1_m Phi+[f0_m-i(Omega.m)f1_m]/[i(Omega.m-omega)] Phi.

Thus the boundary-spin observable is genuinely frequency dependent whenever
its weight has a nonzero projection on a nonzero residue in this expression.
Transported homogeneous initial data are retained separately. If a spin
weight annihilates every residual projection, that particular observable
can be local even when the full tag is not; tag memory alone is not treated
as a proof of nonlocality of every moment. At a forced resonance the finite-
time solution has secular growth rather than a finite constant modulus.

An exposing exact example is one torus harmonic with f1=cos(theta), f0=0
and spin weight cos(theta). Its frequency multiplier is
`omega²/[2(omega²-Omega²)]`. Replacing the missing commutator column by
`f0=L f1=-Omega sin(theta)` instead gives exactly 1/2 for every nonresonant
frequency, with the apparent poles removable. The next paragraph explains
why that repair is an actual material displacement construction, not a
chosen observable coefficient.

For a divergence-free material collar displacement xi=Xi Phi, the exact
Eulerian velocity relation is

    v'=Xi Phidot+[u0,Xi]Phi.

Its transported surface normal density satisfies f0=L f1 automatically.
Choosing Xi with prescribed zero-total-flux normal trace is possible by a
surface Hodge potential and a smooth collar extension of its vector
potential; the curl of that extension is divergence free and can vanish
near the observed core. It therefore repairs the tag kinematics without
changing the core angle jet. Its volume, boundary spin, geometric kinetic
and pressure-Hessian terms are all part of the full material action.

This collar changes the trial displacement and thus its coefficients.
Fixed Kelvin circulation additionally imposes the actual material
circulation-momentum condition; it does not follow from the scalar normal
transport identity alone. The coordinated construction in 0084 owns this
collar/Jacobi and Kelvin-reconstruction block. The present calculation has
executed its exact transport/observable half, identified the precise
commutator column needed for locality, and supplied the complete history
response when a local collar has not yet been selected.

The actual material lift also makes the full spin local with an exact
identity supplied and derived in that coordinated calculation. At the
stationary invariant reference parcel, integration by parts using u0.n=0
and integral_D u0=0 gives

    delta S_D=rho integral_D [r cross xi_t+2 xi cross u0].

This equals the bulk velocity-spin variation PLUS its moving-boundary
term. For xi=Xi_a z_a it is A zdot+C z, with
`A_a=rho integral r cross Xi_a` and `C_a=2rho integral Xi_a cross u0`.
The paired opposite stationary flows leave A unchanged and reverse C.
It is therefore possible to obtain a genuinely local physical spin law
from a compatible material lift. Its A and complete Jacobi metric must be
computed from that lift; they are not assumed equal to the old coadjoint
bulk-spin rows of Q/S when the collar or circulation reaction changes them.

## Route result

Established: exact hybrid momentum/multipole map; the affine fixed-physical-
angle KKS identity using ACTUAL tube spin; and the complete same-R
cotangent/kinetic hybrid centering. This theorem shows why physical tube
spin can survive hybrid centering without contradicting 0072.

Established also: the exact transported material-tag solution, its physical
boundary-spin history, and the necessary-and-sufficient local collar
criterion f0=L f1 for the prescribed velocity columns. The incorrect frozen-
tag subroute is replaced by that transport construction. The remaining
positive two-field action construction is its actual collar/Kelvin/full-
Schur implementation in 0084, not a missing formula for the tag response.
No extra rigid inertia, missing ambient mass, isolated inverse or all-k
closure has been used. These route-scoped results leave the original
objective active and supply executable boundary data for the next lift.
