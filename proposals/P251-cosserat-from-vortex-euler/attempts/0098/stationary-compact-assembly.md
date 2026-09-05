# One stationary Euler law, finite compact pairs, and two different momenta

## 1. Freeze a uniformly bounded positive-probability event

Take0085's actual analytic EPS Beltrami prototype with its invariant solid
torus D and its finitely many disjoint interior support balls. The final core
pair, any finite number of gradient/strain pairs, the selected nonzero minors,
cutoffs, amplitudes and finite carriers are fixed before defining the event.
Every positive matrix margin is strict; every support has positive distance
from the tube boundary and every other support. Choose a finite derivative
order m large enough to bound all determinant-formula derivatives appearing
in the forces, their curls and the required moment estimates.

The construction is continuous in the C^m field jet on these fixed balls:
inversion denominators stay bounded below, the finite coefficient norms stay
bounded above, and the strict Hessian, KKS and core-observation margins persist
in a sufficiently small neighborhood. Use the existing EPS persistence
neighborhood to keep the selected tube and its supports. Intersect these
finitely many neighborhoods. Its radius is strictly positive. In particular,
one does not select carriers afresh near arbitrarily singular random minors.

The isotropic Gaussian Beltrami law constructed in0071 has positive probability
in this compact-open neighborhood. Denote that probability by p_good>0. The
local reconstruction uses the same pivot indices and cutoffs as the prototype;
its coefficients are continuous finite-jet functions on the event. Thus it is
measurable. The actual invariant domain can also be selected measurably, as
the following construction shows; a domain mark alone would not prove that.

Fix an integer b>=3, large enough for the boundary observations used here,
and a reference torus embedding f0. The sourced persistence theorem
(1210.6271, Theorem7.6 and its construction of the spatial torus) applies
with b+1 and a prescribed small embedding bound. Let K be the C^b closure
of that bounded C^(b+1) neighborhood of f0. It is compact by the derivative
equicontinuity/Arzela--Ascoli argument. Choose the neighborhood small enough
that every embedding in K bounds a solid torus, contains the fixed response
supports with a uniform margin, and stays inside the coherence ball. Its
topology and a positive lower volume bound are then fixed. The boundary
selected below has C^b regularity; the Euler field and interior profiles
remain C-infinity. No boundary analyticity is needed for their moment or
pressure-flux formulas.

For a field u in the good event, the compact fiber

    F(u)={f in K: u(f).(partial_1 f cross partial_2 f)=0}

is nonempty by persistence. The displayed residual is continuous in u and
f, and its zero set is closed. Cover K at stage n by finitely many closed
balls of radius2^(-n), and select the first ball whose intersection with
the previously selected compact set and F(u) is nonempty. Each test is
measurable: projection of this closed condition through the compact f
factor is closed. The selected nested compact sets have diameters tending
to zero and yield one f(u); the corresponding ball centers converge to
it measurably. This is an explicit measurable selection, not a premise
about random marks. The enclosed domain, volume and polynomial moments
depend continuously on these C^b embeddings. Implementing the selection
in the auxiliary local frame preserves simultaneous rotation covariance.

Assign an independent proper orthonormal frame to each Poisson candidate.
Apply the fixed prototype rule in that frame. This removes the apparent
preferred axes of pivot selection without changing the Gaussian field law.
These independent local frames are not a single globally rotated anisotropic
sample. Translation and simultaneous rotation of field and candidate frame
commute with the reconstruction.

With candidate intensity tau and isolation radius2R, the retained intensity is

    nu=tau exp[-tau Vol(B_(2R))] p_good>0.

Every retained R-ball is disjoint from the others. All patches are restrictions
of ONE smooth stationary Euler field, not independently glued solutions.
The Gaussian/Poisson law is stationary, isotropic and mixing as in0071, so its
marked factors are ergodic. Each actual field is stationary Euler; the
statistical assertions are about its spatial marked law. Marks are frozen at
the reference state and transported in the material variation, not reselected
to minimize deformed energy.

The finite bounds yield, for example, P>=p_*>0, |B|>=b_*>0,
H_QQ-N^2/P>=k_*>0 and finite upper bounds. Thus

    j_bar=nu E_Palm[B^2/P]/3>0,
    kappa_bar=nu E_Palm[H_QQ-N^2/P]/3>0

are finite for the isotropic scalar-axis population. These are actual
Euler-integral coefficients and declared geometry/law inputs, not a fitted
frequency or the former filament formula. Time-reversal paired actions are
reduced before the pairing, so their odd current term is handled as in0085.

Locality here is special and exact: for any two distinct compact response
supports, BOTH their induced velocities and generators are disjoint. Both
terms of the complete coadjoint Hessian and their KKS cross therefore vanish.
This justifies the local internal blocks, but does not remove a genuine
macro/shape/circulation reaction block; that identification belongs to0097.

## 2. Keep the complete physical hybrid observation

Use0082's actual tube-centroid distribution plus continuous ambient momentum

    p_H=sum_a P_a delta_(X_a)+rho 1_A u.

There is no ambient centroid. Whole-fluid point momentum is p_E=rho u.
For a compact induced-velocity profile v supported strictly inside a selected
tube, integration by parts gives exactly

    integral v_i=0,
    integral(r_j v_i+r_i v_j)=0,
    S=rho integral r cross v.

The actual tube mean response and ambient response vanish, hence p_H(v)=0
for that internal profile. This holds at every slow amplitude assigned to
the patch; it is not a small-wave-number approximation. The same statement
is false for p_E(v), and that distinction is retained.

For Fourier convention exp(-i k.r), define the exact profile response

    V_i(k)=rho integral v_i(r) exp(-i k.r) dr,
    T_ijl=rho integral v_i(r) r_j r_l dr.

Taylor expansion and the antisymmetric first-moment identity give

    V(k)=i k cross S/2 - (k_j k_l/2) T_ijl + R_3(k),
    |R_3(k)| <= rho |k|^3 integral |v(r)| |r|^3 dr/6.

This is an observable identity with an explicit finite remainder, not an
equation of motion. Campbell averaging preserves the estimate because all
profile moments are uniformly bounded on the good event. The exact Fourier
integral is retained when an observable requires terms beyond the displayed
jet. Slow amplitudes are evaluated at the patch centroid; they are not
multiplied through a divergence-free profile in microscopic space without
its solenoidal correction.

The base invariant-domain identity gives the same tensor relation from
actual Reynolds transport. For these compact boundary-compatible tangents
the symmetric moment and the linear shape-rate term vanish. If an additional
reconstruction direction moves the boundary, its Reynolds/shape term from0082
is added rather than assigned to this compact-profile formula.

## 3. Parity and time reversal have different jobs

Under spatial inversion r->-r a polar velocity transforms to -v(-r), while
the angular momentum S is axial. Consequently its first-moment spin stays
fixed and its second velocity moment T changes sign. Pair the full marked
law with its inverted law, including the reversed curl eigenvalue, to cancel
T in the coherent isotropic axial response. A single fixed-helicity sample
is not claimed to contain its parity partner. The remainder remains O(k^3).

Time reversal instead changes the background velocity and KKS sign. Its
action-level pairing cancels the static angle-spin term -BN q/P in0085
after reaction elimination. Assuming the completed action supplies the
derived leading actual spin S_bar=j_bar Phi_dot, the compact part of the
physical observation therefore reads

    rho U_Edot=rho U_Hdot+curl(j_bar Phi_dot)/2+O(partial^3)

with all other reconstruction/shape observations separately retained.
For constant j_bar and prepared zero integration constant,

    U_E=U_H+(j_bar/(2rho)) curl Phi+O(partial^3).

The two fields describe different physical measurements. An action expressed
in U_H earns a centroid/ambient translation claim; substituting U_E requires
transforming BOTH kinetic and potential forms with this same map. In
particular the induced order-k^2 kinetic terms are not optional. No vanishing
point-mean response can refute a nonzero centroid response without this map,
and no canonical coordinate may be called either observation without it.

## 4. Scope and next action

Established as stated: bounded stationary positive-density compact-pair
assembly and its exact finite-radius physical momentum observation, using
the declared sourced law and local construction. The action-dependent spin
law in section3 is an explicitly typed implication; it is not a replacement
for0097's material/cotangent reconstruction. Continue that construction and
0096's complete second-gradient reduction before freezing the parent claim.
