# Exterior-fluid response theorem and exact symplectic centering

## 1. The actual field and physical mean functionals

Use the actual finite-spherical-harmonic EPS field of 0048, with
omega=lambda u, lambda nonzero, and its closed invariant material tube D.
Put X=|D|^-1 integral_D x and K=e cross (x-X), with nonzero rotation tangent
vK. The base mean vanishes: integral_D u=0 by u.n=0. Rotate D together with
the field for the common direction, so its complete mean variation is zero.

For compact generators supported strictly inside OR strictly outside D,
the geometric domain does not change. Its physical mean-velocity functional
and common angular moment are

    c_i(xi)=integral_D P(xi cross omega)_i,
    l(xi)=Omega(K,xi).

Use the actual R³ Leray projection, and define

    w_i=P(1_D e_i),
    f_i=omega cross w_i  (i=1,2,3),
    f_0=-rho vK.

Then `c_i(xi)=integral xi.f_i` and `l(xi)=integral xi.f_0`. The second
identity follows because K cross omega-vK is a global gradient and xi is
compact solenoidal. All constructed generators below are curls of smooth
potentials compactly supported away from the tube boundary. This avoids
assuming that a curl-free field on a solid torus has zero period.

## 2. Exact independence in every exterior open ball

The connected exterior R³\closed(D) is the domain of this argument. Each
w_i is harmonic, curl-free and analytic there. For a constant vector a,
write w_a=sum a_i w_i. The dipole expansion of the Newton potential gives

    w_a(r n)=|D| [3 n(n.a)-a]/(4 pi r³)+O(r^-4),
    n.w_a=|D| (n.a)/(2 pi r³)+O(r^-4),
    grad w_a=O(r^-4).

The signs follow from `P f=f-grad Delta^-1 div f` and
`Delta^-1 1_D=-|D|/(4 pi r)+O(r^-2)`. The leading coefficient depends only
on the positive volume of D, not on a spherical approximation to its shape.

The actual EPS construction in 1210.6271, Theorem 8.3, is a finite sum of
regular spherical Bessel modes, followed by a curl polynomial. Thus omega
has an expansion

    omega(r n)=[A(n) sin(lambda r)+B(n) cos(lambda r)]/r+O(r^-2),
    partial_r omega=lambda[A(n) cos(lambda r)-B(n) sin(lambda r)]/r
                    +O(r^-2).

A and B are smooth finite angular sums and are not both zero for a nonzero
field. To see the last statement without an assumed radiation condition,
expand each Cartesian component into regular `j_l(lambda r)Y_lm` modes.
Its sine/cosine leading terms cannot all vanish unless every coefficient
vanishes, by spherical-harmonic orthogonality. The same statement applies
to vK: rotations and translations preserve finite regular Helmholtz sums,
and vK is nonzero. Therefore curl f0=-rho lambda vK has a nonzero 1/r
far-field coefficient and cannot be O(r^-4).

For the mean response, divergence-free omega and w_a give

    curl(omega cross w_a)=(w_a.grad)omega-(omega.grad)w_a
      =lambda |D| (n.a)
          [A(n) cos(lambda r)-B(n) sin(lambda r)]/(2 pi r⁴)
          +O(r^-5).

Here the transverse derivative in (w_a.grad)omega carries an extra 1/r;
the displayed leading term is the radial derivative. If this curl vanished
throughout the exterior, its sine and cosine coefficients would imply
`(n.a) A(n)=(n.a) B(n)=0`. For a nonzero a, n.a is nonzero on an open dense
set, forcing A=B=0 by continuity. That contradicts the nonzero EPS field.
Consequently the three mean-response curls are linearly independent.

Now suppose a linear combination of curl f0,curl f1,curl f2,curl f3
vanishes on any nonempty exterior open ball. Analytic continuation makes it
vanish in the connected exterior. The 1/r versus 1/r⁴ orders first force
the coefficient of curl f0 to vanish. The preceding dipole argument forces
the remaining three coefficients to vanish. Thus all FOUR response curls
are independent on every such ball. This proves the formerly missing
centered common-momentum independence; it is not a new rank hypothesis.

## 3. Four compact responses with an isotropic KKS span

Choose four pairwise disjoint bounded exterior balls, separated from the
physical core jets and from the raw cages to be used below. In ball j choose
a nonnegative smooth bump chi_j, positive on an open subset, and form

    G^(j)_{ab}=integral chi_j curl f_a.curl f_b,  a,b=0..3,
    g_b^(j)=curl(chi_j curl f_b).

The theorem above makes every G^(j) strictly positive definite. These are
finite exact integrals, and there is no small-eigenvalue numerical decision.
Define ONE response in each ball by

    eta^(j)=sum_b g_b^(j) [(G^(j))^-1]_{bj}.

Integration by parts gives the exact right-inverse property

    F_a(eta^(j))=delta_aj,  F=(l,c_1,c_2,c_3).

Units may be fixed once when normalizing the four dual coordinates; that
normalization is a coordinate choice and not a physical inertia or modulus.
Most importantly, their mutual KKS pairings vanish:

    Omega(eta^(i),eta^(j))=0 for all i,j.

For i different from j this follows by disjoint support; for i=j it follows
from antisymmetry. Therefore their span is isotropic. This additional support
construction is what makes the following moment projection symplectically
controlled; a generic Gram projection alone would not do that.

Let `R F(xi)=sum_j eta^(j) F_j(xi)` and `Pi_F=I-R F`. These are fixed
linear maps computed from the actual field and D. Then F(Pi_F xi)=0 exactly.

## 4. Centered core/cage fields and the FULL KKS matrix

Choose the two opposite physical rotation jets Q_R of 0048 inside D.
Choose its circular internal cage pair C1_k,C2_k in a region disjoint from
those jets and all four response balls, with omega_z of one sign. Choose
a single negative-helicity body cage A_k in a different region, disjoint
from the internal supports and the response balls. Every selected support
lies in one finite coherence region of the same smooth Euler field.

Use the physical curl-generated formula of 0045 and k/lambda>0. With a fixed
nonzero normalization b0, set

    r=b0 eta^(0)+Pi_F A_k,
    Q=Pi_F(Q_R+C1_k),
    S=Pi_F C2_k.

One can take b0 from a positive source Gram entry; it simply fixes the
normalization of the body reaction coordinate. Then

    c_i(r)=c_i(Q)=c_i(S)=0,
    l(r)=b0,  l(Q)=l(S)=0.

All corrections are outside D, so both physical core jets are unchanged.
The actual section angles remain B+q and B-q. The response subspace is
isotropic and disjoint from every raw generator support. Raw body and
internal supports are disjoint. Expanding every cross term therefore gives

    Omega(r,Q)=Omega(r,S)=0,
    Omega(Q,S)=Omega(C1_k,C2_k)=:c_k.

There is NO asymptotic KKS correction here: the selected right inverse
preserves the internal pairing exactly. Its finite-k expression is the
0045 integral

    c_k/rho=(1-lambda/(2k)) integral omega_z phi²
             +k^-2 integral phi_z omega.grad phi.

Its explicit bound gives c_k nonzero at a finite carrier. Together with K
the final KKS form is exactly `b0 dB wedge dy+c_k dq wedge ds`.

## 5. Full positivity, with cell mean already excluded

Every F_j evaluated on a raw circular cage is O(|k|^-2), by integration by
parts in its compact carrier amplitude; all f_j are smooth on that support.
F_j(Q_R) is fixed. Hence the projection adds only uniformly bounded fixed
compact fields, with supports outside the physical tube.

The principal projected cage H is
`rho(1+|k|/|lambda|) diag(A_body,A_internal,A_internal)` with positive
A values. Sum the projection-error constants D,E and fixed-attachment norms
as in 0048 to obtain a finite C_H. The SAME argument gives

    H_compact >= rho[(1+|k|/|lambda|) A_min-C_H] I.

Increasing the carrier to a finite value makes this complete three-by-three
matrix positive. Its projected-energy off-diagonal entries are retained;
they have not been set to zero merely because the generators have disjoint
supports. The global K energy row is exactly zero by relative rotational
symmetry. All fixed response norms are finite because their supports are
bounded; the far-field theorem proved existence of the inverse but does not
put an infinite tail into a response generator.

For the selected physical cell D, each compact tangent now has
`delta V_D=c/|D|=0`, and the common direction has zero complete moving-domain
mean variation. The entire cell-centroid energy Gram in 0052 vanishes on
these directions. Thus the corrected rotor does not contain that cell's
translation energy. The same full KKS/H integrals feed the 0049 Routh API;
the independent-momentum time-reversal pair and the general physical-angle
map are then exactly those of 0048.

This statement centers D. It does not erase kinetic means or energy of
surrounding material parcels in the coherence region. Their pressure,
momentum, and assembly accounting remain explicit in the joint material
construction (0057); a whole-R³ Galilean boost is still not a finite-mass
direction of the EPS relative-energy orbit.

## 6. The ambient angular impulse is retained in the physical current map

The canonical common moment is the full relative impulse of the finite
coherence perturbation, not the bulk impulse of D alone. Split its compact
vorticity change into the tube and surrounding fluid:

    J_coh=-rho/2 integral_R3 |x-X|² e.delta omega=J_D+J_ext,
    J_ext=-rho/2 integral_R3\D |x-X|² e.delta omega.

For common rotations the integrals use the same relative prescription as
0048. For compact responses J_ext has bounded support, even though its
induced velocity has the legitimate Leray tail. The physical cell spin is

    delta L_D=J_coh+B_surface-J_ext.

Thus the precise current correction in 0052 is now
`B_effective=B_surface-J_ext`, not just B_surface. With
`A_ij=epsilon_ijk B_effective,k/2`, use
`P_can=P_centroid+div A`, `S_can=S_physical-B_effective` and retain
`sigma_can=sigma+partial_t A` (and the full convective flux bookkeeping).
The same exact angular-current improvement applies. The surrounding
angular response has been moved into an explicit orbital/couple current,
not discarded or called the physical spin of the core parcel.

For momentum tangent v_j, j=r,s, its coefficient is the exact finite row

    t_eff,j=(rho/2) integral_boundary(D) r² e.(n cross v_j)
              +(rho/2) integral_R3\D r² e.curl v_j.

After the independent ± momentum elimination, the time-even current response
is `t_eff P_momentum^-1 diag(b0,c_k) (Bdot,qdot)^T`. This is the explicit
observable filter carried into the material assembly and normal form. The
finite spatial coherence/assignment of the surrounding response is an
ensemble premise; its kinetic and moment integrals stay in the same action.

## Route verdict

**Established as stated, independent review pending.** Exterior fluid gives
an exact finite compact response right inverse, proves centered common
moment independence without a rank assumption, and restores the full KKS
matrix while preserving the physical jets and positive high-frequency
Hessian. The selected cell's translational mean is exactly absent from the
rotor, and the surrounding angular impulse is included in the current map.
The parent continues its one-action assembly and spatial-gradient joining;
this theorem does not claim those independent calculations are already done.
