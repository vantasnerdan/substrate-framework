# Literal Euler torque, its variation, and the curvature representative

This construction concerns the actual material-tag/continuous-ambient
coarse-graining, not a finite-wall vortex model. Density rho is constant;
p is physical pressure, not pressure divided by rho, so rho D_t u=-grad p. All fields and
material boundaries have the regularity needed for the displayed
integrals; distributional indicators are permitted. Time intervals and
the averaging convention are fixed. No constitutive modulus is inserted.

## What the parent actually asks

The original P251 N3 license specifies a periodic box; N4 specifies
C² periodic fields and compactly supported variations. The question
asks for isotropic constitutive relations and angular-momentum/couple-
stress balance. It does not prescribe a free-boundary couple-traction
experiment. Accepted C-CST-009 likewise states the second-gradient
conditional action and bulk balance, retaining physical current maps.

Consequently the one-dimensional local curvature ambiguity identified
by 0227 is not by itself a missing third *bulk* physical modulus. It is
important to exhibit a stress representative and the boundary/current
change between representatives. It would be a different statement to
predict the literal traction on a specified free boundary without that
change. This proof makes that distinction constructive rather than
assigning the undetermined coefficient a desired numerical value.

## Exact tagged force and mechanical spin

Let chi_a=1_{D_a(t)} be disjoint material indicators and
chi_0=1-sum chi_a the complete continuous ambient complement. Each obeys
partial_t chi+u.grad chi=0. Define

    M_a=rho integral chi_a,
    X_a=(rho/M_a) integral chi_a y,   V_a=Xdot_a,
    P_a=M_a V_a,
    S_a=rho integral chi_a r_a cross (u-V_a), r_a=y-X_a.

Then M_a is constant, integral chi_a r_a=0, and Reynolds transport and
Euler give EXACTLY

    Pdot_a=F_a=-integral chi_a grad p=-integral_boundary p n_a,
    Sdot_a=tau_a=-integral chi_a r_a cross grad p
                =-integral_boundary r_a cross (p n_a).           (1)

The second identity uses D_t r_a=u-V_a. The term
rho integral chi_a r_a cross Vdot_a is zero by centering, not by a
rigid-tag approximation. There is no advective exchange across a
material boundary, but pressure traction is retained there.

The corresponding local point identity makes this distinction visible.
With ell_a=rho chi_a r_a cross (u-V_a),

    partial_t ell_a+div[ell_a tensor u+chi_a p (r_a cross I)]
      =p r_a cross grad chi_a-rho chi_a r_a cross Vdot_a.  (2)

Here (r cross I) has column j equal to r cross e_j. Since
grad chi_a=-n_a delta_boundary, integrating (2) gives (1). The companion
point momentum identity is

    partial_t(rho chi_a u)+div[chi_a(rho u tensor u+p I)]
       =p grad chi_a.                                  (3)

Thus pressure and its full tag-boundary reaction are present before
any constitutive approximation or spatial expansion.

## A fixed hybrid stress and couple-stress representative

Use the ACTUAL hybrid momentum and collapsed intrinsic spin

    J_H=sum P_a delta(x-X_a)+rho chi_0 u,
    S_H=sum S_a delta(x-X_a).

No ambient parcel mass is omitted or replaced by a filling fraction.
Its transport fluxes are

    T_H=sum P_a tensor V_a delta_a+rho chi_0 u tensor u,
    C_H=sum S_a tensor V_a delta_a.

For each point y of a tag boundary, let the force measure on the tag be
f_a(dy)=-p(y)n_a(y)dS_y, let r=y-X_a, and set
delta_s=delta(x-X_a-sr). Fix the straight segment from X_a to y as the
localization convention. Define the exact distributions

    sigma_* = sum integral_boundary integral_0^1
                   f_a(dy) tensor r delta_s ds,
    mu_* = sum integral_boundary integral_0^1
                   (1-s)[r cross f_a(dy)] tensor r delta_s ds,
    sigma_H=-chi_0 p I+sigma_*.                         (4)

These formulas do not assume the tags are convex or impose a physical
force along a bond. The segment localizes the already computed boundary
force. On the periodic universal cover, choose the boundary and centroid
images consistently and periodize the resulting distributions.

Because r.grad_x delta_s=-partial_s delta_s,

    div sigma_* =sum[F_a delta_a-f_a],
    div mu_* =sum tau_a delta_a
               -sum integral_boundary integral_0^1
                       (r cross f_a) delta_s ds.        (5)

For tensor convention (force component, flux direction), define
ax(sigma)_i=epsilon_ijk sigma_jk. Then
ax(f tensor r)=f cross r=-r cross f, so (5) implies

    partial_t J_H+div T_H=div sigma_H,
    partial_t S_H+div C_H=div mu_*-ax(sigma_H).          (6)

The ambient equation (3) with chi_0 has boundary source
p grad chi_0=-sum f_a. It cancels precisely the endpoint forces in
(5). This is why the continuous ambient cannot be silently removed.
Shared tag faces also cause no problem: their opposite force measures
cancel at the physical endpoint even if no open ambient layer lies
between them.

Finally x cross div sigma=div(x cross sigma)+ax(sigma). Adding orbital
and intrinsic angular momentum in (6) cancels the antisymmetric force
stress and gives the exact total conservation law. Convolution with a
fixed smooth coarse kernel commutes with all distributional derivatives,
so (4)--(6) immediately supply a smoothed physical representative.
Changing localization paths or the orbital/intrinsic bookkeeping can
change a representative, but not its retained total virtual work and
balance. A comparison region cutting a tag measures different point-fluid
and centroid-localized angular momentum: its exact cut-cell correction,
or the corresponding controlled finite-radius moment expansion, accompanies
the flux comparison. The hybrid traction is not silently identified with
the bare point-fluid pressure traction on that cut. No symmetry of the
microscopic pressure stress prohibits the
asymmetric hybrid stress or nonzero couple stress in (4).

## Actual Euler/Lin variation, including the moving boundary

Let xi be an actual volume-preserving material displacement,
w=delta u the Eulerian velocity perturbation, and p_1=delta p. Their
kinematic and dynamical equations are

    div xi=div w=0,
    w=xi_t+(u.grad)xi-(xi.grad)u,
    rho[w_t+(u.grad)w+(w.grad)u]=-grad p_1,
    Delta p_1=-2rho tr[(D u)(D w)].                    (7)

The pressure equation carries its actual ambient boundary/periodic
condition, not an independently adjustable normal traction. The
indicator and centered observables vary as

    delta chi=-xi.grad chi,
    delta X=(rho/M) integral_D xi,
    delta V=(rho/M) integral_D [w+(xi.grad)u]=delta Xdot,
    delta S=rho integral_D [
       (xi-delta X) cross (u-V)
       +r cross (w+(xi.grad)u-delta V)].               (8)

Equation (7) also implies the actual tag transport equation
(partial_t+u.grad)delta chi=-w.grad chi. It is not legitimate to keep
the tag fixed merely because a chosen instantaneous orbit generator
vanishes near its original boundary. Formula (8) retains the material
domain, centroid and velocity changes simultaneously; it applies even
when the reference tag is moving and noninvariant.

For a material surface the area-vector variation is

    delta(n dS)=[(div xi)I-(Dxi)^T]n dS=-(Dxi)^T n dS.

Define the material traction-variation tensor

    Q=[p_1+xi.grad p]I-p(Dxi)^T.                       (9)

Then delta f=-Q n dS. The complete force and torque variations are

    delta F=-integral_boundary Q n
           =-integral_D [grad p_1+(Hess p)xi],
    delta tau=-integral_boundary [
          (xi-delta X) cross p n+r cross Q n]
       =-integral_D [(xi-delta X) cross grad p
                     +r cross(grad p_1+(Hess p)xi)].   (10)

For the first equality, div Q=grad p_1+(Hess p)xi because div xi=0.
For the second, the derivative of the lever arm and the antisymmetric
part of Q cancel exactly. Dropping either the moved normal or the
lever-arm variation breaks this identity. There is no unexplained
boundary spin or shape term hidden in (10).

The complete stress variation follows without inventing a constitutive
response. In (4), vary the force measure by delta f above, the lever
arm by delta r=xi(y)-delta X, and each segment distribution by

    delta delta_s=-[delta X+s delta r].grad_x delta_s.  (11)

Apply the product rule to f tensor r for sigma_* and to
(1-s)(r cross f) tensor r for mu_*. Also vary -chi_0 p I and the
actual transport fluxes in (6). Equations (9)--(11) are an explicit
linear boundary/current supplier for any constructed xi,w,p_1; they
retain all endpoint and path-motion terms. In particular they are not
the fixed-domain integral of w alone or a renamed canonical momentum.

## The one-dimensional isotropic gradient freedom

Let G_ij=partial_j Phi_i and use the canonical curvature density

    W_curv=c_tr(tr G)²+c_s||sym G||²+c_a||skew G||².

Its variational couple stress, in the same flux-index convention, is

    m=2c_tr tr(G)I+2c_s sym G+2c_a skew G.              (12)

The bulk coefficients are gamma_T=c_s+c_a and
gamma_L=2(c_s+c_tr). If both are fixed, any two such isotropic local
representatives differ by exactly

    (delta c_s,delta c_a,delta c_tr)=(e,-e,-e),
    delta W=e[tr(G²)-(tr G)²],
    delta m=2e[G^T-tr(G)I].                            (13)

For compatible G=grad Phi, the divergence is identically zero:

    partial_j delta m_ij=0,
    delta m_ij=partial_k A_ijk,
    A_ijk=2e(delta_ik Phi_j-delta_ij Phi_k),
    A_ijk=-A_ikj.                                     (14)

This is an explicit couple-current superpotential, not a new restoring
force. No change to the force-stress axial torque is needed for this
particular improvement. The energy is also an explicit divergence,

    delta W=e div[(Phi.grad)Phi-Phi div Phi].           (15)

For a smooth bounded comparison region O, its exact virtual-work change
is

    delta_variation integral_O delta W
      =integral_boundary(O) delta_variation Phi . delta m n. (16)

The left side is the variation of the boundary functional furnished by
(15). It is generally NOT zero. On a periodic box or under compactly
supported variations it is zero, which is precisely the original
N3/N4 action scope. For a free boundary, replacing the bulk density by
W+delta W while subtracting that SAME boundary functional preserves the
complete action and physical boundary virtual work. Replacing the bulk
density alone changes a traction experiment and is not an equivalence.

The representative from (4) records one literal physical localization
of the torque. Equations (13)--(16) explain why matching its bulk balance
to (12) cannot identify e: divergence annihilates that row. If a future
claim asks for the traction associated with a fixed boundary/localization
convention, (9)--(11) supply the additional response that must be
evaluated on its actual Euler preparation. One does not obtain that
response from the optical dispersion or prescribe e as a fitted modulus.

## Positive conclusion and remaining supplier interface

The actual Euler force/torque and full boundary variation give an exact
hybrid stress/couple-stress representative with all ambient reaction.
The original periodic action and angular balance determine the
isotropic curvature energy as an equivalence class with the explicit
one-dimensional current/boundary improvement (13)--(16). This is a
complete answer to the boundary-representative issue in 0227 at that
scope: it adds no unnecessary third-bulk-modulus requirement.

It does not compute an unprovided free-boundary traction coefficient,
nor replace the actual optical/acoustic/current supplier of the parent.
Any derivative field normalization used in that supplier must transform
its physical current and boundary functional together. The actual
coarse constitutive closure still comes from its same Euler action;
the exact torque identities alone do not manufacture that closure.

Route A: established literal localization and its complete variation.
Route B: established periodic/current-boundary equivalence. Evidence
scope: exact calculus and variational identities; not a new microscopic
modulus or a completed Euler/Cosserat/EPS construction. The next parent
achievement remains the active physical supplier, with no additional
free-boundary experiment inserted into its original success contract.
