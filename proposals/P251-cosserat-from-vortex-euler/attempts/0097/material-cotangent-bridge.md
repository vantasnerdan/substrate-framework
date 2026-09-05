# A single material action: the earned centroid block and reconstruction

This is a joining calculation, not a second construction of the compact
pair in 0085. In particular, an instantaneous moment identity is kept
distinct from the identification of the velocity of a material path.

## 1. Exact tagged material decomposition

Let the selected disjoint invariant tube parcels be `D_a`, with reference
mass measures `dm=rho da`, and let `A` be their entire ambient complement.
For an actual volume-preserving material map `g`, set

    M_a = integral_Da dm,
    X_a = integral_Da g dm / M_a,       r_a = g-X_a,
    P_a = integral_Da v dm,            pi_a = v-P_a/M_a.

Here `v` denotes material momentum per unit mass in the canonical action;
its equality to `gdot` is an equation of that action, not a definition
inserted after a phase-space restriction. The identities

    integral r_a dm = integral pi_a dm = 0,
    Theta = sum_a P_a.dX_a + integral_Da pi_a.dr_a dm
                              + integral_A v.dg dm,
    H = sum_a |P_a|^2/(2M_a) + integral_Da |pi_a|^2 dm/2
                              + integral_A |v|^2 dm/2

are exact algebraic changes of variables in the SAME material action.
Incompressibility, equality of neighboring traces, and pressure tractions
are still constraints on these variables. Thus these identities do not
assign each tube an independently movable boundary or an independent
canonical momentum after its ambient fluid has been eliminated. Such an
elimination generally changes the centroid metric through added-fluid
motion. The ambient integral cannot be omitted or counted a second time.

For an actual Euler path `v=gdot`, the exact observable balances are

    P_a = M_a Xadot,
    Padot = - integral_boundary(Da) p n,
    Sadot = - integral_boundary(Da) (x-X_a) cross p n,
    S_a = integral_Da (x-X_a) cross rho u.

The pressure terms on a shared material face have opposite signs in the
two adjacent parts of the full action. Only their complete sum cancels;
the torque on an individual nonspherical EPS tube need not vanish.

At the stationary base, invariance `u0.n=0` gives
`integral_Da u0 = integral_boundary(Da) x (u0.n) = 0`. Its centroid and
second mass moment are fixed. This is why the actual invariant tubes,
rather than frozen arbitrary material cubes, are the reference parcels.

## 2. A fully reconstructed positive translation block

There is an exact common translation/boost family, not merely a formal
mass appended to an orbit action. For any stationary Euler field `u0,p0`,

    u_V(t,x) = u0(x-X0-Vt) + V,
    p_V(t,x) = p0(x-X0-Vt)

solves Euler. The material flow is

    g_V(t,a) = X0+Vt+g0(t,a),

and every tube and the entire ambient phase are transported by this SAME
map. The extra constant velocity one-form is exact on R3, so it does not
change any closed material-loop circulation. This statement is not an
import of a freely variable harmonic circulation on a periodic torus.

For the stationary finite-energy-density law of 0098 with `E[u0]=0`,

    E[rho |u_V|^2/2] - E[rho |u0|^2/2] = rho |V|^2/2.

The density is the total fluid density, not the selected tube volume
fraction times that density. The common translation in the canonical
one-form is paired with the total momentum density `rho V`. No finite
whole-space energy or whole-space boost integral is asserted: the
identities are finite-window densities followed by the declared stationary
ensemble average. Translation covariance gives the same base energy at
`x-X0-Vt`.

For the compact internal columns of 0085, the induced velocity is compact
inside one tube and has zero integral there. Its cross energy with `V`
is therefore zero before ensemble averaging. Different tube columns have
zero kinetic cross because their actual induced velocities have disjoint
supports. The complete internal kinetic energy is retained, not appended
to a duplicate centroid fluctuation energy.

This proves the zeroth spatial-gradient coherent mass block. A slowly
varying boost is not a Galilean symmetry; its full affine and gradient
blocks require the joint variational calculation below.

## 3. What the compact pair proves at an affine reference jet

Write `v_i=xi_i cross omega0`, with both `xi_i` and `v_i` compact inside
the same invariant parcel. For any such divergence-free velocity,

    integral v_i = 0,
    integral (r_j v_i + r_i v_j) = 0.

For the pair constructed in 0085, the actual spin rows and KKS pairing are

    L(Q)=0,           L(S)=B n,           Omega(Q,S)=B.

Consequently, for an affine incompressible field `U=a+h r`, with
`beta=axl(skew h)`, the physical pairing is

    rho integral U.v_S = beta.L(S) = B (n.beta).

The lifted affine generator

    Gamma_U = U-(n.beta) Q

has `Omega(Gamma_U,S)=0`. Its centroid displacement is unchanged, since
the compact divergence-free Q has zero integral. Together the common
affine part and the physical relative angle implement rigid co-rotation
correctly. These are exact tangent/current identities, including the
ambient contribution (zero for these particular compact internal columns).

They do not yet identify an independent macro momentum on a fixed Kelvin
leaf. Nor does a zero KKS cross by itself prove that a prescribed finite
phase-space trajectory has material velocity equal to its Biot--Savart
velocity.

## 4. The reconstruction term which an exact joining retains

Let `C xi = curl(xi cross omega0)` and, on the stated full-space sector,
`V xi=P(xi cross omega0)`. A material displacement `eta` about the actual
stationary flow has Eulerian velocity

    v_material = eta_t + curl(eta cross u0).

The same fixed Kelvin datum also requires `curl v_material=C eta`, with
the prescribed periods/zero modes. Hence its reconstruction equation is

    eta_t = A eta,
    A eta = (lambda P-curl)(eta cross u0)
          = V eta-curl(V eta)/lambda.

This is a linear identity at each actual retained cell, including macro
wave number zero. It is not an extra all-wave-number nonlinear invariant
manifold demand. Compactness of `VQ,VS` makes the tag normal source zero,
but does not set the interior reconstruction residual to a relabeling.

For a retained orbit embedding `eta=E z`, the KKS/Hessian restriction gives
its own finite constrained action. If its solution is denoted `z(t)`, the
material lift associated with that SAME phase embedding has residual

    R = E zdot - A E z.

Its actual velocity would be `V E z+R`, rather than just `V E z`. Preserving
the old energy/current requires at least `C R=0` with the matching
circulation periods and an explicit reconstruction of that gauge. It
cannot follow merely from `Omega(E,R)=0`. If R is retained instead,
the kinetic, spin and material observable corrections must also be
retained. The complementary-variable construction in 0095 does this
without assigning an inertia by hand.

For the full coadjoint Hessian H and the genuine material Jacobi stiffness
K at a Beltrami base, 0091 establishes the particularly exposing identity

    H(xi,xi)-K(xi,xi) = rho ||A xi||^2.

Thus replacing the orbit closure by a restriction of actual material
configuration paths changes the quadratic form in a computed, generally
nonzero way. The difference is not a new elastic constant.

## 5. Reduction order and the next executable same-action repair

The right relabeling momentum of the full material action is the pullback
one-form `[g^*(v flat)]` modulo exact forms. Keeping the complete value
fixes the Kelvin circulations; reducing it yields the Euler orbit and
its reconstruction. Keeping only the subgroup which preserves selected
tags fixes only its restricted momentum map. That is a useful tagged
cotangent representation, but not permission to import the same full-leaf
H while introducing extra free momentum directions. Conversely a finite
configuration ansatz need not be preserved by the full right relabeling
group, so restriction and full Kelvin reduction need not commute.

For clarity, even a finite canonical action exposes the issue. Adding
phase constraints `chi(x,p)=0` with multipliers changes reconstruction to

    xdot = H_p + (partial_p chi)^T lambda.

The restricted Hamiltonian is still computed from H, but its derivative
need not be the original physical velocity. The reaction term is part of
the action; naming the construction a Dirac reduction does not erase it.

The immediate constructive repair is to retain the missing material/
Kelvin complement alongside the exact centroid/ambient split above and
eliminate it from the SAME action with its observable map. Attempt 0095
supplies the full operator pencil and frequency-dependent Schur formula.
The zero-frequency and slow-frequency expansion, where an actual inverse
and its bounds hold, then computes the corrected inertia, stiffness and
physical spin rather than copying the two-column values. Alternatively,
the material-first positive cages of 0091/0093/0094 provide a distinct
configuration closure whose own Kelvin reaction is kept explicitly.

The N3 affine Cauchy--Born exclusion of nonaffine relaxation is compatible
with a declared constrained energy or action. What remains to be specified
in a physical same-action theorem is whether the retained phase constraints
preserve the physical material reconstruction, or which computed reaction
and observation terms accompany them. No unrestricted microscopic
trajectory closure is silently added here, and neither representation is
declared a parent no-go.

## Source scope

The material split, boost, moment and reconstruction equations above are
derived directly here and in the cited campaign attempts. For the general
distinction between cotangent/semidirect-product reduction and its
momentum-map values, the primary reference is Marsden, Ratiu and Weinstein,
*Semidirect Products and Reduction in Mechanics* (1984),
https://www.cds.caltech.edu/~marsden/bib/1984/01-MaRaWe1984/MaRaWe1984.pdf.
Its abstract and reduction description were inspected; no infinite-fluid
analytic theorem or finite closure is imported from that abstract.
