# Exact induced-mean and cotangent repair of the retained angular moment

## 1. What the same force moment necessarily induces

For a modulated patch force F=xi cross omega and Fourier convention exp(i k.x),
the slowly coherent force is its Fourier transform at the macroscopic wave
vector. With zeroth moment zero, its leading value is -i k_j M_ij. The new
0070 moment is

    M_ij=c delta_ij-epsilon_ijm L_m/(2rho),

where L is the patch angular-momentum response, summed/averaged with the
same stationary intensity as the action. Leray removes the scalar-gradient
part. Consequently the induced physical mean velocity is

    v_s = curl L/(2rho)+higher-gradient terms.

The mean is not zero merely because every unmodulated compact patch has zero
integrated momentum. The zeroth patch moment controls uniform modulation;
the first patch moment controls its gradient. Reflection pairing does not
remove curl of an axial spin, which is a polar velocity.

More fundamentally, for every divergence-free macroscopic test displacement
U, the full Euler identities give

    Omega(U,xi)=rho <U,xi cross omega>
               =rho <U,P(xi cross omega)>.

Thus the retained affine symplectic source and the induced polar mean are
the SAME mixed pairing. This statement uses the full self-adjoint Leray
projector and the full spatial/ensemble inner product. It does not require
isolated-cell kinetic-energy factorization, an all-k invariant subsystem,
or a zero mean assumed from symmetry. At first gradient the pairing becomes

    rho <v_s,Udot>=<L,curl Udot/2>=<L,betadot>.

These equalities include periodic integration by parts or the recorded
boundary-current terms. They explain why the leading source cannot be
retained in the symplectic form while its identical mean velocity is omitted
from the kinetic energy.

## 2. Full cotangent action, before and after physical centering

Let P0 denote the independently supplied mean-momentum coordinate before
the reaction's induced mean has been added. Let H_s be the FULL internal
Euler quadratic energy, including the kinetic energy of that induced mean.
The physical total mean momentum is

    p=P0+rho v_s.

The actual kinetic-square decomposition is

    H_total=|P0|²/(2rho)+P0.v_s+H_s
           =|p|²/(2rho)+H_s-rho|v_s|²/2.

The cross P0.v_s is required by the same Euler kinetic density, not an
optional new constitutive coefficient. Omitting it double-counts the freedom
to specify the physical mean independently of the reaction mean.

The cotangent one-form of 0070 in its uncentered variables is

    Theta=<P0,dU>+<L,dPhi>.

Use the physical angle Phi=q+beta and the exact mixed-pairing identity above.
Changing to the PHYSICAL mean momentum p gives

    Theta=<p,dU>+<L,dq>.

Indeed the +L.dbeta part of dPhi equals +rho v_s.dU, precisely the shift
from P0 to p. This is a change of the full one-form, not just a Hamiltonian
square completion, and therefore preserves the symplectic structure.

Equivalently retain the old variables and eliminate P0 in the full first-order
action:

    L=P0.Udot-|P0|²/(2rho)-P0.v_s+L.Phidot-H_s,
    P0=rho(Udot-v_s).

The result is

    L=rho|Udot|²/2+L.qdot-[H_s-rho|v_s|²/2].

The apparent new beta-rate coupling cancels exactly at its retained order.
It has not been removed by changing the original physical mean displacement:
U still obeys p=rho Udot after the correct mean momentum is varied. This is
the one-action interpretation required by 0057 section 2, which already
states that gradient-induced polar means and their subtraction must remain.

The result is independent of an arbitrary normalization of the response
moment. Multiplying its affine angular moment changes BOTH v_s and the
cotangent cross by the same factor. Tuning that moment cannot evade the
identity. The cancellation is a route-class fact about this point-mean
Eulerian/material-ensemble split, not a no-go for a finite material parcel
with a distinct physical centroid/within-parcel moment decomposition.

## 3. Full reaction operator and the retained second-gradient changes

Write the leading induced mean as v_s=C s, with C of first spatial order,
and the angular response as L=D* s. Here

    C=(curl D*)/(2rho)

to this order. The original reaction energy is

    H_s=<s,P s>/2+<s,N q>+q.H_QQ q/2.

Physical mean centering changes the reaction operator to

    P_center=P-rho C* C.

This subtraction has second gradient order. It does not change the strictly
positive zeroth-order P0 or the leading positive internal action. In the
original slow-gradient neighborhood, positivity persists by the bounded
second-jet estimates; no positivity outside that neighborhood is claimed.
If coordinate-like internal directions also induce a mean E q, retain the
entire block Gram subtraction:

    P_center=P-rho C* C,
    N_center=N-rho C* E,
    H_center=H_QQ-rho E* E.

The eleven-orthogonal physical angle direction has no first-gradient mean;
its first allowed contribution therefore does not alter the leading mixed
rate cancellation proved above. Its higher-order contributions are kept in
this full block formula when computing the complete gradient coefficients.

Vary the time-reversed fluid reactions independently, with the same physical
U,q. Both the induced mean and KKS moment reverse sign before variation;
the mean subtraction and Hessian agree. Their reduced even action has

    J_center=D* P_center^-1 D,
    K_center=H_center-N_center* P_center^-1 N_center,

and kinetic term qdot.J_center qdot/2, NOT Phidot.J_center Phidot/2.
The full reaction-space inverse is retained. All second-gradient corrections
must be computed after the mean subtraction, not copied from 0070's uncentered
coefficient set.

For an exposing scalar isotropic constant-P example, ignoring other gradient
blocks only for this exact test, J0=j and C=curl D*/(2rho) imply

    J_center(k)=j+j² k²/(4rho)+O(k⁴).

Thus even this repair changes a retained gradient coefficient. In physical
(U,Phi) coordinates it gives b=-j/2, m_U=j/4 and m_Phi=j²/(4rho), so
`m_Phi-b²/rho=0`. With g=-kappa/2, the general 0066 normal form gives
`C_eff=C-kappa*j/(2rho)` for this example, rather than the uncentered
separable specialization. If N or the original gradient blocks are nonzero,
their complete Schur derivatives also contribute; this illustrative formula
is not substituted for that full calculation.

## 4. Physical verdict and continuation

The 0070 compact angular-momentum identity and response construction remain
correct: its sphere moment really is physical angular momentum of the
induced velocity, and its finite-carrier Hessian is positive. What fails is
the inference that this source can give an independent absolute-angle
inertia while the physical point mean is separately assigned rho|Udot|²/2
without the induced-mean kinetic/cotangent cross. Restoring that term gives
relative-rate inertia and restores zero leading centroid optical transfer
for the otherwise separable sector:

    g=-kappa/2, b=-j/2,
    l=g-kappa*b/j=0.

This is a specific joining defect repaired by the full centered action,
not a refutation of the original smooth Euler/Cosserat objective or a demand
for exact arbitrary-wavelength closure. The normal-form equivalence remains
available with the correctly recomputed gradient coefficients, but the
requested physical mean transfer remains an unsatisfied construction.

The failure-generated replacement is an actual FINITE MATERIAL parcel
centroid and its internal absolute angular momentum, with all within-parcel
affine velocities, shared-face pressure forces and ambient returns retained.
That coordinate is not the pointwise ensemble mean used here. Its exact
cotangent/kinetic split may retain an independent internal moment after
centering. The parent has activated this distinct route in 0075; it must
derive its split rather than relabel the mean corrected above. Compact
induced velocities can simplify its boundary flux but do not by themselves
evade the present point-mean identity.

Route verdict: the uncentered 0070 positive-coupling join is refuted by its
omitted mean kinetic cross; its corrected full action and operator jets
are established as stated. Evidence scope: exact slow-affine mean-cotangent
repair, with the parent physical coupling objective active.
