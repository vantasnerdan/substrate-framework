# Actual joint history residual and the missing hybrid observation

Use the same physical Fourier sign as euler_observation: coefficients of
exp(-i K.x), so curl of a spatial amplitude is C=i[K cross]. The material
moment comparison below is the Fourier transform of the actual momentum
measure. Constant rho is the complete fluid density.

## 1. Evaluate the complete point-to-centroid difference

For one actual positive material tag let X be its centroid, r=x-X,
P=integral rho chi v, A_ij=integral rho chi v_i r_j, and
B_ijl=integral rho chi v_i r_j r_l. Its resolved momentum minus its
collapsed centroid momentum is exactly

    exp(-i K.X) integral rho chi v [exp(-i K.r)-1].       (1)

At K=epsilon k the first two Taylor coefficients are

    d1=-i A k = i k cross S/2-i I_dot k/2,
    d2=-B:(k tensor k)/2-i(k.X)d1.                      (2)

Here S=integral rho chi r cross v and I=integral rho chi r tensor r;
I_dot=A+A^T follows from material transport and centering. Velocities in
B are absolute: replacing them by v-X_dot drops X_dot tensor I. Even a
uniform Euler boost of a symmetric pair has d1=0 and nonzero d2. The
last term in d2 is the moving centroid Fourier phase. Sum these rows
and retain the full continuous ambient: it is unchanged by collapse,
so it cancels in this *difference*, not in the fluid mass or action.

Taylor's real-phase remainder inside the centroid phase is bounded by
|epsilon|^3 integral rho chi |v| |k.r|^3/6. A linear response or time
derivative needs its corresponding differentiated material-moment bound.
This is a moment identity; it does not assert that point particles prove
a smooth-Euler existence or Sobolev estimate. The executable API computes
(2) from positive finite material data, and the tests differentiate the
literal defining Fourier sum (1), including its time derivative.

If actual perturbed moments have their own K jets, *their* first
variations and coefficients enter (2). In particular an acoustic
first-order shape-rate response contributes at total order two, as
does the zero-order second momentum moment. A leading optical G0 or
an acoustic axial spin is insufficient to evaluate that row.

Define the actual physical hybrid velocity by J_H/rho and its
displacement by its independently specified initial value plus its
time integral. If X_pt denotes the similarly initialized point-fluid
observable, then

    U_H,tt=X_pt,tt - rho^(-1) partial_t Delta J,         (3)

where Delta J is the complete point-minus-hybrid difference, including
(2), reference phase and the material variation. Thus the required
second-order preparation correction is the actual derivative of the
whole Delta J coefficient, not merely an assigned spin inertia. The
initial U_H and U_H,t are separate rows.

## 2. Exact evaluated branch-to-equation defect

For a transverse helicity C=tau with tau=+k or -k set
beta=j/(2rho), alpha=j nu^2/4, d=1+j k^2/(4rho), and

    T=[[1,-beta tau],[tau/2,1]],
    D=diag(a k^2,nu^2+cT k^2), M=diag(rho,j).

The second-order physical stiffness is the canonical micropolar matrix

    K=[[ (rho a+alpha)k^2, -2alpha tau ],
       [ -2alpha tau, j nu^2+j(cT-alpha/rho)k^2 ]].    (4)

For actual branch amplitudes z and actual physical observation
Y=Tz+e, define r=z_tt+Dz. Direct substitution gives the residual

    M Y_tt+K Y = M T r+(K T-M T D)z+M e_tt+K e.        (5)

The exact mismatch in parentheses is

    [[0, beta(rho cT-rho a-alpha) tau k^2],
     [j(cT-a-alpha/rho) tau k^2/2, 0]].               (6)

Thus it starts at cubic spatial order, with no hidden second-order
mismatch. Longitudinal Phi is separate and has exact residual
j[B_L,tt+(nu^2+cL k^2)B_L]. The finite positive rho,j,nu make the
coefficient weights nonzero; to earn an o(k^2) dynamical statement,
the actual physical observation and its first two time derivatives,
and the branch defects, need the corresponding o(k^2) bounds. An
initial action match alone does not supply them.

Equations(3),(5) identify the concrete acoustic source row for repair:
the point mean of C016 is controlled, but the full hybrid acceleration
contains the derivative of its shape and second-moment correction.
The optical spin source0234 and normalizer0228 do not erase (3).

## 3. Keep the full current and its initial charge

For the actual0235 representative Q_ij=q(t)epsilon_ijk U_H,t,k,
set S_int=S_full-div Q and mu_int=mu_full-Q_t. The total angular
balance residual is exactly unchanged:

    S_int,t-div mu_int+ax sigma
       =S_full,t-div mu_full+ax sigma.                 (7)

The accumulated spin current is defined with its own initial value:

    G_full-G_int=Delta G(0)+div integral_0^t Q(s)ds,
    integral_0^t q U_H,t
       =q(t)U_H(t)-q(0)U_H(0)-integral_0^t q_t U_H.    (8)

This accumulated current is distinguished from a literal displacement
dipole whose Euler/Lin identity can include the extra advective terms
in0223 equation(4). Its use does not set those terms to zero.
Matching the source history in (5), the actual measured spin, and the
literal force/torque divergence of0232 is still required for the joint
constitutive supplier. A superpotential identity is not that history.

## Dependency, error and current table

| Row | Source and actual object | Error/cost retained | Effect on join |
| --- | --- | --- | --- |
| Point acoustic D,V | C-CST-015/016, fixed C016 field | finite control norm, cubic Euler remainder, fixed-time diagonal | positive a and point acceleration |
| Acoustic angle | 0231 covariance Gram on same field | finite smooth first-jet norm; odd K parity | Phi_ac=C A/2 at leading order |
| Full hybrid acoustic mean | (1)-(3), actual moved tags and ambient | full shape-rate, second-moment and centroid-phase derivatives | explicit extra acceleration target |
| Optical angle/full spin | 0234 positive tag fractions and Kelvin fields | leading clock o(k^2), correction epsilon_N k^2, full pressure costs | review0242 assesses the actual optical supplier |
| Complete phase/energy | 0228, reviewed0230 | source cross jets, intrinsic jets, polynomial-cost diagonal | action normalization after new controls |
| Physical spin/flux | 0232/0235 | initial G, q_t memory, ambient and moved boundaries | complete current residual (7), not a constitutive assumption |
| Periodic tube | 0236 on C016 | actual finite core and twist at stated periodic scope | same-cell geometry; distinct from Euclidean rings |
| Joint dynamics | (5)-(6) | actual r,e and two time derivatives | remaining source correction is explicit |

Route verdict: established exact joint residual and complete moment interface.
Evidence scope: exact finite algebra and material calculus, not completion of
the actual supplier. The next executable route is both parity passive controls
for (3), preserving the elliptic tag at first order and matching the full
new action with the already-reviewed normalizer.
