# Same-field histories imply a complete periodic current/action bridge

This proof has explicit source inputs: the actual prepared acoustic supplier
of0243 (independent review0244 pending), the optical supplier0234 at its0242
reviewed scope, the0228 complete form normalizer at its0230 reviewed scope,
and0232/0235's exact physical momentum, torque and current representative.
It derives the new join; it does not deem a pending source review complete.
The periodic geometry supplier is separately reviewed in0245. The compact
Euclidean geometry obligation is retained unchanged.

## 1. One physical observation, one preparation, one error diagonal

Use one stationary C016 field and the positive whole-field O(3)/TR law,
including0234's two tag fractions. Whole-state rotations move the field,
lattice, source, observation and common input representation together. Let
A be transverse acoustic branch amplitude and B an unrestricted axial optical
amplitude. All sources are finite linear initial-data maps from their common
configuration/rate amplitudes. Their sum is therefore a single actual linear
Euler/Lin preparation on each realization. Cross terms belong to its complete
second-variation action and are retained in0228, not lost by summing energies.

The physical hybrid U is obtained from its actual J_H/rho, with independently
specified initial displacement; Phi is the fixed covariance Gram observation.
The acoustic supplier repairs the full point-to-hybrid coefficient, not just
the point mean. Its first angle is C A/2. The optical supplier has leading
hybrid U=-j C B/(2rho), with actual positive measured j, and its prescribed
second-order Phi response. Inversion symmetry removes an optical-to-polar
even coefficient and an acoustic-to-axial even coefficient. Whole O(3)
covariance makes the transverse U condition compatible with both branches;
it does not remove longitudinal spin.

In the common ordered prepared sequence these statements give

    (U,Phi)=T(A,B)+e,
    T=[[I,-j C/(2rho)],[C/2,I]],
    A_tt+a|K|^2 A=r_A,
    B_tt+[nu^2 I+cT|K|^2 P_T+cL|K|^2 P_L]B=r_B.          (1)

Here e and its first two time derivatives, r_A and r_B, are o(|K|^2)
in the corresponding physical operator norms on each fixed time interval.
The leading optical clock error is separately made o(|K|^2); it is not
confused with a second-coefficient error multiplied by |K|^2.

Choose both finite acoustic controls before the leading optical width h.
Include their full source norms in the finite constants of0234's reviewed
ordering and0228's normalizer. Apply the same ordering to the finite list
of actual moment and current derivatives, initial observation inverse and
all cross forms. The extra source norm can be arbitrarily large at finite
accuracy; it is fixed before h and does not need a polynomial bound in the
acoustic band's low frequency. This is the single ordering stated in0243.

## 2. Substitute the actual histories in the canonical equations

Let M0=diag(rho I_T,j I) and define

    mu=rho a, alpha=j nu^2/4,
    gammaT=j(cT-alpha/rho), gammaL=j cL.                 (2)

The source targets choose cT>alpha/rho and cL>0 before solving actual output
controls. This gives positive bulk spin curvatures. It is an explicit
prepared-law choice, not an inferred unprepared spectral stiffness. Let K2
be the canonical incompressible micropolar second-order stiffness with these
coefficients, including both off-diagonal blocks -2alpha C.

Equation(5) of0241, independently evaluated against the canonical API, gives

    M0 (U,Phi)_tt+K2(U,Phi)
       =M0 T(r_A,r_B)+(K2 T-M0 T D)(A,B)
                           +M0 e_tt+K2 e.              (3)

The second term is exactly cubic in K. Thus the actual physical histories
satisfy the canonical coupled equations to o(|K|^2) on the common fixed-time
sequence. This conclusion uses the actual source histories in (1); full
initial phase or energy matching on its own would not imply (3).

For completeness, the branch state map from independent initial branch
positions/rates to (A,B,A_t,B_t) has uniformly bounded inverse on each fixed
time interval at K=0: acoustic solutions are 1,t and optical solutions are
cos(nu t),sin(nu t)/nu, each with unit Wronskian. T has determinant
(1+j|K|^2/(4rho))^2 on the five physical positions, so it and its inverse
are bounded in the long-wave neighborhood. The controlled C^1 errors keep
the actual observed state map invertible there. Consequently (3) is an
operator statement across all retained initial amplitudes, not an equality
verified along one selected solution.

## 3. Retain literal full transport, pressure and angular reaction

The exact0232 balances use the actual collapsed-tag plus continuous-ambient
momentum J_H, full intrinsic material spin S_full, and their transport fluxes
T_H,C_H:

    J_H,t=div(sigma_H-T_H),
    S_full,t=div(mu_H-C_H)-ax(sigma_H).                 (4)

T_H is symmetric: every parcel term is M V tensor V and the ambient term
is rho chi0 u tensor u. Its linear variation stays symmetric. Define the
complete force/couple flux representatives F_full=sigma_H-T_H and
N_full=mu_H-C_H. Then ax F_full=ax sigma_H; no convective term was erased.
The force/torque sources and moved boundaries are exactly0232 equations(7)-(11).

With the same actual U define Q_ij=q(t)epsilon_ijk U_t,k. Use the actual
whole-law q of0235, evaluated after the law and tag weights, and set

    S_int=S_full-div Q,       N_int=N_full-Q_t,
    J_int=J_H,               F_int=F_full.              (5)

The acoustic first-spin identity and the explicit optical full-spin target
in0234 give S_int=j Phi_t+o(|K|^2), including the chosen time-derivative
norms. This is where the actual optical controller is used: (5) alone is
not a detector identity. All initial integrated-current rows retain their
independent lower endpoint and q_t memory as in0241 equation(8). A literal
displacement dipole with Euler/Lin advective terms is not silently replaced
by the primitive of S_full.

The exact balance after (5) is therefore

    rho U_tt=div F_int,
    j Phi_tt=div N_int-ax F_int+o(|K|^2).              (6)

The same-density momentum remains unchanged by the exact boundary action
in0235, including its U_tt and gradient-U_t momentum cancellation.

## 4. Compare actual current divergences to the variational representative

Take any canonical isotropic gradient-energy representative satisfying
c_s+c_a=gammaT and 2(c_s+c_tr)=gammaL. Its density is

    W=mu||sym grad U||^2+alpha|curl U-2Phi|^2/2
              +c_tr(div Phi)^2+c_s||sym grad Phi||^2+c_a||skew grad Phi||^2.

The canonical force stress F_can=partial W/partial(grad U), couple stress
N_can=partial W/partial(grad Phi), and local derivative satisfy
partial W/partial Phi=ax F_can. Combining (3) and (6) now evaluates the
literal current difference:

    P_T div(F_int-F_can)=o(|K|^2),
    div(N_int-N_can)-ax(F_int-F_can)=o(|K|^2).           (7)

A longitudinal force gradient is the incompressibility pressure multiplier.
Unlike conservation alone, (7) uses the measured histories, their complete
coupled residual, and the full measured spin normalization. By the bounded
observed-state inverse above, it holds in the coefficient operator norm on
the retained physical state, not only for one prepared amplitude.

For transverse compact/periodic virtual displacement v and arbitrary virtual
rotation psi, integrate the actual stress/couple difference by parts:

    integral [Delta F:grad v+Delta N:grad psi+ax Delta F . psi]
      =integral [-div Delta F . v
                    +(-div Delta N+ax Delta F).psi].    (8)

It is o(|K|^2) by (7). This is a complete periodic bulk constitutive
representative at the claimed order. It is not an equality of literal
tractions on a free surface: (8) retains the corresponding surface terms
on a cut domain, and0232/0235 supply their exact localization/superpotential
bookkeeping. The known one-dimensional curvature representative freedom
is precisely the divergence described in0232, not an extra bulk modulus.

## 5. Inherited action and scope

The same actual sources, including both acoustic parity controls and optical
current targets, enter0228's full phase/Jacobi-energy correction. Its output
errors are retained in (1), while its prescribed forms are T* M0 T and the
corresponding full branch stiffness through degree2. Transforming BOTH forms
by the same physical T returns M0 and K2, as independently reviewed in0230.
Together with the actual state histories (3), this gives the inherited
prepared continuum action, the physical equations and the complete bulk
current equivalence (7)-(8) in the same physical variables. Time-dependent
current memory remains in (5) and its boundary action, not in a silently
redefined autonomous material coefficient.

This is an actual prepared linear-response/second-variation statement on
compact time windows, conditional on the exact source/error inputs explicitly
listed at its opening. It does not assert an unrestricted Euler invariant
manifold, acoustic-time uniformity, finite-amplitude nonlinear stability,
arbitrary prescribed traction, or the unconstructed stationary compact
Euclidean ring ensemble. The periodic geometric supplier has its separate
actual density and tube statement; replacing it by a Euclidean knot still
requires that distinct construction.

Current route status: proof assembled, pending the new acoustic supplier
review0244 and a review of this new joint implication. No accepted claim or
terminal campaign completion is inferred before those source licenses hold.
