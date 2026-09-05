# Actual-observation pullback and the full current implication

This continuation consumes the explicit local-block hypotheses of
`exact-gram.md`, the actual finite-window joint histories of 0250, and the
universal material balance/current identities of 0246. It does not use
the special C016 source formulas as a construction on a new field.

## 1. Normalize in the actual measured initial variables

Let E_B map the finite baseline initial amplitudes to the literal initial
physical state (U,Phi,U_t,Phi_t). Let E_C be the same map on the auxiliary
columns, including their generally nonzero pressure/observation tails.
Let H_can and Omega_can be the canonical energy and phase forms on this
physical state, with the actual positive rho,j and the chosen prepared
moduli. The total measured state map after adding coefficients L Y is

    E=E_B+E_C L Y.                                      (1)

The correct normalization equations are

    H_phys(B+C L Y)=E* H_can E,
    Omega_phys(B+C L Y)=E* Omega_can E.                  (2)

Both sides are quadratic in the same actual source coefficients. Subtract
the observation pullback from each full physical form on the combined
baseline/control space. In particular replace the auxiliary form by

    H_c-E_C*H_can E_C,
    Omega_c-E_C*Omega_can E_C,                          (3)

and the baseline and cross blocks by the corresponding differences.
Equations (2) become precisely the affine simultaneous Gram equations in
`exact-gram.md`, with zero target. Invariant-gap smoothing makes the
corrections (3) arbitrarily small in the retained finite parameter norms.
The explicit right inverse there therefore solves (2) exactly. This
retains the actual initial observation map rather than replacing it by
its limiting map before normalizing the action.

The baseline finite-window state inverse is bounded on each chosen
inventory. Select the auxiliary observed tail smaller than its inverse
margin. Then E is invertible on the retained physical state, and (2)
gives the canonical forms in **those actual measured variables**. Costs
of this finite inverse enter L and the auxiliary carrier choice. The
normalizer is still not a proof of the history equations: those are
separately supplied by 0250 and preserved by the finite-window smoothing
estimate.

## 2. Whole-law symmetry and the measured histories

Use the same positive whole-field orientation/reflection law, material
tags, and fraction labels as the physical supplier. Rotate the auxiliary
field, patch, source, observer, and input representation together. The
finite construction is on these actual sources, and its conjugation
symmetry at opposite wave numbers is retained. A common compact set of
rotated source patches has the same gap and derivative bounds.

There is no change to material density or to probabilities. Positive and
negative **energy signatures of perturbations** are not negative material
fractions. The reflected law continues to cancel the forbidden even
polar/axial blocks. Its allowed odd curl/current blocks remain those
derived in 0241/0246; they are not canceled by calling the law isotropic.

After choosing the finite physical histories and the exact form targets,
choose the auxiliary carrier to make all required observed time derivatives,
initial-map errors, and full spin/current rows smaller than the remaining
accuracy budget. This choice comes after the baseline carrier and its
finite source costs. It leaves the actual measured equations

    (U,Phi)=T(A,B)+e,
    T=[[I,-j C/(2rho)],[C/2,I]],  C=iK cross,
    A_tt+a|K|^2 A=r_A,
    B_tt+[nu^2 I+cT|K|^2 P_T+cL|K|^2 P_L]B=r_B,         (4)

with the same o(|K|^2) finite-window bounds, including e_tt and the
separately controlled leading optical clock error. Longitudinal covariance
is the actual 0257 common-law observable, not a missing transverse chart
component. All full material G and S rows required by (4) are retained.

## 3. The universal balance calculation transfers unchanged

Set mu=rho a, alpha=j nu^2/4,
gammaT=j(cT-alpha/rho), gammaL=j cL. Their positivity comes from the
explicit chosen prepared targets. Substitution in the canonical coupled
equations gives exactly the algebraic identity from 0246:

    M0(U,Phi)_tt+K2(U,Phi)
      =M0 T(r_A,r_B)+(K2 T-M0 T D)(A,B)+M0 e_tt+K2 e.   (5)

The middle matrix is cubic in K. This identity depends on the physical
gain map and coefficients, not the C016 microscopic geometry. The actual
state inverse turns it into an operator estimate for all retained initial
amplitudes.

For the same sources the literal material equations are

    J_H,t=div F_full,
    S_full,t=div N_full-ax F_full.                       (6)

Here F_full and N_full include transport and pressure; the symmetric
convective tensor does not change ax F_full. Use the actual whole-law
current improvement Q_ij=q(t)epsilon_ijk U_t,k from the material identities:

    S_int=S_full-div Q,
    N_int=N_full-Q_t,  J_int=J_H, F_int=F_full.          (7)

The full measured spin control, rather than (7) alone, gives
S_int=j Phi_t+o(|K|^2). Lower-endpoint constants and q_t memory stay in
the actual current primitive and boundary action. Auxiliary sources are
included in the same finite spin/current error list.

Subtracting the canonical variational equations (5) from (6)-(7) gives

    P_T div(F_int-F_can)=o(|K|^2),
    div(N_int-N_can)-ax(F_int-F_can)=o(|K|^2).            (8)

For transverse virtual displacement v and arbitrary virtual rotation psi,
the actual bulk difference therefore obeys

    integral [Delta F:grad v+Delta N:grad psi
                                +ax Delta F dot psi]
      =integral [-div Delta F dot v
                       +(-div Delta N+ax Delta F)dot psi]
      =o(|K|^2).                                        (9)

This is the same complete bulk constitutive representative proved in 0246.
It retains the incompressibility pressure gradient and the boundary terms
on a cut domain. It does not assert equality of freely prescribed literal
surface tractions. Parameter-dependent initial preparations do not change
the exact local equations (6); no derivative of their chosen coefficient
matrix is silently substituted for a material current.

## Scope of the join

Equations (1)-(3) establish exact normalization in the actual initial
observation variables under the normalizer's block/gap hypotheses.
Equations (4)-(9) transfer the full history/action/current implication under
the actual physical supplier hypotheses. This discharges the algebraic
current join without assuming that energy matching creates the histories.

An actual periodic compact carrier with the local constant-curl annulus,
twist, separated source patch, and positive measured geometry is still
required to instantiate all these hypotheses on the requested one field.
The exact fixed 0211 ring is the response prototype; its isolated geometry
does not establish that periodic compact assembly. A small C^4 deformation
licenses transfer of nonzero leading gain determinants on a compact band,
but the actual curved cohomological/full-pressure source construction and
its cost bounds must be instantiated on the resulting field. This file
does not replace those source obligations with continuity of a matrix.

`route_verdict: established as stated for the actual-observation form
normalization and conditional full material current implication`

`evidence_scope: exact pullback and balance identities using the named
physical supplier hypotheses; one-field geometry/response instantiation
and independent review remain open`
