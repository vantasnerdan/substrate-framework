# Source and authority audit

## Independently reviewed inputs

- P253/0068 with P253/0071 supplies the transported real charge tag, joint
  Euler--Maxwell current/action/constraints and local smooth flow.  It does
  not supply a charged carrier or a chiral current.
- P253/0077 with P253/0081 supplies the exact traveling Lorenz equations,
  `Phi=phi-c A_z`, the toroidal Ampere primitive, modified
  Grad--Shafranov/Bernoulli identities and regular-band tag stabilizer.  The
  final 0081 review and verdict hashes are
  `09f8580b2ff4120ef85046b22d85875ea1371e6b2a028913fe888a28cbcb00ce`
  and
  `27ae74404344242b78b1757d7ffb00af307248dfbbc87b6464289cd81f95e099`.
- P253/0080 with P253/0084 supplies the local charged Cao branch only in its
  declared finite uniformly subluminal window.  The final 0084 review and
  verdict hashes are
  `3fb7d2a69489668f6168d7f57eb9b4e86e35621064dd5f55172ca24270279aac`
  and
  `7793688a23156d1b015ed82505f5e48b88223e19e8a36563c2f57ef63f0a74c7`.
- P253/0097 with P253/0099 supplies conditional two-state propagation,
  coherence and shared-action kinematics while leaving a physical Cao flavor
  current blocked.  The final 0099 review and verdict hashes are
  `ae3d5495baca1e3fbf6c50a56903b13f9fdcde9712117a4297431e2b703e58b7`
  and
  `d9102986bb05ccb150aaf0de50a8b4e9fe41ab4849abddc66ad1538362775b75`.
- P253/0100 with P253/0101 supplies the forced Ertel current, exact
  forced-lock condition, regular Cao zero density, magnetization
  superpotential and punctured azimuthal equilibrium.  The final 0101 review
  and verdict hashes are
  `1ed8266d2a6c8487fbdaf5b629fd2d22ac77232c9912c84feeb71b7c0b569ebb`
  and
  `1a9b89788bbe5458b134445e708de660125cdaf6408922513246151562f0cbe3`.
- P253/0085 with P253/0086 supplies the strong-saddle and Maxwell-continuum
  boundary.  It prevents the equilibrium identities below from being used as
  a persistence theorem.

These are active-proposal inputs, not accepted release authority.  The
accepted release remains `v0.183.0`; no P253 result is promoted here.

## Exact calculations derived in this attempt

The following rows are derived rather than attributed to a source:

- the two-label quotient law
  `chi_c D_t lambda=(curl f) dot grad Theta` and its exact continuity flux;
- the automatic variable equilibrium lock from the reviewed modified
  Grad--Shafranov first integral;
- the support counterexample to a global constant lock for the existing
  band-supported tag;
- the order-`g^2` contourwise zero-mean compatibility row;
- the signed-charge-factored Lorenz source
  `chi_c/epsilon_EM*(1-c^2/c_EM^2-c W_z/c_EM^2)`;
- the exact Ampere-eliminated electric/magnetic lock response;
- the radial `m=1` Volterra Green identity;
- the global straight-column exponential-profile cancellation and its explicit
  Liouville realization solving both the Poisson and exponential profile laws;
  and
- the independent inner/outer source-moment control of the annular
  homogeneous dipole coefficients.

The first scratch sign argument for the straight-column response was rejected
during author audit.  Positive Cao vorticity fixes the physical sign of
`P_s/(R s)` oppositely to the scratch `Omega>0` convention, so the electric
and magnetic terms can oppose.  No universal sign or Cao fixed-profile
refutation is retained.

## Domain and interpretation boundary

The exponential law is a global straight-column representation with
noncompact vorticity support.  Its two-dimensional electric and fluid
potentials have logarithmic affine tails; only the globally matched derivative
identity is used.  It is not a compact Cao free-boundary profile.
The annular determinant is a source-level range theorem and does not prove
that its two source variations arise from one admissible Cao material leaf.
That lift must carry `delta chi_0=delta(zeta_0/lambda_0)`, the induced
`delta P`, free boundary, circulation, mean-radius and Casimir rows.

The phase `Theta=theta` removes the entire noncompact symmetry axis and has a
bare gradient norm singular at the axis and nonintegrable without a declared
outer treatment.  A compact ring defect requires a new circle map and does
not inherit `q_A=zeta`.  Nothing here constructs its moving domain, finite
action, Lorentz chiral representation, flavor interaction, P2, P4 or P6.
