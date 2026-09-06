# P253/0072: finite-energy degree-minus-two Euler tails and Coulomb-power interaction

## Frozen parent and route boundary

This is a failure-derived bare-Euler continuation of the independently
reviewed 0069 escape route. The parent P253 objective remains an enduring
same-substrate particle mechanism, with electron structure before neutrino
structure. This attempt asks one narrower but load-bearing question: classify
the smooth-core, homogeneous-degree-minus-two divergence-free velocity tails
that can produce a `d^(-1)` kinetic cross energy, and determine exactly which
SO(3) representations, signs and extra structure can or cannot turn that
power law into an orientation-independent signed scalar interaction.

The attempt does not infer charge, a persistent carrier, P5, an electron, a
neutrino or parent completion from a power law. The noncompact-vorticity route
is compared with the compact-carrier and constrained U(1) routes; it does not
replace them.

## Fixed input and conventions

The reviewed input is the 0069 result at its final corrected scope. In
particular,

    u_a(x)=(a cross x)/(1+|x|^2)^(3/2)

is smooth, divergence free, finite energy and in every finite `H^s`; its
velocity is `O(r^(-2))`, its vorticity is generically non-`L^1`, and two
translated copies have an anisotropic `d^(-1)` kinetic cross energy. That is
evidence for this candidate, not authority for the classification below.

Use physical kinetic energy

    E(u)=rho_0/2 integral_R3 |u|^2 dx

and Fourier convention `fhat(k)=integral exp(-i k dot x)f(x)dx`. A repulsive
like-sign coefficient means a positive `C/d` cross energy and force
`-grad_d(C/|d|)`. Every zero-mode, distributional source and boundary term is
retained explicitly.

## Route A: exact vector-spherical-harmonic classification

Freeze the asymptotic class

    u(r,n)=r^(-2)[V_r(n)n+V_T(n)]+o(r^(-2))

with an explicit weighted `H^s` remainder and a smooth divergence-free core
completion. Decompose the leading field into radial, poloidal and toroidal
vector spherical harmonics, including its `O(3)` parity. Derive
`div_S V_T=0`, the flux row `integral_S2 V_r=0`, the smooth-core extension
condition and finite-energy domain. In particular, decide rather than assume:

- whether the only SO(3)-trivial polar degree-minus-two tail
  `q n/r^2` is a distributional point source and is
  excluded by global smooth incompressibility and zero sphere flux;
- the smallest allowed SO(3) representations after that row is removed;
- how the example `u_a` sits in the `l=1` toroidal sector; and
- which weighted regularity makes the homogeneous tail plus a smooth core an
  admissible local Euler datum.

Transform the classified tail at low Fourier frequency and derive the exact
translation cross bilinear for every harmonic, with uniform control of the
weighted remainder and smooth-core overlap. A `d^(-1)` statement must include
its angular tensor, sign, remainder hypotheses and the condition under which
its coefficient is nonzero. Electric sign is earned only when the coefficient
paired with `q_1 q_2` is positive for every separation direction. A scalar
Coulomb channel requires more than absence of an `l=0` vector: its leading law
must be independent of every allowed internal or spin orientation. Apply
Schur/representation analysis to the bilinear; a nontrivial irrep normally
retains relative-state invariants unless a background lock, independent
channels or a derived average removes them.

## Route B: explicit isotropy escapes

If the one-field classification has no nonzero scalar channel, execute rather
than merely list the principal escapes:

1. an internal vector/director retained as a physical carrier observable;
2. a deterministic triad or multiple sector construction, checking whether
   linear superposition collapses it back to one vector and whether genuinely
   dynamically orthogonal/invariant sectors are being added;
3. a common-frame-locked composite, with its locking energy and same-Euler
   status stated;
4. a time or invariant-measure average, with the averaging dynamics, time
   scale and error needed to obtain an isotropic coefficient; and
5. a structured background that supplies a director, with broken SO(3)
   recorded rather than hidden.

Each route receives one verdict. A statistical isotropic mean is kept
distinct from a deterministic scalar interaction. No route may import the
three independent transverse fields of 0067 while calling the result a
single bare-Euler field.

## Route C: Euler evolution, localization and carrier consumption

For at least one smooth-core representative, derive the local Euler evolution
of the leading `r^(-2)` coefficient from the full pressure projection. Test
the positive invariant suggested by power counting: `u dot grad u=O(r^-5)`
and the pressure monopole cancellations may give `u_t=o(r^-2)`, making `V`
locally conserved. Derive the pressure moments and surface terms rather than
assuming this conclusion, and state the weighted Sobolev/moment hypotheses
and local time interval. Audit energy, absolute momentum, angular momentum,
impulse, helicity and KKS finiteness; finite energy together with non-`L^1`
vorticity can invalidate the campaign's existing moment/action joins.

Then test what two translated, globally overlapping tails mean for a
two-carrier state: exact initial energy splitting, pressure cross terms,
center/impulse observables, and the absence or presence of a same-family
restoring theorem. A `C/d` energy becomes a force law only after a persistent
translated family, moduli variation from the same Euler action, and controlled
deformation/radiation errors are constructed. Test opposite-tail neutral
cancellation as its own candidate.

Compare the result with the campaign's actual particle needs:

- compact or quantitatively localized energy/vorticity;
- a persistent same-family carrier and restoring mechanism;
- reciprocal source and force from one Euler action;
- an orientation-independent signed scalar observable; and
- compatibility with the retained-state/history construction.

A finite-energy algebraic tail may remain a valid particle candidate even
when compact vorticity is lost, but the cost must be explicit. Conversely,
failure of this tail family does not refute compact, topological, background
or foundation-extension routes.

## Verification contract

The exact oracle must derive the spherical divergence and flux identities,
the low-frequency transverse Fourier tensor, the inverse transform producing
the `d^(-1)` kernel, and at least one sensitivity mutation that changes the
SO(3) or sign conclusion. Analytic proofs carry the functional-domain,
remainder and Euler-evolution statements. Production numerics are outside this
attempt unless a separately preregistered scalar remainder survives the exact
analysis.

## Success and continuation

This attempt succeeds by delivering the complete classification and route
verdicts above. If it finds a deterministic scalar repulsive channel, the next
achievement is to place it on a persistent carrier and prove its same-field
interaction and restoring control. If only oriented channels survive, the
next achievement is the strongest explicit isotropy escape or a proof that
the candidate's orientation is a measurable particle degree of freedom. If
localization or persistence fails, preserve the exact Coulomb-power atom and
continue the compact/gauge candidates; do not convert a route verdict into a
global no-go.

Any time average must be autonomous and supply a scale-separated error and
time window. Any background or locked-triad construction must name its added
scale and broken or restored symmetry. The elliptic Euler pressure remains
instantaneous, so this route does not create a strict causal cone.
