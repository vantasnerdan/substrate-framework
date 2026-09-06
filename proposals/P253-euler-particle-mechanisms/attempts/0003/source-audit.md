# Slobodeanu source audit for P253/0003

## Pinned source

- Radu Slobodeanu, *Steady Euler flows and the Faddeev--Skyrme model with
  mass term*, arXiv:1405.3469v3 (16 August 2019),
  <https://arxiv.org/abs/1405.3469>.
- The PDF retrieved from <https://arxiv.org/pdf/1405.3469> on 2026-09-06 has
  SHA-256
  `fdc453ac26b0568e5b1344940c537aa547829f8542d3a0da015ff226aa1a81bf`.
- This audit uses the corrected v3 source, especially equations (3), (5)--(7),
  (10)--(13), Proposition 2, the paragraph following Proposition 2, and
  sections 4.2--4.4. No statement below is inferred from the abstract alone.

## Exact source theorem and conventions

Let `(M^3,g)` be an oriented Riemannian three-manifold with volume form
`nu_g`, let `(N^2,h,J,omega)` be a compact Riemann surface with its area
two-form, and let `phi:M->N` be `C^2`. For a nonnegative `C^1` potential
`P:N->R`, the paper's strong-coupling static functional is

```
E_sigma2,P[phi]
  = (1/2) integral_M ( |phi^* omega|_g^2 + 2 P(phi) ) nu_g.       (S1)
```

Its independent compactly supported map variations give

```
J dphi [ (delta phi^*omega)^sharp ] + grad_h P(phi) = 0.          (S2)
```

Proposition 2 defines

```
B = phi^*omega,
i_V nu_g = B,
V^flat = star_g B.                                                (S3)
```

For a critical `phi`, `V` is a **steady** incompressible Euler field and its
Bernoulli function is the basic function `P_E=P o phi`:

```
div_g V = 0,
i_V d(V^flat) = -d(P o phi),                                     (S4)
```

equivalently `V x curl V = grad(P o phi)` in the paper's Euclidean
convention. Pointwise `|B|=|V|`, so (S1) is

```
E_sigma2,P[phi]
  = (1/2) integral_M |V|^2 nu_g + integral_M P(phi) nu_g.         (S5)
```

The first term is the unit-density kinetic energy. The second is a
solution-dependent Bernoulli-potential term in the static map problem; it is
not the incompressibility multiplier in Hamilton's principle for Euler.

The converse direction has narrower scope. Near a point where a steady Euler
field `V` is nonzero, one takes the local quotient along its streamlines.
`star V^flat` is basic, descends to an area form on a local surface, and the
target metric is then chosen compatible with that descended form. A global
field `phi` follows only for the paper's `S`-integrable flows, whose possibly
singular one-dimensional streamline foliation is simple. Thus the forward map
is global when its field data are global, while the converse is only local
without a separate global quotient theorem.

## Domain transfer to an isolated Euler object

For an isolated Euclidean Faddeev--Skyrme configuration the source imposes
`phi(x)->phi_infinity` as `|x|->infinity`, compactifying physical space to a
based `S^3`. Finite map energy then gives finite kinetic energy through (S5),
provided the potential is integrable. These conditions do not themselves
construct a freely translating localized Euler carrier or a nonlinear Euler
stability neighborhood.

The paper's Proposition 4 is specifically an `R^3`, finite-energy,
axisymmetric, no-swirl, `C^2` nonexistence result for potentials
`(1-phi_3^a)^b`, using an external theorem with decay and asymptotically
constant pressure. It is not a no-go for swirl, non-axisymmetric fields,
other potentials, compact domains, or all Euler particle candidates. The
explicit examples following it use `R^2 x S^1` or `S^3`; their transfer to an
isolated Euclidean particle is therefore a separate construction.

## Full-model distinction

Equation (13) adds the quadratic Dirichlet term `kappa |dphi|^2`. Proposition
3 maps that **static** full Faddeev--Skyrme critical point to a **forced**
steady Euler equation, with

```
F = kappa ( div C_phi - (1/2) grad |dphi|^2 ),   F perpendicular to V. (S6)
```

Consequently the unforced incompressible Euler correspondence is to the
strong-coupling quartic-plus-potential static functional, not automatically
to the standard full dynamical Faddeev--Skyrme theory.

## Topological normalization in the source

For a based map `phi:S^3->S^2`, with standard target area form and
`d alpha=phi^*omega`, the source uses

```
Q(phi) = (1/(16 pi^2)) integral_(S^3) alpha wedge phi^*omega.     (S7)
```

Under (S3), the paper's flux-helicity convention is

```
H_flux(V) = (1/pi^2) integral alpha wedge i_V nu_g = 16 Q(phi).  (S8)
```

This is a genuine exact topology/flux dictionary at the stated global
normalization. It must not be confused with the usual Euler kinetic helicity
`integral V^flat wedge dV^flat`. Rescaling the target area form or the
physical velocity normalization rescales the dimensional flux helicity, so
the integer `Q` alone does not select an absolute physical action unit.

## Source-scoped verdict

The source establishes a global field-to-steady-flow map, a local converse,
static stress balance, energy-density equality for the quartic term, and the
global Hopf/flux-helicity dictionary when the domain and quotient hypotheses
hold. It contains no time-dependent Euler action equivalence, no transfer of
admissible variations or stability, no physical quantization rule, and no
Lorentzian propagation theorem. Those are tested rather than presumed in
`derivation.md`.
