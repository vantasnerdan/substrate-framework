# P247 attempt 0006, gate F1 — derivation of the leading two-clock interaction law

## 1. Z2 theorem: the boost channel carries no classical static interaction

Measured this run (`f1-routing-z2.json`): the extended fixed-J functional
satisfies

    E[S, +chi] = E[S, -chi]

to machine precision, simultaneously on the committed symmetric root and
on randomized S backgrounds, for every angular channel probed. The chi
dependence of the functional enters only through even combinations
(s_m ~ (u^T eta xi)^2 with u = (cosh chi, sinh chi n_hat): cosh even,
spatial-u odd; all contractions quadratic in chi and grad chi). Hence:

  (T1) dE/dchi = 0 at chi = 0 for EVERY S — no configuration class can
       source the boost field at linear order (supersedes the parity
       argument of attempts 0005/0006-partial, which needed background
       symmetry).
  (T2) With the chi-quadratic operator positive (kinetic 1/16 > 0, mass
       Lambda^2 > 0), chi = 0 is the exact classical minimum for every
       background: static clocks do not interact through the boost
       channel AT ANY ORDER. Refutations: F_pair_exchange (as a
       classical channel) and route G's boost-tail mechanism.

The C1 mixed term (linear in d_mu chi on its own) must therefore pair
with an equal opposite linear term on every solution — consistent with
the measured exact zeros: the source is identically co-closed, j is a
curl-type current with no gradient content.

## 2. The interaction law lives in the massive S/response sector

With the boost channel closed, the leading classical two-clock energy is
the overlap interaction of the massive static fields. Write the static
energy of one field sector as E = integral [ (grad phi)^2/2 + m^2 phi^2/2
- V(phi) ] with vacuum phi_vac, phi_i a single-clock solution, and define
the effective source

    j_i(x) = (-nabla^2 + m^2) phi_i(x)   (per canonical normalization),

which is core-localized because the tail solves the linear equation.
The standard lump-lump expansion (kept to leading order in the tail
overlap; higher vertices are O(e^{-2 m d})) gives

    E_int(d) = integral integral j_1(x) G_m(|x - x' - d zhat|) j_2(x')
               dx dx',      G_m(r) = e^{-m r} / (4 pi r),

a Yukawa-type law with the field's own mass and multipole-dependent
power prefactors. Sign: positive (repulsive) when the source overlap is
positive.

Application to the committed root: the S rows are radial profiles
(f_1, f_2, f_3 channels of the hedgehog-like configuration). Each row
with its own fitted tail mass m_a contributes its own pairing term; the
slowest decay dominates the asymptotic law. Measurement plan (all
single-centered, no two-centered solver needed for the leading law):

  M1. Fit m_a from the large-r profile: ln|phi_a - phi_vac| vs r.
  M2. Compute j_a = (-d^2/dr^2 - (2/r) d/dr + m_a^2) phi_a on the grid.
  M3. Evaluate E_int(d) over a separation ladder by the double quadrature
      with the analytic G_m; fit the exponent against m_a.
  M4. Attribute the residual box-growth E(R) ~ +0.38 R (attempt 0003
      ladder) to sectors by evaluating the extended functional's
      components at chi = 0 across R.

Caveat recorded: if the profile fit shows a power-law (massless) tail
rather than an exponential, the same formula holds with m -> 0 and the
pairing becomes Coulombic (~1/d), which would also attribute the linear
box-growth. M1 decides which regime the model is in.

## 3. What remains for F2

The full-functional moderate-separation check (two-centered solve) is a
separate build; the in-run validation is the internal consistency of the
law (single-exponent decay, source core-localization, refinement
stability) plus the C-M5S-008 confined-sign comparison.
