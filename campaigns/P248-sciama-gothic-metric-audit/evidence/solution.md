# P248 corrected construction

## What survives

The gothic relaxed formulation remains a useful exact chart of Einstein
dynamics once its density weight, harmonic coordinate condition, total source,
and pseudotensor convention are kept explicit.  The issue paper's optical
block formulas also remain exact, but only on the determinant-slaved slice
they actually parameterize.  The dimensional Planck-cutoff numbers and the
conditional retarded recurrence pass their respective checks.

## Minimal complete metric map

Use coordinates `x0=c*t`, choose a flow orientation `s` once, and introduce an
independent positive lapse `N`:

```text
v^i = V^i/c
g_00 = -N^2 + gamma_ij v^i v^j
g_0i = s gamma_ij v^j
g_ij = gamma_ij
```

For positive-definite `gamma`, this has

```text
det(g) = -N^2 det(gamma)
g^00 = -1/N^2
g^0i = s v^i/N^2
g^ij = gamma^-1^ij - v^i v^j/N^2.
```

Conversely, any metric with a positive-definite spatial block and negative
temporal Schur complement yields

```text
gamma = g_spatial
v = s gamma^-1 g_mixed
N^2 = -(g_00 - g_mixed^T gamma^-1 g_mixed).
```

The ten-component Jacobian is `-2*s*N*det(gamma)`, so this is a genuine local
bijection.  The paper is recovered by imposing
`N=(det gamma)^(-1/6)`, which is a constraint, not an inverse theorem.

## Correct harmonic/material relation

Let `nbar=(det gamma)^(1/3)` and
`M=partial_t nbar+div(nbar V)`.  The exact gothic time component is

```text
H_s = partial_t(nbar^2) - s div(nbar^2 V)
    = 2 nbar M
      - 2 nbar (1+s) V.grad(nbar)
      - nbar^2 (2+s) div(V).
```

For the material-flow convention `s=-1`, simultaneous material and gothic
continuity is exactly material continuity plus `div(V)=0`.  It does not require
`nbar` to be spatially constant.  For the paper's printed `s=+1`, the joint
condition instead contains `4 V.grad(nbar)+3 nbar div(V)=0`; switching between
these equations without switching the metric sign is the source error.

## Weak conformal and shear repairs

At leading weak order the paper's potential map gives
`|grad Phi|^2=(c^4/4)|grad log(nbar)|^2`.  Its positive gradient energy is
therefore positive Newtonian field energy, not the negative energy it earlier
requires.  The energy ledger needs the opposite sign; a common material action
must still be supplied before the spatial stress is called constitutive.

For positive principal optical indices, replace the additive expression
`n_i-nbar` by

```text
lbar = (log n_1 + log n_2 + log n_3)/3
E_i = log n_i - lbar.
```

Then `sum E_i=0` and `product exp(E_i)=1` exactly, with no conformal double
counting.

## Action-level completion

For any covariant action pulled back as `S[g(q),fields]`, variational chain
rule gives `E_q=J^T E_g`.  The nonzero Jacobian makes `E_q=0` equivalent to
`E_g=0`.  Applying this to accepted `C-STG-001` and supplying the same metric to
accepted `C-WLN-001` and `C-WLN-002` gives a complete conditional optical
Einstein-scalar theory together with the massive square-root/Hamiltonian and
massless null/affine-geodesic worldline sectors.

This solution deliberately does not rename that field redefinition a
microscopic derivation.  A substrate action, universal matter realization,
spin-two emergence, physical value of `G`, preferred foliation, and empirical
equivalence remain excluded and therefore create no hidden debt inside the
promoted theorem.
