# Polish report (beacon 0118)

## Hypothesis under test

The bordered stall (0117) is under-iteration or a shallow artifact
dislodgeable by standard continuations. Killed below: four independent
probes return the same signature (monotone ray to just-above-cur,
frozen active set, detS O(10–10³)).

## Probe ledger

| Probe | Prediction if artifact | Result |
|---|---|---|
| 64-halving deep search | accept below 2⁻³⁰ | asymptote above cur → TRUE ascent |
| Warm reg-chain 1e-3→1e-4 | floor ∝ reg moves | zero motion → floor ≠ smoothing |
| p-chain 2→6 cold per rung | easy-p win carries | p3 3.3e-2 best; p≥4 stall; NO carry |
| Fine mesh 80×40 | floor drops | same stalls 0.07–0.21 → NOT discretization |

## Standing state (unchanged best)

Coarse bordered REG rung: res 2.0e-2, kap=1.018, rbar=0.996,
iz=3.23 (≈π post-hoc), umax≈1.31, μ=0.44, c=0.066. Feed numbers stand
(λ_ω=4.72 exploratory). c<0 ACROSS all runs — flagged: rows may select
off-physical branch; physical-branch selection is part of the next rung.

## Specified next rung (not attempted here)

1. Nested secant/Broyden on (μ,c) using FD row-Jacobian (measured O(10)
   sensitivities, smooth), PDE-Newton inner — decouples focusing
   stiffness from parameter sensitivity.
2. Re-audit J/B/Crow AT the stall state (all audits at bump so far).
3. Resolve c-sign: constrain c>0 in the border or prove the branch.
4. Then: fine mesh → Maxwell O(g) → production feed → G-a2 DONE.

G-a2 remains BLOCKED (polish rung open). No verdict altered.
