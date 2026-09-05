# Independent review: complete noncommuting Schur jets

Reviewer: `/root/smooth_core_review`, distinct from the author of 0062
and its API addition. Date: 2026-09-05. One bounded exact-algebra/code
review under AGENTS.md and the physics skill. Earlier orbit functions
and the dependent 0061 regression are outside this transaction.

## Decision

Established as stated: `hermitian_schur_jet` computes every coefficient
through k² of the inverse momentum block and the complete Hermitian
Schur complement, with noncommuting order and complex adjoints preserved.
No load-bearing correction is requested.

The inputs are coefficient jets, not a derivation of differentiability,
microscopic locality, positivity of the final action or independence of
fluid cells. A corresponding bounded-operator use requires norm-C² jets
and a bounded inverse at zero (for a self-adjoint momentum operator,
a uniform coercive lower bound suffices). Merely positive quadratic
forms on an infinite-dimensional space would not supply that inverse.
The README appropriately leaves the operator license to its separate
full-fluid construction; the API itself implements finite matrices.

## Exact ordered algebra and implementation

I read the frozen README, new dataclass/helper, its new test and the saved
ten-test passing receipt. Exact coefficient matching and independent
matrix differentiation are the strongest oracles; the existing receipt
was inspected and reused.

Write P=P0+kP1+k²P2+o(k²), with no factorial in P2, and I0=P0^-1.
Multiplying P by its inverse coefficientwise gives

```
I1=-I0 P1 I0,
I2=I0 P1 I0 P1 I0-I0 P2 I0.
```

These products solve both inverse identities in their specified order;
there is no commutativity premise. In particular the positive repeated-P1
term cannot be omitted merely because one is interested in second order.

For Q=H-N^* I N, the code subtracts every ordered term
`N_left^* I_middle N_right` with left+middle+right equal to the desired
order. At order two this includes the six distinct triples (2,0,0),
(1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2), with the correct left
adjoint. This includes every mixed N1 contribution and the full inverse
response; none is replaced by a diagonal or single-cell approximation.
Since k is real and the P,H jets are Hermitian, the resulting reduced
coefficients are Hermitian as they must be.

The new test constructs a positive P0 and genuinely noncommuting complex
Hermitian higher jets, then obtains the inverse and Schur complement by
direct matrix inversion. Its derivatives at k=0, divided by factorials,
match every returned coefficient. Thus the test checks the coefficient
convention against an independent rational-matrix calculation, not by
copying the API's ordered sums. The omitted-P1-square mutation is exposed.

The helper enforces three coefficients, nonempty compatible dimensions,
known finite entries, Hermitian P and H, and decidable positive-definite
P0 via its leading principal minors. Undecidable symbolic positivity
remains a documented caller hypothesis. The returned matrices and tuples
are immutable. No claim that a positive P0 alone ensures positive reduced
H is made.

## Provenance and disposition

- `README.md`: `1377b6db6dab3f9ddeb4dae0bdb815c344e474334e8f02bec16eed8f620c9ec0`
- `src/substrate_framework/euler_orbit.py`: `89001e02762a3bcfd79805e8e8ebff999ea0d5e8101e70a2846278ce1069a39b`
- `tests/test_euler_orbit.py`: `87081c678a5c06795b5fb5bee82c9808f51c5bebc1edac72a317dc2544f18af1`

Acceptance is recommended as conditional importable exact coefficient
algebra. There is no numerical small-ratio remainder, no Gaussian inverse
factorization claim and no parent microscopic-closure conclusion.
