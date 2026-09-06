# One bounded 0053 target correction

The independent `0053` review found that `0048` had earned its exact Cao cell
recurrences and finite spectral-window asymptotics but had not constructed the
fixed-domain full Euler graph/Riesz theorem. Root applied the requested
claim-domain correction only to the four claim-bearing `0048` artifacts. The
verifier and its receipts were unchanged and were not rerun.

The correction:

- replaces the resonance location by
  `l*L_Phi/(pi*c_l*delta*L)+O(1/(delta*L^2)+1)`;
- fixes the KKS sign convention by declaring `A e=-i sigma e`, retains the
  positive Krein sign and square-root frequency magnitude, and leaves exact
  density/Fourier/real-complex normalization open;
- makes the moving-boundary full Green/Leray graph jet conditional on `HJ2`;
- retains the `m=0` coupling-window and Schur-tail algebra while naming the
  missing all-`m`, axis, collar, interface and exterior-pressure complement;
- makes the avoided-crossing matrix conditional on that graph realization;
  and
- scopes equations (52)--(53) to a conditional formal fixed-order residual,
  rather than a constructed Euler branch expansion.

| Artifact | Before SHA-256 | After SHA-256 |
|---|---|---|
| `0048/derivation.md` | `3b35b7cb56089a4fb4e78fff535f5f72128fc39c9d2cba9d1299f2d0a32560fd` | `381c59e5debe1886847acd296c795192605576cdf55afc721688d2aea9783b25` |
| `0048/result.yaml` | `13569f910e968012815e68be56cdeb8eff6f5286e03a525e6055ce9692687b29` | `e6b09d0f5e01d80c03aeba2987f0406695ca83f7abf546524d16e837d1cf894b` |
| `0048/source-audit.md` | `d465aab71ef242dd62597086f8401816b4a7ac40eda6e4799d73cb09d23e2577` | `4db69ddded6fd2d84e749271e2f13efb293c00645a8c8efe311cb35820971b03` |
| `0048/validation.md` | `3107ddb505a0cf16273d182176e186115b4e9c9ed21ce727d192bf46c3db2202` | `9ac3c21c514e4024a48f57b0fcc30599f64316aec037008357a87b15ab6490b5` |

The stale-claim scan finds no remaining assertion that `0048` established the
graph-domain Riesz family, exact physical KKS normalization, or a fixed-order
Euler branch. `git diff --check` passes. Active `0052` remains separate and is
not evidence for this review check.
