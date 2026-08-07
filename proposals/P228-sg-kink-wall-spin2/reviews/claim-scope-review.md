# P228 claim-scope review

## Review surface

This review checks the four proposed statements against the implementation,
tests, primary verifier, direct algebraic rederivation, accepted dependencies,
and the failure reasons recorded on PR #6. It is an author-side scope audit,
not independent maintainer acceptance.

## C-EGR-009

Verdict: ready for maintainer review as a conditional microscopic-model claim.
The six wall variables are explicit, local transverse gradients have a derived
kink normalization, the tensor map is invertible, and repeating an axis breaks
the oracle. The one-wall-per-`ell` icosahedral cell remains a declared model
choice and is named as such.

## C-EGR-010

Verdict: ready for maintainer review only with its assumption list intact. The
polarization count is the nullity of a rank-four constraint matrix, the kinetic
form is positive, the Fierz--Pauli ratios are the nullspace of a separate
rank-three gauge solve, and both steps have load-bearing mutations. Constant
volume, force balance, and collective relabeling redundancy are model premises;
the code does not claim the unconstrained scalar sine-Gordon action proves them.

## C-EGR-011

Verdict: ready for maintainer review as a low-energy conditional completion.
The coefficient `B` is upstream of the gauge solve, `G` exists only in the
gauge-completion ledger, and common coefficient rescaling changes it. Deser's
theorem is cited only after its locality, universal coupling, conservation, and
self-coupling premises are stated. No observed `G`, bare inverse coupling, or
regulator enters the API.

## C-EGR-012

Verdict: ready for maintainer review at the line-source ceiling. The static
result now has a second-moment error enclosure instead of a point-source claim
at one profile length. The radiation calculation reuses the landed exact-source
full-retardation transform. Both consumers obtain `G` from the same completion
ledger. This is not a prediction for a three-dimensionally compact body.

## Prior-review regression

P228 does not recreate the three blocking shortcuts from PR #6: it has six
geometric wall separations rather than 24 phases chosen to diagonalize a target
matrix; mode count and parameter count are rank-derived and mutation-sensitive;
and the static line-source approximation carries an explicit radius-dependent
error bound. Scientific promotion and issue closure remain maintainer decisions.
