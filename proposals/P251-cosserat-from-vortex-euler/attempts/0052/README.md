# 0052 — physical material momentum and Euler-orbit spin joining

Parent: P251 / issue #198, same conditional affine smooth-Euler objective.
Owner: Codex `/root/construction_review`; this directory only.

Frozen positive object: identify the common-angular momentum of 0048 with
the physical material-cell angular balance of 0051, retaining its surface
term, and derive translation mass/cross terms from the same Euler action
before joining them to the reduced spin sector. The reusable four-coordinate
Schur/ensemble algebra of 0049 is an input, not an independently supplied
kinetic term. Physical cell centers remain mass centroids unless a derived
and invertible field/current redefinition explicitly replaces them.

Candidate A: mean-center the full material Euler action before Routh
reduction; track the exact angular-current boundary improvement and its
pressure flux. Candidate B: extend the relative Euler orbit by the
Euclidean/Galilean momentum sector and compare its moment map to the material
centroid sector, retaining every boundary and mean-flow term. Register any
failure-generated repair append-only. Selection is one Euler action, exact
center/spin observable, no postulated mass, correct force/couple reaction,
and agreement with 0051/CST4. Equal time-reversal and spatial-parity partners
may cancel only explicitly computed odd cross terms, with independent fluid
reaction coordinates as in 0049.

Analytic oracle: exact material-map kinetic decomposition, Reynolds
transport, Noether moment maps, integration by parts including moving-cell
boundaries, and symbolic matrix variations. No numerical remainder,
empirical comparator, or fitted coefficient. A missing matching identity is
a construction task, not permission to identify two inequivalent spins.

Status: active; contract frozen before the joining derivation.

The exact material/action and boundary-current joining is recorded in
`material-joining.md`. It also derives the finite cell-mean response Gram,
its fixed-response projection, and the precise centered common-momentum
residual, rather than silently reusing the uncentered orbit matrix.

`verify.py` calls the shared 0049 Euler-orbit API and checks the exact
centroid, current-improvement, Gram projection, and paired surface-response
identities. Execute with `env PYTHONPATH=src` in this worktree so the newly
added API is imported rather than the older editable installation.
`first-run.txt` preserves that initial import-path failure; `stdout.txt`
records the repaired invocation. No equation or tolerance changed.

The material/current bridge is established as stated. The full selected-EPS
centered rotor remains active at the explicit response/KKS reconstruction
described in the derivation; this attempt makes no full parent-goal claim.
