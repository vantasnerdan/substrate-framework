# Independent payload and materialization correction review

Reviewer: Herdr geometry-review pane `w3:p4`, fresh non-author.  This is a
bounded post-0269 materialization check.  It does not reopen the scientific
review of C-CST-018.  The reviewed payload has SHA-256
`d91c1e42c1880ea2fa53e2efab45fd8cc4fd382737e668b39dca3411d9afd9b7`;
the pinned 0265 claim statement has SHA-256
`f94d4b0e1d83fe545e65fa69a319a067ce6f3e07f2835dda1e1532dba245a046`.

## Verdict

`route_verdict: established as stated for the bounded payload,
materialization, and terminal-closure correction`

The payload faithfully materializes the statement accepted by the 0269
scientific review.  The repaired reusable API now has the literal cylindrical
Euler residual and the correct independent energy/moment normalizations.  Its
direct and consumer tests expose the original `E`-versus-`Q` false green.
The proposed terminal graph is also an honest objective-level reconciliation:
it closes the new issue200 construction only through its four actual parent
obligations and does not turn the historical N1--N7 route statements into
established claims.  I find no remaining scientific or implementation
correction inside this bounded review.

## Claim text and metadata

Removing the Markdown title and only the final pending-status sentence from
`0265/claim-statement.md` gives the YAML block statement exactly, modulo the
serialization's terminal newline.  No quantifier, coefficient, source-cost,
current term, boundary exclusion, or prepared finite-window limitation has
changed.

The metadata has the right scope:

- `C-CST-018` is an additive `compatible_extension`, with `supersedes: []`.
  Thus it neither rewrites C-CST-017's periodic noncontractible-tube result nor
  claims that the old straight-tube route succeeded.
- Dependencies `C-CST-009` and `C-CST-017` are used only for their accepted
  canonical action/current and branch-substitution algebra.  The assumptions
  say this explicitly and forbid importing C-CST-017's noncompact periodic
  microscopic source as the compact source.
- The compact geometry, response, normalization, and current chain is carried
  as direct evidence through 0263/0268, 0265--0267, 0262/0264, and the 0269
  joint review.  The 0269 review itself hash-pins and checks the load-bearing
  0250 two-polarization and 0260 vector-lift inputs, so omission of those two
  source paths from this already closed evidence list is not a dependency
  defect.  Adding them would be provenance convenience only.
- `symbolic_verified` is honest because the assumptions identify the analytic
  construction and independent reviews as the existence proof and restrict
  symbolic checks to exposing finite algebra, metric, centering, inertia, and
  sign errors.  The prospective `accepted_in: v0.183.0` and adjudication path
  belong to the pending promotion transaction, not to present scientific
  authority.

## Literal reusable boundary

The reviewed module and direct test file have SHA-256 hashes
`f17c251b370e53d6dc2cd512ff0ea9eb0610f246e88bde8475eb5cf49bcf705b`
and
`bdf2c7e0449c28b6118b54c7cc821aa04f965b961055ce2c745607e329636e52`.

For supplied symbolic fields on `r>0`, `compact_ring_fields` returns

    u_r=-psi_z/r,       u_phi=I/r,       u_z=psi_r/r,

and all three components of `(u dot grad)u+grad p`, including
`-u_phi^2/r` in the radial row and `u_r u_phi/r` in the azimuthal row.
Pressure is explicitly specific pressure.  The routine neither derives a
Bernoulli label nor claims finite-radius existence.  Its coordinate checks
reject an explicit nonpositive radius and coincident coordinates, while
unknown symbolic signs remain documented caller hypotheses.  That is the
correct conditional cylindrical API boundary.

The repaired normalizer returns

    j = 2 rho f Q/(3 V_cell),
    a = E/(3 V_cell),
    ell_tag^2 = Q/V_tag,

where `E=int |u|^2 dx` and `Q=int chi |x-X|^2 dx`.  Units, known-sign checks,
and the immutability of the returned matrix agree with the payload.  In
particular, varying `E` leaves `j` unchanged and affects only `a`; this is the
load-bearing repair of the earlier copied-answer test.

The captured affected replay has SHA-256
`595865a1f82687a9d67fa77f3427d716d58da9493b73b97c3efcfae9cf95e05f`
and reports `21 passed in 7.79s`, with exit zero.  The direct probes are
appropriately exposing:

- solid rotation checks centrifugal pressure and fails when that pressure is
  removed;
- a nonzero-poloidal example produces the exact azimuthal residual `-2`, so
  both cylindrical metric contributions matter;
- an independently integrated unit ball gives axial inertia
  `8 pi/15=(2/3)(4 pi/5)`, verifying the `2/3` factor without copying the API;
- the resulting measured `j`, `a`, and `ell_tag^2` feed the canonical joint
  symbol, reproduce its mass block, retain a cubic first mismatch, and satisfy
  the canonical phase-generator identity.

These checks validate the reusable algebra and its immediate consumer.  They
do not independently certify the compact Euler existence or prepared response,
which remain licensed by their reviewed analytic evidence.

## Terminal closure mapping

I additionally reviewed the now-concrete reconciliation at these frozen
materialization surfaces:

- `proposal.yaml`, SHA-256
  `31b1e43fab99606bb0d86fe3c55309aaa9a276d8eead6a65da03d3edce70417f`;
- `0271/historical-route-state.yaml`, SHA-256
  `e19d7c77e8cf9c85ecb541ed86fa43393c9cad028a42a426418140f989afc087`;
- `0271/closure-map.md`, SHA-256
  `b7ba2aef78907ef72f8d5bb2e66f57bd9229156ef0afa8d68f003e500373d840`.

The operational P200 graph in `proposal.yaml` is a sound replacement for the
issue200 execution ladder, with the more precise supplier map recorded in
`closure-map.md`:

| Parent obligation | Actual closure | Status represented |
|---|---|---|
| `P200-G` geometry/density | 0263 finite-`R` carrier on the 0261/0256 inputs, independently reviewed in 0268; 0265 disjoint assembly/action chart and 0258 twist | established for one fixed compact Euclidean ring ensemble, pre-limit tag law, and positive measured `j` |
| `P200-R` same-field histories | 0250 field-changing response, 0260 Kelvin lift, 0257 whole-law angle, and 0266 centered full-oblique Schur/cellwise phase, reviewed in 0267/0269 | established for the full finite-`K` inverse, explicit `O(|K|^-2)` acoustic cost, pressure, and finite-window histories |
| `P200-A` action/current | 0262 full-form supplier/current bridge, reviewed in 0264, then the 0265 observed-state normalizer, connection, and current, jointly reviewed in 0269 | established for both inherited forms, all mixed crosses, the same measured `j`, and literal endpoint/`q_t` memory |
| `P200-P` promotion | 0269 individual claim review, repaired 0270 reusable algebra, this 0271 payload/reuse check, then registry/release/generated-memory agreement and terminal validation | **active** until promotion and full-suite receipts exist |

The archived historical graph remains a different route record.  In
particular, N2's regularized Biot--Savart bending/twist formulas and N3's
resulting explicit filament moduli remain `active` in
`historical-route-state.yaml`; their licenses remain unearned.  The new compact
prepared-continuum proof uses different microscopic and response suppliers.
Issue200's explicit operational supersession therefore changes which route
answers the unchanged common compact-continuum objective without asserting
that N2, N3, or the rest of the old route-specific ladder has been solved or
exhausted.  Keeping `supersedes: []` in the scientific claim payload is
consistent with that distinction.

The preservation file contains the complete prior top-level question,
objective, N1--N7 graph and licenses, route frontier, and candidate history.
It retains N1's historical established status, N2--N7 as active route nodes,
and the old Biot--Savart/modulus licenses as unearned.  In particular,
`alpha=L_v*T/6` remains unsupported; the reviewed compact preparation instead
realizes `alpha=j*nu^2/4`.  This is route replacement at the issue200 parent
level, not retroactive success or an exhaustion certificate for unused
historical candidates.

The ten-achievement map is also internally honest at this checkpoint.  It
pins the object, dependency closure, framework fit, exposing oracles,
candidate history, independent reviews, and reusable code to their concrete
artifacts.  For replay it records only the completed 21-test affected set and
defers the branch-selected full-suite result to `campaign-validation.stdout`.
Registry, release, generated canon/memory agreement, and the final debt/status
receipt remain work of `P200-P`.  Therefore the map records reviewed
scientific closure plus a promotion execution plan, but **not** a terminal
campaign assertion at the present state.

The many earlier `parent_state: active` fields elsewhere in `proposal.yaml`
are dated append-only observations of their respective attempts.  They need
not be rewritten or reinterpreted as current competing parent graphs.  They
remain historical until the forthcoming final record supplies the single
authoritative operational supersession; neither their presence nor the three
earned P200 scientific licenses permits marking the campaign terminal before
`P200-P` earns `LP200-4`.

Remaining achievement: execute the central C-CST-018 registry/release,
generated-record, full-suite, and terminal validation transaction represented
by `P200-P`, then use the distinct non-author operational merge path.  No
further scientific dependency is exposed by this correction check.
