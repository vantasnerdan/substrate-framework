# P253/0061 joined independent review of 0055 and 0060

## Frozen transaction and correction provenance

This transaction independently reviews `0055` and `0060`. The reviewer
authored or implemented neither target nor its APIs, tests, verifiers, repairs,
or receipts. The frozen README has SHA-256
`94d1a9c7dedf1080bbb07f6ad03b9870ba6eaf2bfaf5f7bf3684c8c1afff2d5a`.
Central activation exited `0`; the command, stdout, empty stderr, and exit
hashes are respectively
`0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The substantive pass requested one joined bounded correction. Its receipt has
SHA-256
`3ad8b9f174c9bb779daed8f8f7ea79ba4226272a214dcba2e1b88d249c002bbf`.
This final check was limited to the changed statements and their direct
scientific dependency edges. No unchanged oracle was rerun and no second
substantive pass was performed.

All receipt post-hashes match the worktree:

| Artifact | Corrected SHA-256 |
|---|---|
| `0055/README.md` | `bf46ea2c9c82d16e26cede0a699447ba79a2e69bdafedf5d1aa38fe6b7cb358a` |
| `0055/construction.md` | `5351f8bf3abd7b622c0902419f6c99200e5f142fd48e7e59e492625df41dca96` |
| `0055/result.yaml` | `938ad472fe8dad612877276da562285a2b6a9d3f5223ccee906a4ab47e3c81ef` |
| `0055/source-audit.md` | `db20e8e969613341d46db73bdc360d92b8d4d6044275586d97bc73d7720117c9` |
| `0055/validation.md` | `745335c8577b32c07319278c4caeae17614ba13301d9173825df31d0054ace08` |
| `0060/README.md` | `44ed68373f5307247bb3b79fd1a72646fd84fdc09a95b760f31ecef9cfd63dc3` |
| `0060/construction.md` | `487891b43cf40b810f6d547909adfdc37a23c64c6a2e868d1b007885ca12ec6e` |
| `0060/result.yaml` | `0eac38f436c384248e2c69fd31a9d8c0a41e892b6fa6d6d1cb02fc6401de1edb` |
| `0060/source-audit.md` | `91afbf5c4a6159cc8dc2cece066377a7c1fa193466de6c64c575ce77c289c06b` |
| `0060/validation.md` | `bef8cde5d00cd77989d2bff8b28802f78154387dd24404d1e05627e85751215d` |
| `0060/multipole-order-correction.md` | `2d72ea85ffa3623cc45adbe57a4e4b6bff0bea8fd33af480c6491f4e6fde3640` |

## Unit A: 0055 scale, pressure, and lattice theorem

The active Euler similarity has the claimed three-dimensional weights. Direct
change of variables gives circulation `A/B`, energy `rho_0 A^2/B^5`, helicity
`A^2/B^4`, and physical angular-momentum/KKS action
`rho_0 A/B^4`, with discrete topology unchanged under positive dilation. The
displayed fixed-label subfamilies retain continuous action freedom. The
bounded repair correctly excludes the action `J` itself from this inference:
prescribing `J` fixes it tautologically. Thus the result is a selection
boundary for topology and conserved non-action labels in bare fixed-density
whole-space Euler, not a theorem about boundaries, backgrounds, constitutive
scales, or every scale-breaking Euler environment.

The compact-swirl calculation uses
`Delta pi=-partial_i u_j partial_j u_i` and the crossed Cartesian contraction.
Its convergent trace-free quadrupole produces a direction-dependent
`pi=O(r^-3)` term and hence `-grad pi=O(r^-4)`, including the checked sign and
physical-density conversion. This is an instantaneous elliptic pressure tail,
not a Coulomb law or an all-time particle interaction.

Multiplying a classical hyperbolic action by an independent positive `K`
rescales its canonical form and action periods without changing its equations
or characteristic speed `c`; changing `c` does not fix `K`. This is an exact
conditional extension argument, not an Euler-derived physical foundation.
Finally, three-dimensional shell comparison gives absolute convergence for
`p>3` and first-spatial-moment convergence for `p>4`.

**Unit A verdict: established at the corrected bare-Euler and conditional
extension scope.**

## Unit B: 0060 stationary cell and neutral tangent

Disjoint compact velocities with flat collars and matched exterior kinematic
pressure constants sum to an exact smooth finite-energy stationary Euler
field. For one stationary component,

    integral u_i u_k = -delta_ik integral pi
      = -(delta_ik/rho_0) integral (p-p_infinity).

The repair therefore states the physical/kinematic convention correctly.
Exterior pressure is constant, so all static exterior pressure multipoles
vanish; the stationary assembly supplies no nonzero intercell pressure tail.

The proper-octahedral tensor identities are exact, but the six signed-axis
sites form `O/C4`, not a free 24-site orbit. The correction now recognizes
that a generic nonzero seed can have zero `C4` projection. A nonzero six-site
trivial-sector tangent consequently remains conditional on an exhibited
nonzero `C4`-invariant compact DA seed (or on one of the stated free-orbit or
overlapping-stabilizer alternatives).

For compact curl vorticity, the zeroth moment vanishes, the first moment is
antisymmetric, and the second moment is symmetric in its two spatial indices.
The corrected translated formulas retain every center term before averaging.
The proper-octahedral rank-three invariant is proportional to the alternating
tensor, so symmetry of the two spatial indices excludes it. Conditional on a
nonzero trivial-sector tangent, this cancels the Hodge dipole and quadrupole
and gives the valid upper bounds `v=O(r^-5)` and `grad v=O(r^-6)`.

Compact generators supported componentwise inside the flat collars, with
circulation, harmonic, centering, and stabilizer rows imposed separately,
give an exact direct-product coadjoint tangent and a split KKS form. Global
Hodge velocity, pressure, and the kinetic Hessian retain cross terms; no
invariant product dynamics follows.

An `r^-5` scalar tail in three dimensions has finite zeroth and first lattice
moments. It does not by itself yield a `C^1` Bloch operator. That continuation
still requires normalized invariant modes, a uniform full-pressure/Leray
`D(A_0)->X` coupling of the same order, graph-resolvent isolation, and
complement leakage control.

**Unit B verdict: established for the stationary composite, exterior-pressure
flatness, translated tensor projector, conditional far-tail implication, and
componentwise product-leaf geometry; blocked for a nonzero six-site invariant
tangent and for the Bloch band by the named constructions.**

## Correction check and joined verdict

The single correction fully closes the review findings. It preserves the
strong exact results while removing the two unsupported promotions: universal
selection beyond bare Euler in 0055, and automatic nonvanishing/operator-band
closure in 0060. The receipt's stale-language scan, YAML parse, and scoped
diff check are consistent with the inspected final text. No further bounded
correction is needed.

The joined result is therefore established as a classical scale/neutral-cell
boundary: bare Euler admits the exact similarity and instantaneous pressure
tail, while the octahedral construction supplies an exact stationary geometry
and a conditional `r^-5` neutral tangent. It does not supply an invariant
positive mode, Bloch band, strict cone, action selection, particle, P4, or
parent completion.

The next construction is a nonzero `C4`-invariant compact DA tangent (or a
free 24-site replacement), followed by an actual positive internal mode and a
uniform full Euler/Leray graph-domain coupling and complement theorem. The
parent campaign remains active.
