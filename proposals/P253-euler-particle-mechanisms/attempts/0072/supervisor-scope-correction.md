# Supervisor scope correction: energy sign versus mechanical force

This append-only correction records the final Sol-High mathematical audit of
0072 before checkpointing.

Pre-correction hashes:

- `derivation.md`: `8742d309527b1bcf8ed67240a66e0831969499faaccecb8680f0b664beec107c`
- `result.yaml`: `75cda69cffe5c182185d332cbb3325acd091f081b392e150715d370bd9592dc5`
- `validation.md`: `f770782233a190cd91b15580501f3bb7ec3f5e072621544970219a347640249a`
- `artifact-hashes.sha256`: `0e68282dd0b1156aad6152f4b9ff814bcd47b6fa276352343afd073d39623a32`

Post-correction claim-bearing hashes:

- `derivation.md`: `f738ebb02d2341d1e51671ac12a6a575eb8f802cec1a9a121ca5d86c53cc965c`
- `result.yaml`: `2f474fd31b157bbcfec5050060ecca7c19b50f125efcf8aeb215addeb64e5703`
- `validation.md`: `1e04402a22c98c1098218f1334f4f00aade5f5401f658b052d458e8363b3d7e6`

The exact positive, decreasing `C/d` kinetic cross energy is unchanged. It is
now called a repulsive-sign effective potential only conditionally on `d`
being constructed as an admissible mechanical coordinate on a finite or
renormalized coadjoint/KKS phase space. The standard identity-at-infinity
compact-carrier orbit does not yet contain an independent translation of one
of the two overlapping noncompact summands, so 0072 does not claim an Euler
force.

The next stationary test is also ordered correctly. The full inverse-VSH
angular field of the fixed-frame symbol must first satisfy the homogeneous
sphere equations. A faster-decaying core cannot repair a nonzero degree-minus-
six residual. Only after that leading test passes may a core be tuned to repair
the global isotropic-stress row and the full steady gluing problem be posed.

No formula, public API behavior, oracle predicate, or focused test changed.
The existing final exact and focused receipts remain applicable and were not
rerun for this wording-only correction.

`artifact-hashes.sha256` is deliberately outside this receipt's self-reference
set. It is regenerated after this receipt and verified independently against
every listed artifact.
