# P253/0061 bounded joined correction receipt

Date: 2026-09-06

This single append-only receipt covers the independent 0061 review correction
applied to the author-owned `0055` and `0060` claim surfaces. No source module,
test, verifier, captured oracle output, central record, or Git history was
changed. The existing scientific oracles were not rerun because every checked
algebraic identity and the repaired `r^-5`/`r^-6` exponent are unchanged.

## Corrected dependency edges

- `0055` now states only the bare fixed-density whole-space Euler boundary:
  topology and the displayed conserved non-action labels do not universally
  select an action, while prescribing `J` itself fixes `J` tautologically.
  Boundaries, backgrounds, constitutive scales, and hyperbolic extensions are
  outside the route verdict.
- `0060` separates the six signed-axis orbit `O/C4` from a free 24-site orbit.
  The 24-element tensor projector remains exact, while nonvanishing of the
  six-site DA average requires a nonzero `C4`-invariant projected tangent.
- `0060` now displays the exact translated zeroth, first, and second vorticity
  moments and excludes the proper-octahedral epsilon tensor using symmetry in
  the two spatial indices.
- Physical and kinematic pressure conventions are distinguished:
  `M_ik=-(delta_ik/rho_0) integral(p-p_infinity)` and
  `M_ik=-delta_ik integral pi` for `pi=(p-p_infinity)/rho_0`.
- The direct-product DA statement now imposes support/collar, circulation,
  harmonic, centering, and stabilizer rows componentwise. KKS splits over
  disjoint cores; whole-space Hodge/pressure and the kinetic Hessian retain
  cross terms.
- Scalar lattice convergence is preserved. Bloch `C^1` remains conditional on
  normalized invariant modes and a uniform full pressure/Leray
  `D(A_0)->X`, `r^-5` coupling plus graph-resolvent/complement control.

## Pre/post SHA-256

| Artifact | Pre-correction | Post-correction |
|---|---|---|
| `0055/README.md` | `b39546f5f3f34bccf1bae2aa6b3993886ff37a6d9085e0db6585a5233018daf1` | `bf46ea2c9c82d16e26cede0a699447ba79a2e69bdafedf5d1aa38fe6b7cb358a` |
| `0055/construction.md` | `372af8de9c86538155f34e7cff597e5bc523b42370caad00eacb75f177f49c3c` | `5351f8bf3abd7b622c0902419f6c99200e5f142fd48e7e59e492625df41dca96` |
| `0055/result.yaml` | `ac67b78a4cba97b9bab9fe40da1f3832365cd0765a556130dac92312df0ea81e` | `938ad472fe8dad612877276da562285a2b6a9d3f5223ccee906a4ab47e3c81ef` |
| `0055/source-audit.md` | `45d69883d8fa7b7e52901a0f760746ed5d56ff83411ac41d749d27523440600d` | `db20e8e969613341d46db73bdc360d92b8d4d6044275586d97bc73d7720117c9` |
| `0055/validation.md` | `ae9c041f4f6684b30cdddb34b7c62cac9fc9b429e5eb925a86ca79d96ff69989` | `745335c8577b32c07319278c4caeae17614ba13301d9173825df31d0054ace08` |
| `0060/README.md` | `74a72c03f156f931888bb7d6b60b13098b7e703a51b2e53aa4a7948f75413475` | `44ed68373f5307247bb3b79fd1a72646fd84fdc09a95b760f31ecef9cfd63dc3` |
| `0060/construction.md` | `1e5d2c769daf81ea688d1a029fa2f2c46ed5f98b5ad57559b10d1a1062cda372` | `487891b43cf40b810f6d547909adfdc37a23c64c6a2e868d1b007885ca12ec6e` |
| `0060/result.yaml` | `46f8a0da69e30a01bc46cee0d946f0fe3ad3411deba00eefd84d41a57bf459f7` | `0eac38f436c384248e2c69fd31a9d8c0a41e892b6fa6d6d1cb02fc6401de1edb` |
| `0060/source-audit.md` | `d7154821f6953cf0b7e661dd78c90abd4f84426a0d762bfcac71b267d15244f4` | `91afbf5c4a6159cc8dc2cece066377a7c1fa193466de6c64c575ce77c289c06b` |
| `0060/validation.md` | `06d3487954e6af9139b5f8a27423170f413821757630f9cd5a2bc02c4feb84e8` | `bef8cde5d00cd77989d2bff8b28802f78154387dd24404d1e05627e85751215d` |
| `0060/multipole-order-correction.md` | `9b3589b0062f502d4ad91c60a0296e7330309e1a21441ec288fcdbb1dd7f181d` | `2d72ea85ffa3623cc45adbe57a4e4b6bff0bea8fd33af480c6491f4e6fde3640` |

The unchanged authoritative algebraic evidence retains these hashes:

- `0055/verify_scale_causality.py`:
  `51c0d4ae72c49d4eba525f9d55ac1a535d998dac7d34878aaabf0d976212f5b0`;
- `0055/extension-ledger.stdout`:
  `ea56f577a3df62732bcb84c6e363a91d3f890b61851084313290a218916ccbe2`;
- `0060/verify_neutral_cell.py`:
  `895cbd35955209b8b56bab3560b59809affec613190ff26b8ca3cf4c49369c0f`;
- `0060/multipole-order-correction.stdout`:
  `f3a908db95d4fd54c845b0f58223111d596199df74b0f2f568a173f0379e0b07`;
- `src/substrate_framework/euler_neutral_cell.py`:
  `b666df746bf7bb96a04fa7abfa3e08b230622d2c8b89cc6bf5d59fb843f4fe63`;
- `tests/test_euler_neutral_cell.py`:
  `a3b59c32618fa2e5847a3df0aa7bd81ccec258fed0cf232fbe66c24c880da56a`.

## Bounded validation

Using the repository interpreter, both result YAML files parsed and retained
`parent_state: active` and `exhaustion_claimed: false`; the repaired 0060
route and direct-product verdict values were asserted exactly. A stale-language
scan found none of the superseded similarity, nonzero-average, free-orbit,
pressure, or unconditional Bloch phrases. `git diff --check` passed on both
attempt directories. These are textual/schema/diff checks only, not new
scientific oracles.
