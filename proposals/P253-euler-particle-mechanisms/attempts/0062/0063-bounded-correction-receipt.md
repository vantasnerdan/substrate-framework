# 0063 bounded correction receipt

Date: 2026-09-06

This is the single append-only correction package requested by the independent
0063 review. It changes no established 0062 algebraic identity and does not
rerun the unchanged oracle.

## Scientific correction

1. Equation (7a) now converts cylindrical coordinate coefficients with
   `D(r')=diag(1,r',1)` before the orthonormal rotation and applies
   `D(r)^(-1)` on output. The abstract Cartesian `B_R3` block and signs are
   unchanged. At the axis, `D(r)^(-1)` denotes the regular physical-harmonic
   compatibility limit, not uncontrolled division by `r`.
2. At the frozen integer index `s=4`, `p>=6`, and fixed `|n|>=2`, an exact
   divergence-free h-by-h collar sequence proves that ordinary ambient
   vorticity and generator graph norms do not bound the displacement quotient.
   The quotient lower bound uses harmonic orthogonality, positive-core
   injectivity of `C_0`, finite-row removal, and support separation from
   exterior stabilizers.
3. The Hodge row uses the exact identity
   `B_R3 curl(xi cross omega_0)=P_L(xi cross omega_0)`, including compact-curl
   zero moment and the decaying low-frequency convention. No unconditional
   inhomogeneous Biot--Savart gain at zero frequency is claimed.
4. Same-character finite rows are removed/restored by finite-dimensional
   profile corrections. Bogovskii repairs only divergence introduced by a
   cutoff or profile correction and is chosen to preserve those rows. Both
   `q_*` and `V_*` remain undefined; the next trace may start from one fixed-n
   seed and its first source-bearing resonance.

The positive-core `zeta_min^(-1)` inverse, weighted orbit topology, and
source-specific sandwiched trace remain separate live routes. The collar proof
uses the exact Cao positive-part source representation to derive finite upper
derivative bounds; it does not differentiate a bare `O(d^p)` statement and
does not assume a nonzero boundary coefficient.

## Pre/post hashes

| Artifact | Pre-correction SHA-256 | Post-correction SHA-256 |
| --- | --- | --- |
| `derivation.md` | `df197caa166fcdd0c862970e2af91eff0c1c20b93964efa81eccc97830dff0ca` | `f8e62434a0b52e37ff89bd5d57bcfdecc4859e7089ab6a64d8e9a4b90101b449` |
| `result.yaml` | `a71d2f1546de2a5b8e44a1d563a7b22c76f93e123c48ebc22fc74acb69217f44` | `d559a4c3b559c245be2994725aff57a6b36307daa6ad140ce379b9c9b5a82d9d` |
| `source-audit.md` | `92d6afd3915c491faeeb5ed489637d9fee2b037f140ab37a35715b77b9546eb0` | `a6119cfe92d1d2ee05e639cd0938720a1082224716260425aadaac86c2ee738b` |
| `validation.md` | `e727081186659d8c2b1ab3d1075be9308ded3516ae734887d1a35e9f69d5db43` | `57cc0ed8e6f47894cdda0a8157d32258d34d2d7b1f3f29b4d89b46c0e6551f77` |

## Unchanged executable evidence

- `verify_block_topology.py`:
  `cd3c2fcc4ab612f73d7b393b4885bc523316e977184aa9b7e51efae0194d62ad`
- `block-topology.stdout`:
  `44bf4ba1c7fc71bb27dd55bc678d7f1b083987557eebbc360b7fa4aef421b14a`
- empty `block-topology.stderr`:
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- `block-topology.command.txt`:
  `f68635603fe5076f740d781d64a41be10747db7dc01a2378292862dc8c40ced6`
- `block-topology.exit`:
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`

The oracle predicate did not change, so no rerun was performed merely to
create a fresh tally.
