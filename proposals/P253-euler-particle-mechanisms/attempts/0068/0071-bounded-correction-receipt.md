# P253/0071 bounded analytic-domain correction receipt

Date: 2026-09-06

The single independent-review correction makes two related evidence/domain
repairs without changing the equations or route verdicts.

1. The resolved-state map now uses a bounded finite-rank idempotent `Pi` with
   smooth range. `Pi` and `Q=I-Pi` are bounded on the declared dimensionally
   typed product Sobolev spaces and reduce the full joint constraint kernel.
   No untyped `L^2` orthogonality across `(u,chi,E,B)` is assumed or used.
2. The exact oracle is described as an algebraic regression check. The action
   and displayed PDE calculation carry the Poynting and constraint-propagation
   proofs. The propagation speed is explicitly the Maxwell-subsystem speed.

Pre/post SHA-256:

| Artifact | Pre | Post |
|---|---|---|
| `derivation.md` | `26db2c0e6d119644ee7ed9308826e104bf91999ac55d21d1c824d64f9be996bc` | `41bfff94e0c72d0c028564b2ee7c04d231fc3ea3eca75d1d7c556facb86d7816` |
| `validation.md` | `c2e573ee17815473d419cc7d28cab40ba2d4d1c621f28f33e408154de95f1c9a` | `15a83b8888a4737264d3f877ea0030fcbf52e400ce79692c8c587ffcd28c6a58` |
| `sol-high-pre-review-precision.md` | `4c826663db94812e565cef9f5dfe9ac0912989207f403fb3a72fb52b4254f523` | `c41c7236be4f8c8db14655805b28dc353139437139d8f05d41c1356aee2161a6` |

The API, verifier predicates and tests did not change, so their existing
captured receipts were not rerun. YAML, stale-language, newline and scoped
diff checks are rerun after this receipt and recorded by the independent
correction check.
