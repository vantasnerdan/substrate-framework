# P253/0093 bounded correction receipt for P253/0092

This single bounded correction repairs the charge/action dimensions and
notation without changing the topology, BF, continuous-tag, compact-character,
or fixed-charge rank verdicts.

The transported-tag coupling is now `g_tag`.  The compact phase uses a
dimensionless number density `n`, imported action coefficient `S_0`, physical
charge unit `g_0`, integer character `m`, and
`kappa_m=g_0*m/S_0` with units charge/action.  Thus

    J_phase=S_0*N_phase,
    Q=g_0*m*N_phase=kappa_m*J_phase,
    J_*=4*pi*epsilon_EM*c_EM/kappa_m^2

in the explicitly imposed common-phase constraint.  No map identifying
`g_tag` with `kappa_m` is claimed.  Under `A'=aA`, both `g_0` and `kappa_m`
scale by `1/a` while `S_0,m` remain fixed.

Pre-correction and post-correction SHA-256 values are:

| artifact | pre | post |
|---|---|---|
| `README.md` | `ad62ed84b161016031c3dc0a6f69fe1267843665752c02c370fed952ddfd705f` | `2e8138c4871f6548fc9f8a215fd747252a33bf2d4f6099afc1d7f04c7ae85dd0` |
| `derivation.md` | `76960dc6092a23db6838098ca4e9ea12cabdb62f34ac0eafd2e588b5aaac06d1` | `3ce9062d93a43ea3d3f65eb817ff2fc3d1c9e98da47f83ded583ea4e184f7cbc` |
| `source-audit.md` | `0c63315bf3fb9207599f1c712c24c32d075cfa8a9422547dd084d5e6bcdf8fe8` | `dd170288f4f7d3b82b214532d9bd808d5b556c2afc15dd0cd9afa5aee48f3daa` |
| `result.yaml` | `47765950a9318b968c5a30b31fb00f2ac58d83a2eaa1748a230cce2e61bb36d8` | `a28ea21fcbc42b383489089ff475176a54d63d98313e2ba513075436f7b5b678` |
| `validation.md` | `41900ddbea87b042f9846a1aacdaac50d536c96cfaeb3a4bb3fa3a4f624a31a9` | `29797a58b37fa5f587769c57c9250f4cdb135f45f623f6dfd2fe759b95edd4e6` |
| `verify_action_selection.py` | `ebd18ba2f2990f333cb237805eae9b7d3437d7b41d992f0029583acd5a295b20` | `86ca3ae3859abb42ce2db0a6c572b96e62d498ebeefc6f97caa79a48257b2509` |
| `src/substrate_framework/euler_gauge_action_selection.py` | `ad00d857f115cb9e5153e02b6e2def6d8dcbc784533d848a381d1b57b31f0560` | `c52dc54d5e51f85cf65ffca95ccac7eed518c6ab3c6717f3c5bb63c173854172` |
| `tests/test_euler_gauge_action_selection.py` | `6a39053c3187287d25d6f4e13fe5ad52e3584b1f88b828447d5fbb4db00590fa` | `807d345b9d311df88a98dd893b4d32f3e147029e7b8d5817ea6b4654a673089d` |
| `exact.stdout` | `81413b850eed30387d2686b93b10a8b41fba6356ea242b1773f58f4c9674ea43` | `92c7c28a986de1da2e5283b2afce79ada72f4b4e08dadeca3e2bf5430ed14294` |

Because production API bytes and the exact verifier wording changed, the
focused API tests and exact symbolic verifier were rerun.  They pass with
eight tests and twelve `PASS` rows respectively, both exit `0`; the exact
stderr is empty.  The repository validator was also rerun and exits `0` with
the standard `WORKFLOW VALID` line.  No production numerics or unrelated test
suite ran.  The final manifest is intentionally outside this receipt's
self-reference set and is refreshed after this receipt.

The correction-only check found one remaining untyped `Q/g` phrase in the
final dependency list.  The same transaction changes it to `Q/g_0` without
altering mathematics or encoded predicates.  The derivation hash immediately
before this clerical closure was
`580256346e30cc520d60fd117329a5f6d592b593848de27e2f65ee99829e4057`;
the final hash is
`3ce9062d93a43ea3d3f65eb817ff2fc3d1c9e98da47f83ded583ea4e184f7cbc`.
No API, test, verifier, oracle, or validation command was rerun for this
wording-only closure.
