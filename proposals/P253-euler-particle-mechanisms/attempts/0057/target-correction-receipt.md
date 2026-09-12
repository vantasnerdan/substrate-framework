# P253/0057 joined target correction receipt

This is the single bounded correction package requested by the independent
0057 reviewer. It changes claim classification and analytic-domain wording in
0052/0054; it does not change the exact formulas, verifier, API, or tests.
Accordingly the unchanged symbolic oracles were not rerun.

## P253/0052 hashes

| artifact | before | after |
|---|---|---|
| `README.md` | `219334fcd5bbf40f69d691dc45596a82926b5f25e63836dc3a9eee5990ba6bd1` | `f534e6ee4b1dabb0ea3e0fdb9a0223139f5deb860571b1ee98229f06093d28ee` |
| `derivation.md` | `d353a4ffb8ea1ce65b419f8152b3292add53671a3d5e94195685722ae95f1d5b` | `6ea25cd99b91686d0d7efada0bc756656b82b4ebe44f310848511ab63da23e76` |
| `result.yaml` | `9bab0def7992099eace135783d77e4e2a84f96f55792089d0f51177912afc8a8` | `1d0781660745266a447bf8ecf2f2e67c7804bfb1572060ad05ce56ec7824cc7a` |
| `validation.md` | `e7a9ee362443ceca42dbb631977ad87218fe6a3353ff9722079007f7810de509` | `054dd0f82cf17afc2b444d545677fbbdf5f6d1f4419cf35548bf4d55fd6eebbe` |
| `source-audit.md` | `1dd942371cf4591ad991dfe9bff8ab67b9a64895ffc84f06442adc04dfd558c2` | `1da32efb5f36267c3e44717461fb178a3c3d3f1ea5c759f8db6f89275a6dace9` |

The corrected resonance inversion, logarithmic scaling, exact local Cao cells,
unequal-volume Hanzawa/Piola representation, toroidal distance, and physical
KKS normalization are unchanged. `GR`, the nonnormal all-sector graph
resolvent/Grushin estimate, is now an independent hypothesis in addition to
`HJ2`.

## P253/0054 hashes

| artifact | before | after |
|---|---|---|
| `README.md` | `a13c902a9798e48897aa787f3eecdaf331843a78d0ade783690d8ca476e96a8c` | `973729e59d2acbdf295eff2ddba3e95a9421bea08d886488bc44256361319e0a` |
| `derivation.md` | `62e3d45cd8c45099c51b47467339a0824b3fdc1ddedcc371645258723575d760` | `6d4835ae109395c5e7c4be4d8d9fa285d634879e1be10ae383f350b464ca54ee` |
| `result.yaml` | `f060d63e92cdf69d414a0a795dcbd35d4e5e17da6ae0557de6d5a7e4d28eebe6` | `ad6e7a19f519294d63a404c274eafe9fa3171a7dc4de65131012200cdeb6c947` |
| `validation.md` | `e5ee164a75b4e75a13a48d2056c9a1f15d867c3a96ae7902881522a9ec88e829` | `0f648c06e2d5221248cb87d7cb399b5e93abe0e6502db374a5ec7e8025494543` |
| `source-audit.md` | `4ba9323c8f2e16f3418d1809d1245b49e9d688f924f7a50d6b6b952e6c508f6f` | `9d5c7de1cf150438a3a610e214636855f0ea2d590340257af747c00c4dd89519` |

The exact pulled-back Biot--Savart first/second shape identities, exact
toroidal distance/phase, local Cao cells, metric-Leray algebra, and limited
`p>=6` no-sheet trace fact remain. The exact augmented carrier branch, common
closed DA graph, symbol/two-index Schur estimates, and norm-resolvent/Riesz
transfer are now explicitly blocked or conditional on their named analytic
constructions.

## Local integrity

Both result YAML files parse. The scoped `git diff --check` passes. The stale
promotion scan finds no claim-bearing statement that 0054 discharged `HJ2` or
that 0052 constructed the full nonnormal graph inverse; preregistration target
language remains historical contract text and is followed by an explicit
post-execution correction boundary.

## Clerical reconciliation requested during the bounded check

The first corrected source-audit hash
`7211f38f9e1a1773868200ed155edd3362bccc94935beb1c318c8c2295acf1ed`
still contained one stale sentence saying differentiability was obtained after
constructing a local branch. That sentence is removed; the final 0054
`source-audit.md` hash is `9d5c7de1cf150438a3a610e214636855f0ea2d590340257af747c00c4dd89519`. This is part of the same correction
transaction and changes no mathematics or oracle.
