# A87 Sigma code review

**Verdict:** correct within the three-file source/evidence review scope.  
**Complete findings:** []  
**Confidence:** 0.90.

No reportable patch-introduced defect was identified. This is an implementation
review verdict, not a completed-simulation or numerical-accuracy verdict.

## Scope and reviewed original bytes

Repository: `/home/administrator/projects/substrate-ni-collective`; branch:
`research/P254-ni-collective-response`; canonical issue: #225.

All paths below are relative to `proposals/P254-ni-collective-response/attempts/`.

| Reviewed file | Bytes | SHA256 |
| --- | ---: | --- |
| `0087/dressed_ee.py` | 10033 | `d21c9b3a73c036de17eae99f71156e751b84aeae27d09e80438f1491bb0415cb` |
| `0087/produce_component.py` | 18795 | `58a59c90f015f707c4843d940441055f0ea0d755fe0a2ce0f483facd3a9d8336` |
| `0087/component_store.py` | 11198 | `8f2267490a61997b6ebd7b5ace6901317fd63d25645562f11e61fb0a1d340e68` |

## Source review conclusions

- EE retains all64 inner-h contributions and independent all128 forward/reverse
  electron callbacks. Greater/lesser orientations and theta(0)=1 retarded masking
  match A60/A46. The active EE path does not substitute Gaussian U separation.
- Direct selection retains the full128-slot bubble. Contiguous/noncontiguous
  static-row selection preserves external-major/full64-h matrix/operand routing.
  EI maps each missing target back to its original source and exported legs.
- Committed outputs are excluded. Resume verifies committed digests and rolls
  ahead-of-checkpoint bitmap bits back to checkpoint ownership. Numeric payloads
  and digests are flushed/fsynced before bitmap ownership and atomic checkpoint
  publication/directory fsync; failed commits prevent subsequent publication.
- The A76 consumer processes the emitted Sigma layout and coverage. No unmatched
  consuming-side dispatch or silently discarded EE/EI output was identified.

## Evidence read

Complete new sources and patches; A60 streamed EE; A46 BLAS/shared contractions
and physical validation; A62 history sources, functional, and banks; A63 actual
`DressedHistory`; A76 emitter/consumer; owning material/Coulomb loaders, factors,
EI edge mapping/monomials, and source-descriptor parser. Also read A87
preregistration, actual-G1 selected-direct oracle, descriptor recovery, relevant
actual-provider metadata, and the first failed EE manager receipt/traceback.
The JSON report names the individual evidence paths.

## Evidence limits and scoped risks

- The existing actual-G1 oracle covers q0, pair (19,7), both contour orientations,
  and reordered/noncontiguous plus tail direct outputs: zero reported error at
  rtol=1e-12/atol=0, with wrong-row mutations rejected. It is not full-exchange,
  full-component, other-q, or crash-injection execution evidence.
- Descriptor recovery reused the owning parser: original bytes/SHA match, legacy
  string/integer counts authenticate, and independent wrong-SHA/wrong-count
  mutations were rejected. No re-pin is evidenced.
- The first EE manager failed before contraction with zero numerical output cells,
  costing **93.64355171797797 seconds**. Corrected EE002/EI001 being reported
  running at review time does not establish numerical completion.
- No tests, imports, builds, simulations, Native runs, source/data edits, or
  old-suite replays were performed. Only these requested reports were written;
  coreutils supplied byte counts/SHA256. Persistence claims above are source
  inspection conclusions, not exercised crash-recovery results.
- No complete collision map, next Dyson update with M2, accuracy, conservation,
  calibrated loss/rate, or closure claim follows from this verdict.
