# Tool receipts - 0125 (beacon)

- Scale sweep: inline heredoc `- <<EOF` driver (exit 0) reusing
  sharp_dQ.build_F/Q_of_F + x18_probes.nodal_gradient + ga_pipeline
  on trust-r3 tensor mesh. Printed: 9.1750 / 3.4895 / 1.5141 /
  0.7030 / 0.3383 at scales 1/.5/.25/.125/.0625; xnodes 1/1/0/0/0.
  (RuntimeWarnings in build_F revolve: pre-existing np.where
  divide-by-zero-evaluated-both-branches pattern, also in 0117 feed;
  values unaffected.)
- Asymptotic slope probe (inline, exit 0): exact = 0.1054 / 0.0261 /
  0.0052 at sc .02/.005/.001 (slopes 5.27/5.22/5.21) -> exposed the
  /6 error in G1 v1 (0.87 story inconsistent with slope 5.21).
- FD-vs-einsum check (inline, exit 0): diff O(eps) exact-theory
  convergence (1.04e-3/1.04e-4/1.04e-5) — pipeline exonerated, error
  localized to the dzeta factor.
- Corrected reruns: sharp_dQ.py tensor+trust-r3 (exit 0): G1 =
  5.2052, margin x2.1. sharp_dQ.py fitted default (exit 0): G1 =
  4.1588, margin x2.6. G2'/sweep numbers unchanged (never used /6).
