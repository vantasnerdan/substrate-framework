# Gap closure (hand-maintained; board renders it, health checks it)

One row per open terminal-PR gap: the exact next artifact, who owns it, and
the STATUS `blocked-on` token(s) it discharges (copied verbatim). Keep rows
current on landing; never delete a row — mark it CLOSED with the closing
artifact. `owner`: beacon/cipher/drift/shepherd/owner(triage: user ruling).

| gap | next artifact | owner | class | waits (verbatim tokens) | status |
|-----|---------------|-------|-------|-------------------------|--------|
| G-a2 numerics | attempts/0111-beacon-ga-field G-a2 charged-branch member build (field arrays + norm certs; frontier ga-status.md) + attempts/0120-beacon-trust/trust-report.md | beacon | agent | G-a2-branch-numerics, bg_7-trust, G-a2-fitted-mesh-or-errorbars | trust rows-met 6.6e-3, floor structural p≈0.4; next: fitted mesh or errorbars; IDEA-03 1/3 closed |
| R-EM2 decision | owner approve/amend/decline of attempts/0112-cipher-rem2-draft | owner | owner | owner-review, R-EM2-decision, R-EM2-import | awaiting ruling; drift 0112 draft-technical PASS banked |
| S9 test | attempts/0114-beacon-s9/s9_probe.py (+design.md) | beacon | agent | S9 | in progress (dir landed, no wait token yet) |
| shadowing scope | attempts/0113-cipher-shadow (README + receipts) + attempts/0120-cipher-m2b1 (A1-H5, A2-H1) | cipher | agent | shadowing, M2-B1-proof, M2-B1-H4, M2-B1-H4-proof | A1 H5 PASS, A2 H1 repaired PASS; A3 queued; current token H4-proof |
| Euler persistence | joined-PR supervision bundle (tracks G-a2 + R-EM2) | shepherd | owner | Euler-persistence | open supervision umbrella, not a build |
