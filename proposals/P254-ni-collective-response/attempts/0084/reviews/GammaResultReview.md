## Verdict: correct_at_scoped_claims

**Complete findings array:** empty — no actionable patch-introduced defects found in this assigned scope. No incremental findings were omitted.

**Confidence:** 0.94.

### Evidence and reviewed scope

- **Gamma-only recovery without electron replay or state repair.** `attempts/0084/finish_gamma.py:83–103,177–193` authenticates the original Gamma invocation/input and the completed electron bank. The recovery opens electrons through the existing read-only loader; its sole solver invocation targets Gamma. `attempts/0061/production_transport.py:304–317` authenticates completed bank bytes and maps them read-only. Raw number drift, initial-state defects and density diagnostics remain visible in the actual readout.

- **Authentic input and preserved failure.** The original and recovery receipts both identify `/var/tmp/p254-a84-intrinsic-full31-001/gamma/.staging/slot-000/input.bin`; its recorded digest agrees with A85 preregistration. The original timeout receipt, invocation and packed input remain in that failed stage. Recovery uses a distinct staging directory and authenticates the protected originals again before publication (`finish_gamma.py:193,278–294,337–340`).

- **Actual complete publication, not merely successful Native execution.** The new Gamma receipt records exit zero and 214.55090089095756 seconds. The existing `ResultView` copies the result into the Gamma bank, followed by commit/finalization; the recovery then publishes `provider.json`, opens that actual path through `SavedHistory`, consumes it, and only afterward writes readout/result (`finish_gamma.py:291–296,315–369`). The saved artifacts report 128 electronic sectors, 114 Gamma positions, genuine nodes 0–30 and integration order one. Counts distinguish 128 reused electron solves, one new successful Gamma solve and the retained failed Gamma invocation: 130 total, 129 successful.

- **Legitimate A83 endpoint contrast.** `attempts/0083/projected_mean_control.py:53–76` checks the original reference, coefficient manifest, electronic preparation and full Gamma preparation; it consumes the archived A78 mean/rho nodes and uses the original Ni1/x probe divided by `sqrt(N)`. The coupled tangent retains density/K feedback and adds the common curvature with its true bilinear transpose. The saved control consumer uses the original ordered normalizer. Its published input-mean identity agrees with A84’s provider. The control is explicitly nodal/interpolated, not a continuous-mean match or full114 control history.

- **Scientific and delivery claims remain scoped.** `SOLUTION.md:47–62,98–102` distinguishes the first Picard update from selfconsistency, the bare reference from the node-matched control, and the Gamma contribution from a full spatial coherence curve. It explicitly retains particle-number nonconservation and disclaims temporal accuracy, irreversible material loss, lifetime, rate and heat. The issue/source PR remain open; no merged-delivery claim is made.

### Accounting qualification

The immutable result’s **968.0730148446378 seconds** uses the recovery worker interval, **264.74176404904574 seconds**, rather than the new manager’s **265.1844534496777 seconds**. Thus it is not the full-manager cumulative duration; summing the recorded manager intervals and compile/oracle charge gives approximately **968.5157042452697 seconds**. The prior failure remains charged, and Native time is separately identified as overlapping rather than added twice. Root’s announced supplementary accounting record was not yet an inspected artifact in this review.

### Exercised verification

**None new.** Review consisted of reading the patch, implementation/consumer paths, existing receipts and saved publications. No files were modified; no build, lint, formatting, test, benchmark, Native execution or numerical re-consumption was performed.