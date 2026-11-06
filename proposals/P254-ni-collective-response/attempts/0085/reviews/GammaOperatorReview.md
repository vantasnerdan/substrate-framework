## Verdict: correct_at_scoped_claims

**Complete findings array: []**  
**Overall correctness:** correct  
**Confidence:** 0.96

No patch-introduced, actionable defect was found in the A85 numerical operator, packed codec, specific oracle, or actual Gamma ABI integration.

### Numerical operator

- **Retarded solution retained.** `attempts/0085/real_axis_keldysh.hpp:127–142` constructs the two independent kernels, \(F=-(Q\Sigma)^R\) and \(F^{cc}=-(\Sigma Q)^R\), and preserves A54’s startup call and subsequent reversed-kernel timestep argument order. This matches `attempts/0054/supplied_keldysh.hpp:39–45`.
- **Both noncommuting companions retained.** The distinct \(A=R\Sigma^R\) and \(A^{cc}=\Sigma^R R\), complementary lesser histories, and complete initial-source sandwich remain in `real_axis_keldysh.hpp:144–164`. The noise term and independently evaluated companion remain at lines 166–174.
- **Startup continuation retained.** The reusable timestep covers `max(n, order)+1` entries; complete companion histories exist before sandwich evaluation (`real_axis_keldysh.hpp:100–115, 151–157`). This agrees with the pinned upstream helpers’ future-node accesses through the integration order in `cntr_convolution_impl.hpp:1245–1318, 1396–1494`.
- **Removed imaginary work is structurally zero for this ABI.** P254VIE1 supplies only R/<. The original codec clears all other channels, and the upstream mixed convolution is \(A^{tv}B^M+A^RB^{tv}\); the extra lesser contribution uses mixed legs (`cntr_convolution_impl.hpp:901–977, 1076–1153`). Thus these zero channels do not supply a missing preparation term. The new type has no imaginary accessors and throws on absent real channels.
- **Unused diagnostics do not feed the result.** A54 accumulates diagnostics and residual products after constructing the physical result (`supplied_keldysh.hpp:58–84`); A61 discards the returned diagnostics. Their removal does not remove a result correction.
- Orders 1–5, statistics ±1, complete uniform histories, and dimension-dependent matrix operations remain supported without a Gamma-value special case.

### Codec and consuming side

`attempts/0085/lean_result.cpp:12–44` retains the shared header validation, exact input extent/caller-cap checks, original times, component ordering, and result-only R/< output. The direct codec preserves row-major triangular storage, finite-value rejection, and lesser-diagonal rejection without repair (`real_axis_keldysh.hpp:54–78`).

The consuming boundary was inspected: `attempts/0084/finish_gamma.py:277–296` invokes the binary, imports through existing `ResultView`, commits, and finalizes. `attempts/0061/production_transport.py:123–166` explicitly validates and consumes both packed channels; neither is silently dropped.

### Existing oracle evidence

The oracle checks exact header/times/extent and finite R/< entries before applying its entrywise predicate (`oracle_real_axis.py:54–76`). It hard-binds the original A61 binary and requires the drop-initial, drop-noise, and wrong-order mutants to fail decisively.

The persisted `/var/tmp/p254-a85-real-axis-001/oracle/oracle-report.json` records:

- **20 accepted comparisons:** both statistics, orders 1–5, startup-only and propagated histories.
- **Zero reported R/< differences** in every comparison.
- **30 rejected mutants**, with minimum maximum-tolerance ratio approximately **1,829,338**.
- Expected rejection of existing output, insufficient caller cap, truncated extent, and nonfinite input.

These are finite-fixture results, not universal numerical qualification. The oracle uses dimension-two fixtures and beta 1.

### Actual Gamma evidence

The retained native receipt records exit 0 in **214.55090089095756 s**. Its stdout reports `nt=30`, dimension 114, statistics +1, and **206,272,792 output bytes**. The Gamma commit receipt binds the same A85 binary as the oracle and the original **412,545,304-byte input**; bank finalization and the completed provider/readout result are present.

This establishes actual ABI consumption and completion—not equality against a completed old Gamma result, which is unavailable after the old timeout, nor a measured speedup or continuum/selfconsistency/conservation claim.

**Reviewed scope:** A85 operator, entrypoint, CMake target, oracle and settled evidence; A54/A61 references; relevant pinned upstream convolution/VIE formulas; immediate result-consumer boundary.  
**Newly exercised verification:** none. No builds, tests, benchmarks, Native runs, or file modifications were performed.