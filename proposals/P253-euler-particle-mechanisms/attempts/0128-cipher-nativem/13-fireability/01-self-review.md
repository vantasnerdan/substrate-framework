# 13-fireability — same-agent review (designer reviewing own protocol)

Status: SAME-AGENT pass, not an independent firewall. Drift rebuttal invited;
any drift finding overturns this verdict. Scope: review ask §§1/2/4/5 in
`00-fireability-protocol.md` (#9B02).

## Findings (5 repairs, all doc-level, none structural)

**F1 (§1 loophole — saturation anchor).** R1/R2/R4 calibrate on the p → 1
saturation asymptote, but the doc never says in which channel saturation is
identified. If the asymptote is located by wave-signal plateau (signal stops
growing with poling field), wave-channel data feeds p* — circular by the
doc's own §1 rule. Repair A1: the asymptote must be identified IN the
independent channel (Δn / Δε plateau), never in wave signal.

**F2 (§2/§4 inconsistency — R0 zero).** §2 claims R0 needs "no measurement";
§4(b) admits it needs an upper bound p < p_lo. The model p is a LOCAL
correlation parameter; macroscopic unpoling does not guarantee local
decorrelation, and residual local p gives a weak signal that confounds the
null. Repair A2: R0 = preparation symmetry PLUS an independent-channel zero
reading (cheap null on R1/R2 at the same preparation, not a precision
measurement). The two sections are reconciled on this form.

**F3 (§2 — R3 common mode).** R3's p-map rides the declared dipole
statistics — the same premise D1's K_n derivation rides. A statistics-model
error shifts formula and p-map in correlated ways: common-mode failure the
tier-1.5 flag alone does not contain. Repair A3: R3 is never the sole route;
every FIREABLE return pairs R3 with R1 or R2.

**F4 (§2 — R4 non-cancellation).** Two-state ratios cancel K, C_c, κ_n, ξ
only at leading order under IDENTICAL mechanical state; poling across
preparations can shift pre-stress/gap terms (C_c, K_p λ_±) that do not
cancel. Repair A4: R4 ladder poles the SAME sample (or bounds ΔC_c/Δgap
across samples inside H); the build charter prices the non-cancellation.

**F5 (§4/§5 gaps).** (a) Budget ledger omits the instrument-noise floor —
detection limits enter O_null alongside H. (b) §5's UNFIREABLE return should
name WHICH leg failed (route precision vs preparation window vs
hygiene-dominance) so the failure constructs the next attempt.
Repairs A5a/A5b accordingly.

## Verdict

CONDITIONAL PASS → all five repairs applied as amendments A1–A5 in the
protocol doc (original §§ frozen intact, append-only). Post-repair verdict:
PASS (same-agent). Mu-ratio CIRCULAR verdict, blinded-insertion bar,
tiering, and the §4 inequality stand — tightened, not weakened.
No FB-D verdict touched; no numbers produced.

Rebuttal surface for drift: F1 (is in-channel plateau sufficient, or does
the asymptote need an absolute standard?); F2 (does the zero-reading
reintroduce the channel it was meant to keep cheap?); F3 (is the never-sole
rule strong enough for common-mode statistics error?); F5b (is leg-naming
sufficient failure-derived continuation?).
