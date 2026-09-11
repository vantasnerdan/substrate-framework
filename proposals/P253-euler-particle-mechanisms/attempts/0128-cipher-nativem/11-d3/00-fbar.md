# D3 F-bar freeze (pre-compute; D-08 form) — multi-period geometric phase

Thesis (SYNTHESIS.md unrouted breaker): if charge is a PER-PERIOD phase, all
short-window response fits (S1/B2) were structurally blind. Construction:
shot orbit integrated over N_per = 5 consecutive periods (section-flow,
N=64); transverse Floquet eigenvector parallel-transported period-to-period
via overlap chaining (N4 Berry machinery along TIME); observable = geometric
phase per period γ₁ (mod 2π) + coherence across periods.
- ALIVE-leaning: |γ₁| > 0.05 rad stable across ≥4 periods (std/mean < 0.5) →
  per-period phase EXISTS (charge candidate revived at multi-period level;
  P2/AB-adjacent leg opens; nothing claimed beyond existence).
- DEAD: |γ₁| ≤ 0.01 rad (floor: overlap-angle resolution ~1e-3 × periods) →
  no per-period phase; short-window blindness moot (nothing to be blind to).
- INCOHERENT-dead: std(γ₁ across periods)/mean > 1 → phase random-walks.
- Gray otherwise → UNRESOLVED + named leg (more periods, N=128).
- STOP (frozen): monitor section return each period; closure res > 1e-6 →
  break + report horizon H (verdict on completed periods only, H stated).
- Overlap gate: adjacent-period eigvec overlap ≥ 0.99 else INVALID link
  (re-examine, no verdict) — same gate as N4.
P2-adjacency (honest scope): a stable γ₁ is AB-ADJACENT (phase per period),
not AB itself (no circuit, no flux quantum) — supports P2-structure at
multi-period level; dynamical charge still missing-5. Run reports
(γ₁ per period, coherence, horizon, verdict), nothing else.
