# Tool receipts — first execution (beacon 0108)

Reuse rule: receipts valid while code/inputs/proposition unchanged (task brief §4).

## T1 — branch / tree / herd (2026-09-10)

Command:
`git status --short --branch; git log -5 --oneline --decorate; git branch --show-current; herdr agent list | head; pwd`

Stdout (abridged, full in job bg_1):
- `## research/203-euler-particles-herdr-resume`
- `?? proposals/P253-euler-particle-mechanisms/herd/`
- `b4731855 (HEAD -> research/203-euler-particles-herdr-resume, origin/research/203-euler-particles, research/203-euler-particles) Adopt main discovery workflow; preserve P253 pause and research (#206)`
- `research/203-euler-particles-herdr-resume`
- herdr list: shepherd (w5:p1), atlas (w5:p3), … (truncated at 768 chars by runner)
- `/home/dan/substrate-framework`

Stderr: none observed. Exit: 0 (job bg_1 completed).

## T2 — memory + graph freshness (2026-09-10)

Commands:
`memory --version` → `memory, version 0.2.0`, exit 0.
`node .gitnexus/run.cjs status` → `Repository: /home/dan/substrate-framework`, `Indexed: 9/9/2026, 1:09:30 PM`, `Indexed commit: 98c155f`, `Current commit: b473185`, `Status: stale (re-run gitnexus analyze)`, exit 0.
Delivered via job bg_2. Graph STALE — named limitation (see README).

## T3 — memory search, particle routes (2026-09-10)

Command: `memory search --base "$PWD/memory" "P253 Euler particle vortex ring" | head -n 80`
Exit: 0. Top hits:
1. `codex/efforts/P253-mathematical-supervision.md > Exact Evidence and Remaining Achievements (73.76)` → `attempts/0037/supervisor-projection-analysis.md` accessible counterexample + bounded-co…
2. Same file `> Goal and Success Contract (43.75)` → pause 2026-09-07; recovery entry `PAUSED.md`.
3. `codex/proposals/P253-euler-particle-mechanisms.md > Sol-High 0088 freeze… (41.78)` → 0095 activation SHA `be4a0c67…`.
4. `… > Sol-High restart recovery… (40.75)` → joint charged Euler–Maxwell Hessian + doublet frontier.
5. `codex/efforts/issue200-continuation.md (40.02)` → 0248 balanced-ring self-frame deficit.
(Delivered as bg_3 completion notice.)

## T4 — memory grep, retained claims (2026-09-10)

Commands: `memory grep --base "$PWD/memory" -F "0042" | head -n 40`; `memory grep … -F "0095" | head -n 40`. Exit: 0 / 0.
Key outputs:
- `codex/proposals/P253-euler-particle-mechanisms.md:145`: "Root0042 constructs an exact finite-rank resolved/complement Euler state … Independent0045 establishes the joined result after one bounded analytic-domain correction …"
- `codex/proposals/P253-euler-particle-mechanisms.md:676,703`: "P253/0095 is centrally activated only against corrected README SHA-256 `be4a0c67…`."
Stderr: none (truncation notice only: "[Some lines truncated to 768 chars]").

## T5 — memory search, carrier suppliers (2026-09-10)

Command: `memory search --base "$PWD/memory" "Hill vortex Choi Cao vortex ring Gavrilov compact" | head -n 60`. Exit: 0.
Key hit: `codex/proposals/P253-euler-particle-mechanisms.md > Hyperbolic compact-carrier verdict… (42.94)`:
"Independent0038 establishes the joined0032/0039 result … On the same fixed compact Gavrilov carrier, the exact whole-space dynamically accessible linear Euler group … obeys `||S(jT_*)||_ess>=lambda_+^j` … adverse evidence for that carrier as an LP2 particle sector; linear semigroup theorem, not nonline…".
Second search `memory search … "Cao thin ring linear growth 0095 0104"`: hits on 0105 lock classification, 0038 verdict, 0053 corrected finite-window theorem, 0041 quantum-boundary checkpoint. Exit: 0.

## T6 — GitNexus query/context (2026-09-10)

- `node .gitnexus/run.cjs query "Hill vortex ring …"` (no `--repo`): ERROR `Multiple repositories indexed. Specify which one with the "repo" parameter. Available: amazon-new-best-sellers, substrate, substrate-framework (/home/dan/substrate-framework), …`. Exit: 0 (process exit; tool-level error).
- Same with `--repo "substrate-framework (/home/dan/substrate-framework)"`: ERROR `Repository "…" not found. Available: (same list)`. Exit: 0.
- `node .gitnexus/run.cjs query "Hill vortex" -r "substrate-framework"`: SUCCESS, `{"processes": [], "process_symbols": [], "definitions": [abelian_higgs_vortex…, vortex_dynamics…]}` — no indexed Hill-vortex flow. Exit: 0.
- `node .gitnexus/run.cjs context "src/substrate_framework/euler" -r "substrate-framework"`: `{"error": "Symbol 'src/substrate_framework/euler' not found"}`. Exit: 0.
- `node .gitnexus/run.cjs query --help`: usage confirms `-r, --repo <name>`. Exit: 0.
Conclusion: graph stale (T2) + empty on this concept (T6) → AST/grep + direct reads used; limitation named.

## T7 — repo grep, Euler modules (2026-09-10)

Tool: `grep` path `src/substrate_framework`, pattern `Hill|vortex_ring|Gavrilov|impulse|helicity`.
Hits (with snapshot tags): `euler_action_locking.py#B6F0` (impulse `pi·rho·kappa·R²`),
`euler_asymptotic_tails.py#78E3`, `euler_cao_schur.py#7B60`, `euler_charge_multipole.py#3848`,
`euler_impulse.py#3FB8`, `euler_p2_principal.py#C40E` ("do not assert … coercive Hill-vortex Hessian … evidence is the active P253/0095 attempt"),
`euler_scale_causality.py#D3EA`, `euler_twisted_carrier.py#8000`, `vortex_dynamics.py#312D`.

## T8 — source reads (2026-09-10)

- `src/substrate_framework/euler_cao_schur.py#7B60` (full, 93 lines): leading jet only; "do not prove the Green-operator Fredholm theorem or construct a charged Euler–Maxwell branch."
- `src/substrate_framework/euler_impulse.py#3FB8` (full, 51 lines): dipole + cross energy; "caller supplies actual impulses…; algebraic evaluation alone supplies no carrier or stability."
- `src/substrate_framework/euler_p2_principal.py#C40E:1–73`: "expose algebraic identities only… do not assert generator domain, finite-curvature Cao expansion, nonlinear persistence, or coercive Hill-vortex Hessian."
- `src/substrate_framework/euler_scale_causality.py#D3EA:1–120`: similarity weights + pressure quadrupole; "do not select an action quantum… do not turn an effective carrier band into a relativistic causal cone."

## T9 — 0104 independent review (resume input, re-read 2026-09-10)

- `proposals/P253-euler-particle-mechanisms/attempts/0104/review.md#B695` (full, 327 lines): Units A–F established at stated scopes; G blocked (total-momentum domain + finite-row submersion); H established (fixed-delta small-charge linear observed growth); I blocked. Combined: linear A–F theorem stands; nonlinear refutation blocked.
- `proposals/P253-euler-particle-mechanisms/attempts/0104/verdicts.yaml#8B78` (full): machine-readable unit verdicts + `next_dependency` (Unit G witnesses + complementary column branch) + `does_not_license` list.

## T10 — memory validate (2026-09-10)

Command: `memory validate --base "$PWD" "$PWD/memory" | tail -n 15` → `All 1068 file(s) valid.` Exit: 0.

## T11 — attempt dirs (2026-09-10)

Command: `ls proposals/P253-euler-particle-mechanisms/attempts/ | head -n 30; ls …/attempts/0042/ | head -n 20; ls …/attempts/0108-beacon-sources/`.
Output: 0001…0030 present; `attempts/0042/` contains `README.md, construction.md, result.yaml, review-correction-0045.md, validation.md, …`; `attempts/0108-beacon-sources/` absent (`No such file or directory`, exit 2) → created by this attempt.
