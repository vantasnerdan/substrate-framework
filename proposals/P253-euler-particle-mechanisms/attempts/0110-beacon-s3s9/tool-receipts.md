# Tool receipts — 0110 (beacon)

## V1 — cache re-hash (2026-09-10, exit 0)

Command: `sha256sum /tmp/.../2206.10165.pdf /tmp/.../2011.06808.pdf
/tmp/.../1810.08020.pdf /tmp/.../hattori-fukumoto-2003.pdf`
Output: `6d90be6b…528ca` (Cao) / `0b2fe633…832783` (Choi) /
`fcaca85f…e963a0ebc3c9` (Gavrilov) / `c6a35c44…905951550a7` (HF) —
ALL FOUR match the hashes recorded in 0095/0107 source audits.
Cache dir listing: 14M, 24 PDFs across P253-0005…P253-supervisor; S3/S4/S5/
S6/S7/S9/CLV absent (gap this attempt fills by fresh fetch, except S5/S6
which were never cached — now read, not cached; /tmp is non-durable).

## V2 — S3 (exit 0)

`read https://arxiv.org/abs/2207.03263` → abs (v4, 11 Nov 2023, math.AP).
`read https://arxiv.org/pdf/2207.03263` → full text, 3089-line extraction
(artifact://35): intro, Theorem 1, reduced systems (1.14)–(1.15),
Kaufmann–Scully ansatz, inner-outer scheme, finite-window limits.

## V3 — S4 (exit 0)

`read https://arxiv.org/abs/2603.21644` → abs (v1, 23 Mar 2026).
`read https://arxiv.org/pdf/2603.21644` → intro + Thm 1.1 + §§1.2–1.3,
6181-line extraction (artifact://40): Hamiltonian/diagonal, Cantor λ set,
degenerate KAM+Nash–Moser, S3-regime adoption, no-stability scope.

## V4 — S5 (exit 0)

`read https://arxiv.org/abs/1405.3469` → abs (v3, 16 Aug 2019).
`read https://arxiv.org/pdf/1405.3469` → full text, 390 lines
(artifact://36): Prop 2 (local converse), Prop 3 (forced Euler), Prop 4
(no-go), Derrick Liouville, helicity/Hopf bounds, R²×S¹ example.

## V5 — S7 companion (exit 0)

`read https://arxiv.org/abs/hep-th/9705176` → abs (v1, 22 May 1997).
`read https://arxiv.org/pdf/hep-th/9705176` → full text, 392 lines
(artifact://38): action (1)–(2), virial (5), Hopf bound (9), gradient-flow
numerics §4–5, G-scale §3. Nature letter (ref [7] therein) not separately
fetched — recorded boundary in s7 note.

## V6 — S6 (exit 0)

`read https://arxiv.org/abs/1909.09054` → abs (v4, 18 Jan 2024).
`read https://arxiv.org/pdf/1909.09054` → full text: S-integrability setup,
Prop 1 (Q=2 analytic), items (i)–(viii) incl. item-iv transport failure,
KKPS comparison, helicity gap.

## V7 — S9 (exit 0)

`read https://arxiv.org/abs/2107.06035` → abs (v2, 4 Feb 2022).
`read https://arxiv.org/pdf/2107.06035` → §§1–2, 1245-line extraction
(artifact://37): Theorems A/B, Corollaries 1.1–1.3, Hill data (2.5),
W_H = 2/15, (L¹∩L²)+impulse class, Pozrikidis confirmation.

## V8 — S8/CLV (exit 0)

`read https://cims.nyu.edu/~vicol/CLV1.pdf` → full text, 529-line
extraction (artifact://39): Grad–Shafranov + localizability (21)–(23),
swirl necessity, Hölder/weak Theorem 2, shell support via cutoff,
multiscale Theorem 3.

## D1–D2 — memory (exits 0)

`memory search "Davila leapfrogging vortex rings 2207"` → no direct
coverage (top hit P026 vortex-plan, 23.64). `memory search "Slobodeanu
Faddeev Skyrme steady Euler correspondence"` → frozen universe + P251
entries only (top 70.22, universe declaration). Both confirm fresh sourcing
was required; no prior claim reused.
