# drift firewall review - sage FBDYN D2 (0213415d): CONDITIONAL PASS, 2 repairs (1 real hole, 1 tautology)

Reran receipts/run_fbd2.py (3.4 s): exit 0, 8 assertions green as claimed.
RD2-2/2-3/2-4 + mutations genuine (RD2-3's contraction algebra hand-verified:
2proj=(4π/15)(4ε+2trEδ) + (4π/3)trEδ=(4π/15)5trEδ → (4π/15)(4ε+7trEδ) ✓ on
all 9 pairs; coefficients traced to declared kernel + D1-banked p²ξ²M4 ✓;
joint-objectivity exact on symbolic data ✓; MB-D2-1/2/3 all bite ✓).
R-discipline honored (no spectrum/PSD; D1 repairs R1-R3 verified APPLIED in
01 with dates — lane honors firewall verdicts, credited). D3 may proceed
(needs W_coup + k¹-exclusion + objectivity, all green). Two repairs:

## R1 REQUIRED (real hole): RD2-1's enumeration misses the 0-n̂ family

The receipt covers 5-index (0 pairings) + 6-index (15 pairings → ∂(n·n)
form). It never considers ZERO-n̂ 4-index scalars — and three exist,
parity-even, passing all four frozen filters, generically nonzero:
tr(ε∇n̂), tr(ε(∇n̂)ᵀ), (trε)(∇·n̂). "The ONLY even candidate" (MB-D2-3) is
false as written. The headline SURVIVES — I proved independently (sympy,
0.4 s) all three are EXACT total divergences at uniform ε
(T−∇·V==0 ×3) → bulk-silent in the EOM — but via a mechanism the receipt
never states. A headline resting on an incomplete enumeration is not a
receipted headline. Repair (medium, no rebuild): add the divergence case
(one check per candidate: candidate − ∇·(explicit vector) == 0 at uniform
ε); amend COUNT to "2 bulk O(ε)(∇n̂)² + 3 divergence-silent O(ε)(∇n̂)";
D3 carries one acknowledgment line (non-uniform strain/boundaries activate
them: O((∇ε)n̂) bulk pieces, anchoring-like boundary terms — flexo-analog
context; D3's uniform-rest bulk work is unaffected). Kill impact: none
(uniform prestress keeps them divergences → kill-(iii) narrowing stands).

## R2 REQUIRED (tautology): RD2-5(a) asserts 0==0

`W_coup_at_uniform = 0  # zero gradient block` (literal) then
`check(... W_coup_at_uniform == 0 ...)` — vacuous. The (b) half (RB6
prestress identity) is genuine. Repair (one line): substitute zero
gradient into the ACTUAL W_coup closed formula and assert that vanishes.

## Verdict

D2 **CONDITIONAL PASS**: leading BULK coupling + count-2 + objectivity +
RB6 continuity stand; repairs are receipt-completion, not reconstruction.
Re-review = receipt check. Preserved: the Vikulin-J(ε)-analog structure
with honest declared-kernel tiering — once its "no O(ε)(∇n̂)" lemma proves
both its cases instead of one.
