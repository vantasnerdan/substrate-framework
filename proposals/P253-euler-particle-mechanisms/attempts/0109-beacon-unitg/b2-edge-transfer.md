# B2 edge-transfer estimate — cipher PoC-1 denominator (beacon 0109)

Request: cipher inbox 2026-09-10T16:01Z — accessible-packet/interface transfer
magnitude + scaling in R/Z at the compact edge; order-of-magnitude + limits
suffice for a capped EXPLORATORY verdict. Algebra pinned by 0109 checks B1–B6.

## Frozen-column exact mechanism (0107 (51c)–(51d), re-verified)

For `C_col = x[[0,−Z],[R,0]]` with `R,Z > 0`, `ν = |x|√(RZ)`:

- `exp(t·C_col) = cos(νt)·I + sin(νt)/ν·C_col` (B1+B2).
- At `t_s = π/(2ν)`: lower entry magnitude `√(R/Z)` exactly (B4).
- Symmetrizer `H_0 = diag(R,Z)`: `H_0·C_col` skew (B5) — uniform
  symmetrizability through `x = 0`; the physical KKS form `|x|·H_0`
  degenerates there, so this is NOT coercivity.
- Center: `K = 4ZR/Ω² → 16` at smooth center (B6); tongue coefficient
  `|b(0)| = 15|U_θ|/256 > 0` (0104 Unit A, cited).

## Edge scaling (fixed `R → R_edge > 0`, `Z → 0`)

| Quantity | Scaling | Limit |
|---|---|---|
| Amplification at `t_s` | `√(R_edge/Z)` → ∞ | No uniform-in-core all-time bound in unweighted amplitude norm |
| Arrival time `t_s` | `π/(2|x|√(R_edge·Z))` → ∞ | Slow channel; fast packets see less |
| At `Z = 0` exactly | Jordan shear `I + t·x[[0,0],[R,0]]`, growth rate `|x|·R_edge·t` (B3) | Linear-in-t, not exponential; nilpotent |

## Hard limits (why EXPLORATORY cap stays)

1. Frozen-column ≠ finite-Cao DA transfer: no accessible localized edge packet
   or interface/free-boundary operator control exists (0104 Unit I blocked).
2. The accessibility weight may remove or renormalize this channel (B2 active).
3. `x → 0` needs graph/resolvent control; KKS form degenerates.
4. PoC-1 normalization use: divide reconnection-barrier numerator by `√(R/Z)`
   evaluated at the PoC's edge parameters; report the ratio as
   order-of-magnitude with (1)–(3) as stated caveats. Not a stability claim.

## Pointers

Verifier: `verify_unitg_b2.py` checks B1–B6 (exit 0, see `tool-receipts.md` U1).
Source derivation: `attempts/0107/derivation.md` §§10–11 (draft, cited scope only).
