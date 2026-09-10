# Exposing checks receipt — drift baseline (2026-09-10)

Two sympy audit cells on 0107 load-bearing identities. 0107 derivation.py verifier was being
edited at interruption, so these are drift's independent probes, not a re-run of its receipt.

## EC-1: Fourier impulse sign + antisymmetry (0107 §3, eqs 11–12)

Command: eval cell "Fourier sign and impulse checks" (sympy 1.14.0, python 3.12.2).
Probes: (a) `qhat_i = -i k_j M_ji` + `k·qhat=0` forces symmetric part of M to zero;
(b) `M_mn = rho^-1 eps_lmn I_l` is exactly antisymmetric; (c) `qhat = +i rho^-1 k×I`
matches the direct `-i M^T k` expansion (diff = 0 vector) under 0107's convention
`qhat = ∫e^{-ik·x}q dx`, `uhat = i k×qhat/|k|²`;
(d) MUST-FAIL mutation: sign-flipped `qhat` still satisfies `k·q = 0`.
Output: `k^T S k` contraction form confirmed; `M+M^T = 0 matrix`; `qhat-diff = 0`;
`k.q_wrongsign = 0`.
Verdict: PASS for (a–c) — 0107 (12) sign consistent with its stated convention. LIMIT (d):
transversality alone cannot certify the sign; any reviewer check that only tests `k·qhat=0`
is a non-exposing oracle for sign errors. Sign must be checked by direct expansion (as here).
Firewall rule recorded in ledger checklist §4.

## EC-2: impulse kernel identity + Maxwell div-B defect (0107 §2, eqs 4, 8–9)

Command: eval cell "Impulse kernel and Maxwell defect checks" (same env).
Probes: (a) pointwise `(1/2)(x×curl A)` expansion: `eps_ilm x_l eps_mjk d_j A_k`
= `(1/2)(x_k ∂_i A_k − x_j ∂_j A_i)` — confirmed componentwise (comp1–3 match the
δδ−δδ contraction exactly, no leftover terms). Integrated by parts this gives
`(1/2)∫x×curlA = (1/2)(−∫A + 3∫A) = ∫A`, i.e. 0107 (4)/(6) hold ONLY with vanishing
surface terms = the "localized class / expanding-ball flux" hypothesis. That hypothesis is
exactly the unclosed part of Route A (pause-state open point), not a formality.
(b) Maxwell: `(1/μ)(curl B)×B − div σ_EM = −(B/μ) divB` componentwise (all 3 rows:
`lhs−rhs = −B_i divB/μ`). Confirms the pause-state corrected-identity note: the unconstrained
defect `−(B/μ_EM)divB` must be exposed, and the physical momentum equation is recovered only
on `B = curl A`. Any verifier run on a non-divergence-free B without this defect term is a
false green — this was the recorded first-failure mechanism of the local Maxwell oracle.
Output: comp1–3 expansions + diff structure + 3 defect rows, all as expected. Exit 0 both cells.
Verdict: PASS as algebraic probes; they bound but do not close the analytic hypotheses
(decay class, flux limit, fixed-T Euler–Maxwell solution regularity). Those stay on the ledger
as Route-A blockers.
