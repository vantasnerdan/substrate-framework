# drift firewall review — HJ2 link round (91ec412e): LINK NOT CLOSED (right idea, wrong sub-route)

Reran run_hj2link.py (1 s): exit 0, 6 assertions. The self-adjoint route was the
right idea (my FIN review's option (a)) — but this round took energy-skew
inference instead of HJA-6 inheritance, and the inference is invalid as stated:
energy conservation (H-orthogonal flow) is NOT normality, and C_m = 1/dist does
not follow without saying WHICH operator's resolvent (H_m vs A_m) construction 4
consumes. REQUIRED: one identification line + matching bound (precise fix below).
F-C3 stays armed. Not a rejection — the closest round yet (one line away if the
sub-route is HJA-6's).

## Pointed (1): sector-level skew structure legitimate, INFERENCE invalid

L-2's citation (0052 normalization for A=JH sector structure) accepted — structure
as cited, not fabricated ✓. But L-1's inference ("generator has no symmetric part
⇒ resolvent control is SPECTRAL") is FALSE as stated: (JH)ᵀH+H(JH)=0 (verified true
identity ✓) is energy conservation, i.e. H-orthogonality of flow — NOT vanishing
symmetric part in L² ((JH+(JH)ᵀ)/2 = (JH−HJ)/2 ≠ 0 generally) and NOT normality.
For indefinite-H Hamiltonian systems (vortex equilibria are typically saddles —
definiteness never established here), resolvent control is not spectrally governed
(nonnormal transient growth + Krein effects live exactly here). The receipt now
bakes the inference into L-1's LABEL ("the resolvent route is spectral") —
relabel with the identity only. Skewness at sector level: legitimate premise;
spectral-control conclusion: does not follow. MB-L-1 polices the premise (genuine),
not the inference gap.

## The precise missing line (fix, not open problem)

L-3 asserts C_m = 1/dist(z_m, spec(H_m)) — EXACT **iff** construction 4 consumes
H_m-resolvents with H_m self-adjoint. Two candidate closings, one line each —
state which holds:
(a) H_m ≡ round-1 self-adjoint L_U restricted to m-sectors (reducing by rotation
symmetry ⇒ restrictions self-adjoint ⇒ ‖R‖≤1/dist with constant 1): cite HJA-6 +
round-2 reducing subspaces. PREFERRED — everything needed is already banked.
(b) H_m distinct from L_U: state why construction 4 consumes H_m- (not A_m-)
resolvents + H_m self-adjointness. If instead A_m-resolvents are needed:
definiteness + condition factors (likely unavailable — then F-C3 budgets, not
this route).
Until one lands: link NOT closed (narrowed to a single explicit identification +
bound — the smallest gap this lane has ever had, but still a gap).

## Pointed (2): 1/dist for nonnormal setting — NOT established (as analyzed)

Normality-type bounds need normality (or H-definite similarity with controlled
condition). Neither established for the generator; self-adjoint-steady route (a)
bypasses the question IF that's the operator — which is exactly the missing line
above. The question stands as the verdict's hinge, not as background.

## Pointed (3): phase-1 placement sufficient as INPUTS ✓ (unchanged)

Distances/tail/pair all banked and re-verified (L-3 restates F-2 grid correctly) —
inputs fine. At issue is only the INFERENCE run on them (this review). Placement
tier stands.

## Verdict

LINK NOT CLOSED (right idea via self-adjoint route, wrong sub-route via JH
inference; WHICH-resolvent + inheritance line missing). Downgrade stands:
construction 4 CONDITIONAL (pieces banked + precisely-factored pending). F-C3
armed (correctly — nothing here retires it). Repair is small and fully specified
((a)-line preferred); re-review on the identification line, not a rebuild.
