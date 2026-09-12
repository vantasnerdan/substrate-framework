# IDEA-05 diagnostic: soft modes vs S9 tail direction (beacon)

Source: cipher IDEA-05 (adopted, owner=beacon): bank soft-mode shapes,
overlap vs S9 diameter mode as diagnostic, never QOI input. External
JD/deflation pointers NOT used (unvetted leads recorded, not consumed).

## Measurement (bank_modes.py, exit 0, trust state, coarse)

8 softest Jacobian modes banked (`soft-modes.npz` + metrics):

| mode | eval | z-dipole | core wt |
|---|---|---|---|
| 0 | ~0.016 | −1.03 | 0.96 |
| 1 | ~0.022 | −1.18 | 0.70 |
| 2–7 | 0.036–0.075 | mixed | 0.49–0.92 |

## Reading

Modes 0–1 are TAIL-elongated (−z dipoles ≈ −1): the solver's soft
cluster CONTAINS the S9 filamentation direction (0114: behind-particle
lag into −z tail). The ×15–60 amplification (0121) acts preferentially
along the exposed channel — the error budget is physically
interpretable, not generic noise.

## Mesh consequence (acted on)

The 0122 fitted band covered the CORE contour only. If the fitted
solve stalls with soft share persisting, the specified next mesh move
is tail-band extension (−z, following modes 0–1), not more core
refinement. Deflation of these modes for QOI purposes is FORBIDDEN
(it would delete the exposed channel from the budget); deflation as a
SOLVER preconditioner is allowed with the unpreconditioned residual
kept as the verdict metric.
