# 0132 verdict: S1-3D KILL (sign + m=1 residual) (beacon)

Build: run_3d.py under frozen 0131 bars, no loosening. A3 code
imported read-only. Gate first: eps-linearity PASS (2.63e-4).

## Measured (printed, exit 0)

- Shared kappa* = -6.6596 → WRONG-WAY SIGN (L = +1 seed) KILL.
- m=1: r = 0.9047 → KILL (> 0.25; 3.6x above line, 9x above pass —
  nowhere near gray zone).
- m=2: r = 1.1158 on response norm 3.9e-5 (500x below m=1's 2.0e-2)
  → leg INCONCLUSIVE (fits noise, not evidence; reported
  non-firing, contributes nothing to the kill).
- F1 verdict: KILL (sign + m=1 residual, either alone sufficient).

## Mechanism (named)

Plane-mismatch persists into m = 1: residual ≈ 0.9 repeats the
non-overlap signature, and the best-fit direction is opposite
(negative kappa) — stronger than the axisymmetric kill (which kept
slope agreement 1.05). Breaking axisymmetry did not open a
fittable toroidal channel: the Magnus-form template remains
structurally orthogonal to the response, with wrong-way best fit.
m ≥ 3 response norms ≤ 1e-7 (no content at this seed).
## Confirmation alignment (drift review-beacon-s1-3d, KILL CONFIRMED)

Over-determined: m=1-alone own-best fit reproduces pooled kappa*
exactly (pooling-proof sign); residual at own-best-fit (stronger
than shared-kappa reading). Discipline adopted going forward:
exclude-empty-channels-from-pool (moot here by 1000x weight).
Surviving addresses (unbuilt, none claimed): m>=2-seeded charge
tests; non-leapfrog carriers; non-Magnus filament back-reaction
forms. Full-v×B template correction carries to any future filament
build.

## Scope (no overreach)

Kills S1-3D as built (leapfrog + m = 1 seed, Magnus t̂×mutual,
two passages). Does not touch S3/S4 lanes. S1 native program:
axisymmetric dead + 3D dead → filament-charge route CLOSED
pending a materially different charge construction (none on the
table; not claimed).
