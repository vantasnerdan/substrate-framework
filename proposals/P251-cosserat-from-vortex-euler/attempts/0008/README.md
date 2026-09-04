# Attempt 0008 — N2 block B2 (twist channel), inside lambda-piece and the soft-sign finding

## Exact inside O(k^2 a^2) piece

Inside pressure with the Poincare radial wavenumber lambda = Om k / w_t
(w_t = w - 2 Om): P = alpha r^2 (1 - lambda^2 r^2/12) (regular J_2(lambda r)
expansion). Momentum solutions (Cartesian-verified polar basis):

  B_r  = 2 alpha r^2 W + alpha lambda^2 r^4 (Om - 2W)/6
  B_th = -2 alpha r^2 W - alpha lambda^2 r^4 (Om - 2W)/6 + alpha lambda^2 r^4 (2W - Om)/3

with W = -w_t = 3 Om - w (positive on the physical branch). Kinematic BC at
the sheet: B_r(a) = W eta a rho D, D = W^2 - Om^2, giving

  alpha = eta rho D / (2a) · (1 - lambda^2 a^2 (Om - 2W)/(12 W)).

Energy to O(lambda^2): B_r^2 + B_th^2 = 8 alpha^2 r^4 W^2 - 2 lambda^2
alpha^2 r^6 W^2, so

  E_in/L = (pi rho W^2 eta^2 a^2/4) · (1 + lambda^2 a^2/18),
  Delta E_in/L = + pi rho W^2 eta^2 a^4 lambda^2 / 72
               = + pi rho Om^2 eta^2 a^4 k^2 / 72    (lambda^2 = 4 k^2/9).

Outside (attempt 0007): Delta E_out/L = - (9 pi/64) rho Om^2 eta^2 a^4 k^2.

## Sum and the soft-sign finding

  Delta E_tot/L (fixed sheet displacement eta) = pi rho Om^2 eta^2 a^4 k^2
  · (1/72 - 9/64) = - 73 pi rho Om^2 eta^2 a^4 k^2 / 576  < 0.

At fixed displacement amplitude the polarization channel's field energy
DECREASES with wavenumber. A negative k^2 coefficient in a fixed-amplitude
energy is the signature of computing the wrong functional: for Doppler-
shifted modes on a rotating background the kinetic energy is not the wave
energy (indefinite Krein signature); the canonical wave-action energy (or an
explicitly constrained ensemble definition) is the licensed object for any
stiffness claim. This is an analytic-closure question inside N2's license
chain and blocks the C_tw extraction from the isolated tube until the
receipt is earned (governance: analytic closure before production claims).

## Consequences for the ladder

1. C-CST-002's C_tw is NOT claimable from the fixed-eta isolated-tube energy.
   Licensed routes (next attempts): (a) canonical wave-action energy for
   Kelvin modes on solid rotation (Krein-signature receipt), or (b) the
   constrained-ensemble definition of the couple modulus (fixed action per
   length), which likely shifts the couple stress to the ENSEMBLE kernel
   (N3) rather than the single tube -- consistent with 0005's gauge finding
   that an axisymmetric tube carries no twist energy at all.
2. The exact pieces above stand (inside +, outside -, both no-log, no fitted
   constants) and feed whichever licensed definition wins.
3. No modulus claimed. N2 stays active.

## Note

A sign slip in the kinematic boundary condition (advected-sheet Doppler)
was caught by re-deriving the B_r(a) sign from w_t = -W; corrected algebra
re-verified in SymPy (cross-term 2(S T_r + S_th T_th) with the consistent
W-form basis).
