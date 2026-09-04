# Attempt 0014 — N2/m=1: interior Poincare wavenumber defect corrected; LIA branch closed

## Diagnosis (implementation defect, route repaired in-run)

The rederivation of the interior problem from source (linearized Euler about
solid rotation, corrected 2 Omega Coriolis, m-general) gives, by two
independent exact routes (direct incompressibility substitution and
div-momentum with curl-momentum):

  P'' + P'/r - m^2 P/r^2 + k^2 (4 Om^2/wt^2 - 1) P = 0,   wt = omega - m Om

i.e. the classical Poincare equation grad^2 p + (4 Om^2/wt^2) p_zz = 0,
which reproduces the plane-wave inertial dispersion wt^2 = 4 Om^2 kz^2/|k|^2
exactly. The campaign's 0008-era lambda-modification carried

  mu^2 = -k^2 (4 Om^2 - wt^2)/(omega wt)   [DEFECTIVE: 1/omega instead of 1/wt^2]

The defect is invisible at the m=2 branch (omega ~ Om there, so 1/omega ~
1/wt^2 at leading order and the I_m vs J_m distinction only enters at
O((ka)^4)), which is why the m=2 channel closed to machine zero through
three probe batteries while m=1 was pathologically unsolvable (mu -> inf as
omega -> 0, manufacturing the 0012-0013 "no exact root" diagnosis, the
0.35k quasi-branch, and the sector-gap enumeration).

Route verdicts on the 0013 candidate sectors: (1) tilt correction of the
sheet-strength transport — EMPTY (all tilt-geometric corrections to the
theta-jump enter at O(eta^2); gamma_z = [v_th] exactly at linear order);
(2) exterior impulse sector — EMPTY (the exterior base flow is the
irrotational Omega a^2/r vortex, omega_0 = 0, so the linearized vorticity
equation imposes no constraint on the exterior potential field; the
campaign's exterior Bernoulli p'_out = i rho wt C K_m(ka) is correct at the
sheet); (3) finite-mu interior — the right neighborhood, but the actual
defect was the interior wavenumber itself.

## Repaired system and result

Corrected interior: P = A I_m(mu r), mu^2 = k^2 (4 Om^2 - wt^2)/wt^2.
The two remaining consistency conditions (pressure continuity with advected
Bernoulli; frozen sheet strength [v_th] = 2 Om eta) satisfy, at leading
small-k order with the exact Bessel expansions,

  r1 = rho a Om · r2          (proportionality verified symbolically and
                               numerically: r1 = r2 to 4 digits at every root)

with the common root

  ** omega = -(Gamma k^2 / 4 pi) ( ln(2/(ka)) - gamma - 1/8 ) + o(1) **

the classical LIA form with log-coefficient exactly 1x and Rankine-core
constant c1 = -gamma - 1/8 [NUMERIC_EVIDENCE, ASYMPTOTIC: root tracked
against the analytic form at k = 1e-2..1e-5; |r| at the root falls to
2.5e-12 at k=1e-5 (vs 6.8e-6 at the wrong-location minima); the extracted
c1(k) = -0.82 ± 0.15 (grid-limited) against the analytic -0.702]. The sign
is the (+theta,+k) helicity; the conjugate helicity gives +.

m=2 control in the corrected system: the branch omega = Om - Om(ka)^2/12
closes at r1 = 4.4e-13 (k=1e-3), scaling ~ k^4 — machine-exact at the
O((ka)^2) order. The recorded m=2 O((ka)^2) results are unchanged (the
defective and corrected interiors coincide through O((ka)^2) at the branch);
only O((ka)^4) structure would differ, which was never claimed. The
single-tube C_tw = 0 verdict (static, k = 0 structure) is unaffected.

## Status

- m=1 bend channel: slow branch CLOSED at leading-log order
  (route verdict: established as stated, with c1 exact to O(1) and the
  o(1) remainder named — the t^{3/2} matched-asymptotics correction).
- The 0012/0013 "no exact common zero", "0.35k quasi-branch
  sqrt(2)/4 Om(ka)", and sector-enumeration records stand as recorded and
  are hereby superseded as defective-input artifacts of the mu^2 defect.
- Next: frozen verifier verify_cst002.py covering BOTH m=1 helicity
  branches and the m=2 channel; mandatory mutations include (m7) the
  defective 1/omega wavenumber must FAIL — it survived three probe
  batteries across attempts 0008-0013, which is exactly the failure mode
  the mutation discipline exists to catch. Then classical corroboration of
  c1(Rankine) from a Kelvin/Saffman-level source (cited-not-recalled).
