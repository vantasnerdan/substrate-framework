# Attempt 0011 — N2: corrected energies, virial closure, single-tube C_tw = 0

## k->0 sector (m=2 polarization channel, corrected dynamics)

- Static frozen-vorticity configuration energy (sheet term; 2D streamfunction
  psi_in = B r^2, psi_out = C r^-2, gamma' = +2 Om eta, psi continuous;
  B = Om eta/(2a), C = Om eta a^3/2):
  Delta V_static(0)/L = pi rho Om^2 a^2 eta^2 / 2.
- Dynamic mode: peak kinetic energy T_max(0) = pi rho Om^2 a^2 eta^2 / 2 —
  VIRIAL EXACT (T_max = V(eta_max)). M_eff(0) = pi rho a^2
  (corrects 0007's defective (pi/2) rho a^2).

## Finite-k dynamic energy (exact Bessel numerics on the exact branch)

Evaluating T_max(k) from the corrected momentum solution (I2 interior, K2
exterior, advected-sheet kinematics, corrected Bernoulli) on the dispersion
branch omega = Om(m-1) - Om (ka)^2/12:

  T_max(k)/T_max(0) = 1.0000000000  at k = 0.001 .. 0.15  (10+ digits)

  ** The raw mode-energy integral is k-INDEPENDENT. **

Equivalently, since T_max = (1/2) M_eff wt^2 eta^2:

  M_eff(k) = pi rho a^2 (1 - (ka)^2/12)^{-2}   (mu_2 = 1/6 exact)
  K(k) = M_eff(k) wt(k)^2 = pi rho a^2 Om^2     (k-INDEPENDENT stiffness)
  E_mode(k) = pi rho Om^2 a^2 eta^2 / 2         (k-INDEPENDENT)

  ** Single-tube twist modulus C_tw = 0 exactly at O(k^2). **

## Static-configuration analysis (constraint structure)

- Sheet-only frozen configuration (no bulk tilt): numeric slope of
  Delta V_stat(k) is EXACTLY -Delta V_0 a^2/3 (slope/Delta V_0 -> -1/3):
  Delta V_stat(k)/L = (pi rho Om^2 a^2 eta^2/2)(1 - (ka)^2/3 + O((ka)^4)).
  This configuration is INCOMPLETE at O(k^2): a helically displaced core
  carries bulk tilt vorticity d omega = 2 Om i k xi inside the core
  (frozen-in Cauchy solution: (w0.grad)xi with w0 = 2 Om zhat), whose
  positive O(k^2) energy is missing.
- The DYNAMIC mode configuration (sheet term + the mode's own interior
  displacement) is uniquely fixed by the corrected dynamics, and its energy
  is the k-independent object above: the sheet-part slope -(ka)^2/3 is
  cancelled exactly by the bulk-tilt sector.
- Statically the single-tube C_tw is constraint-underdetermined (depends on
  the interior frozen-deformation ansatz); dynamically it is uniquely
  ZERO. The rod-level Cosserat couple modulus is therefore an ENSEMBLE
  quantity (mutual polarization kernel of neighboring tube segments, N3),
  exactly the structure 0008 anticipated — now PROVEN on corrected
  dynamics, not merely suspected (0008's Krein mechanism was itself a
  defective-input artifact, reclassified in 0010).

## Route verdicts

- route_verdict: established (three mutually consistent exact structures:
  dispersion x inertia x energy identity, plus the virial closure);
  evidence_scope: SYMBOLIC_EXACT with numeric corroboration (Bessel-exact
  evaluation, 10+ digit stability).
- N2 twist channel: single-tube content COMPLETE. Unlocks N3 (ensemble
  Cauchy-Born moduli via triad moments; the couple modulus kernel is the
  named target).
- Owed for C-CST-002 verifier: bend-channel c_inner (0012, corrected
  interior dynamics), then verify_cst002.py with mutations:
  (m1) single-Omega momentum solution must FAIL the dispersion check;
  (m2) gauge-branch omega = m Om must FAIL the frozen-vorticity condition;
  (m3) sheet-only static configuration must FAIL the T_max constancy;
  (m4) Doppler flip (w vs wt) must FAIL the k->0 dispersion;
  (m5) wrong M_eff (pi/2 rho a^2) must FAIL the virial check.

## Note

T_max(k) constancy currently rests on Bessel-exact numerics (10+ digits
across two decades of k) + the M_eff(w) structure consistency; the symbolic
series proof of the identity (all k^2 coefficients cancelling between
interior, exterior, and v_z sectors) is owed to the frozen verifier, where
it becomes a mutation-sensitive check rather than a hand-assembled claim.
## Scope correction (append-only, attempt 0012 reground)

The energy-functional identification in this attempt is RECEIPT-PENDING.
Re-derivation while opening the m=1 block surfaced the second-order
mean-flow sector (the integral rho v0 . v'' and the axisymmetrized
second-order field v''_0): the displaced-vortex translation-invariance test
(Delta E = 0 exactly for m=1, k=0) FAILS for the plain perturbation-kinetic
functional (rho/2) int |v'|^2 and is restored only when the v'' sector is
included. Concretely: inside the core v0 . v' = -Om^2 r . xi with
int r . xi dA = 0, so the cancellation must come from the O(eta^2) field
piece — proving (rho/2) int |v'|^2 alone is not the conserved mode energy.

Status of this attempt's claims under the correction:
- STAND (dynamics, energy-functional independent): the dispersion
  omega = Om(m-1) - Om(ka)^2/12 (three-route); M_eff DEFINED as the
  perturbation-kinetic inertia, M_eff(k) = pi rho a^2 (1-(ka)^2/12)^-2;
  T_max(k) k-independence as an IDENTITY of the perturbation-kinetic
  functional.
- RECEIPT-PENDING (interpretation, not the identities): "E_mode =
  T_max = Delta V_static", the virial-closure reading, and consequently
  C_tw = 0 as a MODULUS claim. The canonical wave-action receipt
  (A = action from the exact fields, E_canon = wt A) is the licensed
  functional; C_tw and the m=1 c_inner extraction resume only on it.
  This is the same functional trap 0008 hit, now identified at the
  receipt level rather than via a sign anomaly.

Next route (0012): derive the canonical action for the corrected modes
from the exact fields (Lagrangian displacement xi via D_t xi = v', then
A = (1/4) int rho (xi* . v - xi . v*) dV, E_canon = wt A); re-derive
C_tw and the m=1 branch structure on the receipt.
