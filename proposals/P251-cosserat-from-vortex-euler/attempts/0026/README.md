# Attempt 0026 — Euler inter-tube interaction: exact crossing energy, parity verdict, and the polarization-locking specification

## Review context

PR #199 review (changes requested) correctly identified that the N3 coupling
`alpha = L_v T/6` was produced by retaining the quadratic norm of the
FIRST-ORDER relative rotation (`grad u - skew(Phi)`), while the exact
relative Green-Lagrange measure cancels for any free rigid rotation pair
(`R^T R = I`). Straight-tube tension alone therefore supplies no
frame-locking interaction; an Euler-microscopic interaction with nonzero
EXACT second variation in the relative rotation is required. This attempt
works the reviewer's two named candidates.

## Part 1 — inter-tube Biot-Savart crossing route: exact receipt, then BLOCKED with mechanism

Two skew straight filaments, circulations Gamma, minimum separation d,
relative angle theta (line 1 along zhat; line 2 through d*yhat along
(sin th, 0, cos th), offset perpendicular to t2). The Biot-Savart cross
energy E_int = rho Gamma^2 ∮_{C2} u1 · dx2 has CONSTANT numerator
(t1 x R)·t2 = −d sin th, so the double integral reduces to one quadrature:

    E_int(theta)/rhoGamma^2 = (−d sin th / 4π) ∫∫ [d² + r² + s² − 2 d r sin th]^(−3/2) ds dr
                            = −(1/2) tan(theta)          [EXACT]

Symbolic reduction + 30-digit mpmath quadrature agree to ≤ 2e-19 relative
(`parts/crossing_energy.py`, exit 0).

Properties:
- **Objective** (structural): E_int depends only on relative geometry
  (theta, d); a coherent rigid rotation of both filaments leaves it exactly
  invariant. Objectivity is not the failure mode here.
- **Parity-ODD**: E_int(−theta) = −E_int(theta); second variation
  E'' = −rho Gamma² sec²theta · tan theta is nonzero at generic theta but is
  itself odd.
- **Verdict: BLOCKED (with mechanism) for the quadratic locking sector.**
  A parity-symmetric isotropic ensemble has <E''> = 0: the crossing
  interaction contributes NO quadratic |rot u − 2 Phi|² locking stiffness.
  The surviving structure is a first-order chiral (parity-odd, Hall-type)
  coupling, which the isotropic ensemble kills at mean level; recorded as a
  frontier observation for chiral sub-ensembles, not a micropolar modulus.

## Part 2 — core-polarization route: specification (to be executed in 0027/0028)

The gauge argument does not extend to the m=2 core polarization: it is a
shape degree of freedom locked to the tube's circulation (Omega_self =
Gamma/(2 pi a²)), not a frame choice. The candidate exact locking energy is
the frozen-vorticity (Cauchy-invariant) energy of a displaced core in an
AMBIENT rotation: in the core frame the ambient vorticity is
2(Omega_amb − Omega_self), the frozen sheet strength is
gamma' = 2 (Omega_amb − Omega_self) eta, and 0011's exact Poincaré
streamfunction solution gives

    E_lock/L = (pi rho / 2) (Omega_amb − Omega_self)² a² eta²   [candidate]

— exactly quadratic in the RELATIVE rotation only, vanishing under coherent
rotation (objectivity satisfied by construction, not by truncation), with
stiffness ∂²E/∂(Omega_amb − Omega_self)² = π rho a² eta² > 0.

Required before claiming (0027): re-derive the two-frame frozen-sheet jump
and the streamfunction energy with the momentum-residual-correct velocity
forms of attempt 0019 (`parts/mode_solution.py`:
v_r = −i(w_t p′ + 2 Om m p/r)/(rho D), v_t = (2 Om p′ + m w_t p/r)/(rho D),
v_z = k p/(rho w_t), D = w_t² − 4 Om²), with mutation probes on the sheet
jump sign and the relative-rotation structure. Then 0028 rebuilds N3 as
single-tube (objective, exact measure) + pairwise/ambient interaction
locking, and replays the ensemble contraction.

## Status

- route: crossing-route BLOCKED (exact, mutation-ready via the parity
  argument); polarization route SPECIFIED with exact candidate structure.
- evidence_scope: EXACT (sympy + 30-digit mpmath) for part 1; candidate
  (declared, not yet earned) for part 2.
