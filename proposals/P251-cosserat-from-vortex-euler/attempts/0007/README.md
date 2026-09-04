# Attempt 0007 — N2 block B2 (twist channel), exact mode energies and finite-k identity

## Exact k->0 energies of the m=2 polarization mode (displacement-normalized)

With the sheet displacement amplitude eta and w_t = w - 2 Om:

  v_r^in  = I m A r^{m-1} (Om(m+1)-w)/(rho D)     (inside, P = A r^m)
  phi_out = C r^{-m},  v = grad phi,  C = I w_t eta a^{m+1}/m

  E_in /L  = pi/4 · rho w_t^2 eta^2 a^2        (m-independent prefactor form:
              E_in = pi m |A|^2 Gp^2 a^{2m}/(2 rho D^2), Gp = Om(m+1)-w)
  E_out/L  = pi/4 · rho w_t^2 eta^2 a^2        (E_out = rho pi m |C|^2 a^{-2m}/2)
  E_tot/L  = pi/2 · rho w_t^2 eta^2 a^2        (m=2)

  Effective inertia (coefficient of (d eta/dt)^2 in the mode energy):
    M_eff = 2 E_tot/(w_t^2 eta^2) = (pi/2) rho a^2   per unit length (m=2).

Physical branch w = Om(m-1)/2 = Om/2 => w_t = -3 Om/2:
  K_eff = M_eff w_t^2 = (9 pi/8) rho a^2 Om^2,  Om = G/(2 pi a^2).

Consistency probes: (i) E_in = E_out exactly at k=0 (equi-partition of the
surface mode between core and exterior kinetic energy) -- a nontrivial check
on both amplitude solutions; (ii) dimensional structure rho a^2 matches the
Comparsi microinertia scale J_i = <rho_0 r_i^2>.

## Finite-k outside energy: exact identity, no log in the m=2 channel

Full perturbation energy of the exterior potential field
phi = C K_2(kappa r) e^{2 i th} e^{i k z}:

  E_out(kappa)/L = (rho pi/2) |C|^2 · ( -x0 K_2(x0) K_2'(x0) ),  x0 = kappa a,

via the Bessel identity  int_{x0}^oo x[K_m'^2 + (1 + m^2/x^2) K_m^2] dx
= -x0 K_m(x0) K_m'(x0)  (integration by parts with the modified Bessel ODE;
boundary at infinity exponentially dead; sign positive since K_m' < 0).

Small-x0 expansion: -K_2/K_2' = (x/2)(1 - x^2/4 + O(x^4)), so with the
kinematic normalization |C|^2 = w_t^2 eta^2/(kappa^2 K_2'(kappa a)^2):

  E_out(kappa)/L = (pi rho w_t^2 eta^2 a^2/4) · (1 - kappa^2 a^2/4 + O(kappa^4 a^4)).

STRUCTURAL FINDING: the m=2 (polarization) channel has NO logarithmic
k-divergence -- unlike the m=1 bending channel, whose line-bending energy
carries ln(1/(k a)). The isolated tube's twist stiffness at long wavelength
is therefore carried entirely by the exact O(k^2 a^2) coefficients:

  - outside: -pi rho w_t^2 eta^2 a^4 k^2 / 16  (from the -x^2/4 term above,
    exact, no fitted constant);
  - inside: the Poincare radial wavenumber lambda^2 = -k^2 w_t-dependent
    modification of J_2(lambda r) vs the power law -- O(k^2 a^2), derivation
    owed (attempt 0008).

C_tw = the net coefficient (difference of the two exact O(k^2 a^2)
energies); no log, no fitted cutoff -- the declared long-wave premise
(k a << 1) does the work the cutoff did in the m=1 channel.

## Status

E_tot(k=0), M_eff, K_eff, and the outside finite-k coefficient: established
(SymPy exact + identity route). Inside lambda-modification and the C_tw
assembly: owed to attempt 0008, then C-CST-002's verifier with mutation
probes (wrong D-sign, J_2-in-Rankine mismatch, flipped sheet Doppler).
