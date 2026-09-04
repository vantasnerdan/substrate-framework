# Attempt 0006 — N2 block B2 (twist stiffness), m=2 polarization-mode dispersion

## Route

Exact linear Kelvin mode of a Rankine vortex tube (vorticity 2Om = G/(pi a^2)
inside r<a, potential outside, vortex-sheet boundary), k -> 0, m-harmonic.

Representation discipline (per problem-deconstruction): after two polar-
coordinate eliminations failed on curvature/contraction bookkeeping
(representation verdicts, not physics), the linearization was redone in
Cartesian form where (u0.grad) = Om d/dth and (V.grad)u0 = Om zhat x V are
exact vector identities; polar components introduced only through the
orthonormal basis, with SymPy carrying every factor.

## Result

Inside momentum solution: v_r = I(Om m r P' + Om m P - w r P')/(r rho D),
v_th = (-Om m^2 P - Om r P' + m w P)/(r rho D), D = (w - m Om)^2 - Om^2.
Incompressibility factors as I (Om m - w)^3 (r^2 P'' + r P' - m^2 P):
regular inside solution P_in = A r^m (Laplace), outside P_out = B r^{-m}.

Boundary conditions (sheet advected by the base flow -- kinematic condition
uses the material derivative w_t = w - m Om, NOT w):
  kinematic in : I m A a^{m-1} (Om(m+1) - w)/(rho D) = -I w_t eta
  kinematic out: -m C a^{-m-1} = -I w_t eta ;  p_out = I w rho C a^{-m}
  dynamic      : p_in(a) = p_out(a)

Dispersion:  A from kinematic-in = -w_t eta rho D / (m a^{m-1} (Om(m+1)-w));
A from dynamic+kinematic-out = -w w_t rho a^{1-m} eta / m.  Equating:

  (w - Om(m+1))(w - Om(m-1)) = -w (Om(m+1) - w)   for w_t != 0
  =>  2 w = Om (m - 1)   =>   w = Om (m-1)/2,
  plus the w_t = 0 branch w = m Om (both amplitudes A = C = 0 there).

  ** omega in { m Om ,  Om (m-1)/2 } **

## Theorem probes (both pass)

1. m=1 translation neutrality: the (m-1)/2 branch gives w = 0 at m = 1 --
   a transversely displaced straight vortex is an equilibrium. PASS.
2. Kirchhoff circular limit: the w = m Om branch is the sheet co-rotating
   with the base flow (zero perturbation fields, zero energy: at w_t = 0
   both amplitude solutions give A = C = 0): pattern speed w/m = Om --
   exactly the epsilon->0 limit of the Kirchhoff elliptical-patch rotation
   rate, i.e. the polarization is an equilibrium OF THE ROTATING SHEET to
   linear order, which is Kirchhoff's O(epsilon) result. PASS.

## Status and next route

C-CST-002 extraction continues on this branch: the physical polarization
wave is w = Om(m-1)/2 (the w = m Om branch is the kinematically advected
sheet displacement, zero perturbation energy); the mode energy's second
variation with respect to the polarization amplitude gives the twist
modulus C_tw, and the mode's kinetic-energy moment gives the microinertia
J_i -- both against the Comparsi target form (issue #198, eqs 1-2).
Two-route check owed: the static-ellipse (quadrupole) energy from the
line-tension functional must reproduce the mode's potential-energy
curvature. No modulus is claimed until that identification and the energy
two-route land (attempts 0007+).
