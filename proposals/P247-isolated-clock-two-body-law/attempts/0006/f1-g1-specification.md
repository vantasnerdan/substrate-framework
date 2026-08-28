# P247 attempt 0006 — F1 refined scope and G1 rotor specification

Preregistered by the attempt-0006 manifest; this document pins the two
deliverables precisely for execution (this run or the next).

## F1 refined scope — source-current extraction and massive pairing

The single-clock measurement (attempt 0005, d1-source.json) showed the
INTEGRATED linear source vanishes exactly. The pairing derivation needs the
source CURRENT density j(x), not its integral: the linear-in-chi action term
has the form

    E_linear = integral j^mu(x) d_mu chi(x) dx,

so the chi equation of motion sources from -d_mu j^mu (a divergence). On a
mu-symmetric single clock, div j = 0 pointwise-or-integrally (parity); for a
parity-broken pair, div j_1, div j_2 != 0 where the clocks' structures
misalign, and the second-order exchange energy is

    E_int(d) = (1/16) integral integral j_1(x) . grad grad G_Lambda(x-x') .
               j_2(x') dx dx',
    G_Lambda(d) = exp(-Lambda d) / (4 pi d),

the massive dipole-dipole pairing (decay e^{-Lambda d} with polynomial
prefactors determined by the multipole content of j).

Extraction plan (numerical, through the verified functional):
1. Nodal-chi probe: relax the modal ansatz for the chi channel to nodal
   values on the (r, mu) grid; the linear response of the total energy to
   each nodal chi value and to its finite-difference gradient gives the
   discrete operator L = d^2E/dchi^2 and the source vector s = dE/dchi at
   chi = 0, directly from autograd on the extended functional.
2. Validate: L reproduces (1/16)(-Laplacian + Lambda^2) on smooth probes up
   to the quadrature budget; s integrates to zero on the symmetric root
   (attempt-0005 consistency).
3. Pairing: place two root copies at separation d (aligned and anti-aligned
   clock axes), evaluate s on each, form the Green pairing above, and fit
   the decay against Lambda.
4. F2: one full two-clock nonlinear evaluation at moderate separation as
   the sign-and-coefficient check.

## G1 rotor specification — fixed-axis clock ansatz

The physical fixed-J clock with an orientation rotates about a FIXED
laboratory axis (say z) rather than the local director:

    S(t, x) = R_z(Omega t) S_0(x) R_z(-Omega t),

breaking spherical symmetry: the relaxed background is axially symmetric
(mu-dependent), so the solver needs

1. angular_modes >= 2 in the S channels (the committed machinery already
   carries the angular_basis and the rotation_z derivative term; the
   extension is the angular-mode count and the fixed-axis response W_z =
   [N_z, S]);
2. the Omega-explicit functional from attempt-0004 (c3_functional already
   takes omega as a parameter) with the fixed-axis response;
3. odd-mu chi modes for the parity-odd boost tail;
4. grids: the committed (r, mu) product quadrature with angular_nodes
   scaled by 2x per added angular mode; aliasing and continuity gates as
   in attempts 0003-0004.

Numerical cost class: 2D modal Newton, comparable to the attempt-0003
driver per Newton step; a full fixed-J ladder is a dedicated attempt.

## Standing gates (unchanged)

F2 confirms or contradicts F1 by sign and coefficient. G1 execution is a
separate preregistration. No PR before issue-178 completes; no fitted
constants.
