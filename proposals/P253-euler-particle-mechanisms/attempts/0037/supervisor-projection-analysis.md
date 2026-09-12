# Mathematical steering: a bounded physical critical coordinate

Author: recovery supervisor, main-model mathematical analysis, 2026-09-06.
Scope: the existing 0037 constrained-projection route; no accepted claim delta.
Inputs: 0037/README.md, mixed-casimir-flux.md and scattering-block-reduction.md;
0030/column-propagation.md; 0027/solitary-wave-construction.md (uniform positive
vorticity near the axis). This diagnoses a proposed coordinate, not the carrier.
The exact flux identities survive. The P2 topology and positive objective stay fixed.

## 1. An explicit accessible sequence tests the coordinate

Work first at the background column, where the proposed critical projection
must already be bounded. Choose a radial cutoff chi in C_c^infinity([0,infinity))
equal to one near zero, with its rescaled support inside the uniform-vorticity
core omega_0 = W e_z, W>0. Let G be a nonzero smooth compact axial function,
g=G', and F(s)=s^2 chi(s). For epsilon>0 define

    psi_e(r,z) = epsilon^2 F(r/epsilon) g(z),
    v_e = (-psi_e,z/r, 0, psi_e,r/r),
    eta_e = -Delta_* psi_e/r^2,   delta xi = 0.

These are smooth Cartesian, compactly supported, divergence-free velocities.
The background circulation label xi=L(r) is unchanged. Direct change of
variable s=r/epsilon gives the exact physical kinetic norm

    ||v_e||_rho^2 = 2 pi rho [epsilon^2 C_F ||g||_2^2
                              +epsilon^4 D_F ||g'||_2^2],
    C_F = integral_0^infinity F'(s)^2/s ds > 0,
    D_F = integral_0^infinity F(s)^2/s ds > 0.                 (S1)

Every integral is finite. The poloidal part of the positive column norm in
0030 is exactly this kinetic norm. Consequently ||v_e||=O(epsilon).

Writing delta omega_theta = -Delta_* psi_e/r, the axisymmetric curl identity
gives

    M_1(v_e)(z) = integral_0^infinity r eta_e dr
       = [-psi_e,r/r]_0^infinity - partial_zz integral psi_e/r dr
       = 2g(z) - epsilon^2 B_F g''(z),
    B_F = integral_0^infinity F(s)/s ds.                     (S2)

Thus M_1(v_e) tends to 2g in axial L2 while the input tends to zero.
Equivalently M_1(v_e/||v_e||) is unbounded. Nonzero evaluation on the critical
eigenvector does not make M_1 a continuous functional on the full state space.

## 2. Accessibility, fixed invariants and the generator domain survive

The counterexample is not an arbitrary nonaccessible velocity. Define the
compact smooth azimuthal displacement

    a_e = (1/W) partial_z^{-1}(delta omega_theta) e_theta,

where the axial primitive is explicit using g=G':

    a_e,theta = -(1/(W r)) [
          (partial_rr - r^-1 partial_r)(epsilon^2 F(r/epsilon)) G
          + epsilon^2 F(r/epsilon) G''] .                  (S3)

Near the axis F(r/epsilon)=r^2/epsilon^2, so its first bracketed term vanishes
and a_e,theta is proportional to r; Cartesian smoothness follows. Pure
axisymmetric azimuthal a_e is divergence-free. On its support omega_0=W e_z,
and

    curl(a_e cross omega_0) = W partial_z a_e
                            = delta omega_theta e_theta.  (S4)

Outside this support both sides vanish. The velocity reconstructed from this
vorticity is v_e itself, by whole-space finite-energy uniqueness. Each state
is therefore an exact compact-generator tangent to the same coadjoint orbit.

Every integrated mixed-Casimir first variation vanishes: eta_e is a linear
combination of g and g'', both with zero axial integral, and delta xi=0.
Pure-label Casimirs have zero first variation as well. Fixing those invariants
does not remove (S1)-(S2).

On these pure poloidal states the translating column generator acts as

    G_c v_e = c partial_z v_e - W v_e^r e_theta.             (S5)

The azimuthal term is divergence-free, so there is no omitted pressure term
in (S5). The velocity and its derivatives are supported in the uniform core.
Equation (S1) gives ||partial_z v_e||=O(epsilon) and
||v_e^r||=O(epsilon^2); hence ||v_e||+||G_c v_e||=O(epsilon).
The graph norm therefore does not cure the discontinuity. This is not an
objection to a stronger Sobolev topology when explicitly justified; it is a
counterexample to the presently proposed energy/column-graph projection.

## 3. Positive repair: a regular-label density, then its spectral dual

Choose a smooth nonnegative nonzero D supported strictly inside the regular
label range, away from the axis and flat exterior label values. On the column,

    M_D = integral r eta D(L(r)) dr
        = integral (psi_r/r) partial_r[D(L)] dr
            -partial_zz integral psi D(L)/r dr.             (S6)

Both radial weights have compact support in an annulus. Cauchy-Schwarz and
radial energy control bound the first term; with a fixed low-axial-frequency
cutoff the second is also bounded in the physical energy norm. On the
threshold profile f_0>0,

    eta_0 = Phi f_0/(c_0^2 r^2),
    M_D(f_0) = integral Phi f_0 D(L)/(c_0^2 r) dr > 0.       (S7)

Thus there are bounded, nontrivial critical coordinates available in the
original energy setting. The local mixed-Casimir flux law applies to this D
and retains its output axial derivative. A regular-label coordinate avoids
the artificial axis trace; it does not discard physical perturbations.

Two constructions still matter. First, this coordinate is not automatically
the left eigenvector of the physical column oscillator, and a poloidal
density alone need not distinguish its positive and negative traveling
branches. Construct the bounded transformation to the actual adjoint
spectral pair, with any flux companion and free off-diagonal terms retained.
Second, an output derivative is not by itself a uniform small parameter:
prove the scaled operator/local-energy estimate in the stated norms. A mere
factor ik does not automatically cancel the separate (c-c_0)^-1 time loss.

## Route verdict and continuation

The D=1 energy-bounded coordinate route is refuted by (S1)-(S5). The exact
local conservation law is preserved. Equations (S6)-(S7) supply a concrete
repair on the same carrier and physical perturbation space; the required
spectral-dual transformation and scaled complement bound remain active.
This conclusion neither refutes the solitary wave nor reduces the P2 goal.
