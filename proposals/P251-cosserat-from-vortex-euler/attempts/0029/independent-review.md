# Independent scientific review: P251 attempts 0029–0034

Reviewer: Codex agent `/root/construction_review`, non-author and
non-implementer of the reviewed constructions. Signed 2026-09-04.
Boundary: `rankine_modes.py`, its tests, revised `verify_cst002.py` and
`verify_cst003.py`, and attempts 0029–0034. Base release: v0.171.0;
branch baseline: 10a0f31. One substantive pass and one correction check.
Review used the actual equations, source, and captured execution outputs;
no tests were rerun and no additional reviewers were commissioned.

## Individual scientific verdicts

| Unit | Strongest supported statement and evidence |
| --- | --- |
| 0029 and `core_velocity` | Established: Cartesian Euler differentiation fixes the original Coriolis signs; 0019's replacement fields fail that operator. The angle/rate dimensional distinction and finite-k contrast-only counterexample are exact. |
| 0031 and revised CST002 | Established: the isolated Rankine mode boundary equation and its branch asymptotics, including the m=1 constant 1/4−EulerGamma and m=2 coefficient −1/6. The recorded 40/60-digit roots corroborate the exact asymptotic derivation with numerical evidence. These are mode results. |
| 0030 | Established: the specified Kida elliptic patch has positive angle stiffness and canonical inertia in its elliptic manifold. Hamiltonian signs, Kirchhoff limit, and frozen-shape mutual moments are consistent. The external linear flow is prescribed. |
| 0032 | Established: the closed six-point-vortex invariant sector has the stated angle action, K=27 rho Gamma²/(8 pi), I=3 pi rho S²/2, and positive oscillator frequency. The independent Cartesian Biot–Savart Jacobian corroborates the normalization. Restored common inertia and cross inertia are correct. |
| 0033 | Established: the tangent-dot-tangent mutual kinetic-energy kernel, its tilt Hessian, and its positive finite-segment contribution for opposite circulations. The statement retains its declared local-interaction scope. |
| 0034 | Established conditionally: the local covariant angle map diagonalizes the complete kinetic form and fixes stiffness from that same action; gradient and isotropic contractions are correct, with L_v=6 L_cell. Macro gradient inertia leaves both transverse omega² branches unchanged through k². The retained mixed gradient terms also first affect higher orders. |
| Revised CST003 | Established conditionally: sphere moments and constitutive coefficient matching for the supplied microscopic energies. Exact Green–Lagrange cancellation correctly prevents deriving angle locking from a free director rotation using line tension alone. |

The exact identities have `symbolic_verified` support; the scalar-root
precision ladder is `numeric_evidence`. These are individual scientific
review verdicts at the stated boundaries, not registry acceptance or a
release promotion. No further scope reduction is requested.

## Correction and correction check

The initial API documentation called the exterior logarithmic derivative a
pressure derivative. With exterior perturbation u'=grad(varphi),
p'=i rho [omega−m Omega a²/r²] varphi. Thus at r=a,

    a p'_r/p' = K + 2 m Omega/sigma,
    K = a varphi_r/varphi, sigma = omega−m Omega.

The implemented determinant and roots already used the correct potential
derivative. The correction check confirms that the module now explicitly
declares irrotational exterior perturbations, the API names the two different
logarithmic derivatives, and the test comment and attempt 0031 specification
agree. No equation changed. This closes the sole requested correction; no
blocking equation finding remains in the reviewed units.

## Original objective and continuation

The strongest new microscopic result is the closed positive relative-angle
action in 0032, with its explicit same-action field map in 0034. It supplies
real angle locking and reaction dynamics. The original exact smooth
stationary Euler-to-Cosserat objective still needs the joint translational
Euler action and affine cell-frame construction, controlled finite-core and
mutual-interaction corrections, and a smooth stationary EPS tube realization
carrying the relevant collective modes. The reviewed point-vortex relative
equilibrium and local-induction gradient model do not supply these dependencies.

Claim promotion: no original C-CST objective claim is promoted by this review.
Goal completion: not established. Scientific exhaustion: not established.
The supported constructions remain substantive campaign progress; the named
missing constructions are the active continuation boundary.
