# Primary source and exact evidence receipt

Primary: Enciso and Peralta-Salas, *Existence of knotted vortex tubes in
steady Euler flows*, archived arXiv:1210.6271 version 2, October 2014.

- `../../sources/1210.6271.pdf` SHA256
  `309808ccd801552587a04a13f8cea3a85894511c3387a3f4fa8c5eb9e6d5e738`.
- Extracted text SHA256
  `854f09f74d8806d5dd694c6f72c85eb8a63c3b37d443d8fc331b94dd41f549a1`.

Theorem 6.8, extracted lines 2036-2089: any fixed eigenvalue bound Lambda,
sufficiently small thickness, exact tangent local Beltrami solution and
the weighted coordinate derivative estimate. This licenses the general
leading-jet calculation, not arbitrary fixed-lambda KAM persistence.

Theorem 7.6, lines 2825-2850: analytic area-preserving return map,
Diophantine invariant boundary conjugacy and nonzero normal torsion imply
persistence under sufficiently small same-measure perturbations. The
new torus derives these hypotheses independently. The proof's measure
identification (7.26)-(7.32), lines approximately 3220-3310, handles the
actual nearby velocity's different flux measure before that theorem.

Theorem 8.3, lines 3551 onward: any fixed nonzero eigenvalue, local
Beltrami field on a compact set with connected complement, arbitrary
C^k approximation by a global same-eigenvalue Beltrami field with O(1/r)
derivative decay. It transfers the new robust solid torus. The conclusion
is an unknot; no topology is changed by this approximation.

No use is made of Theorem 7.8's weak-lambda twist coefficient. The exact
axisymmetric Grad-Shafranov equation, Bessel limit, true transit integral,
flux-action twist, parameter selection and gyro/moment comparison are
derived in this attempt, rather than imported from a different MHD or
filament dynamics.

`fixed-lambda-torus.md` SHA256:
`1e9469026859fb062e005fab4b241c290d2f49d668cb5aac7c2bc7deb0bd3193`.
`verify.py` SHA256:
`21c5c107a7cb1f6563ae6554e9faf18bdb3e0f39bce84ea253dd64ee9d8cf5d0`.

First execution exposed a symbolic Bessel quotient-series limitation;
`first-run.txt` preserves its location and classification. Individual
exact Bessel series followed by polynomial division repaired that
representation. The corrected and final strengthened executions both
passed all 20 exact checks. Initial Ruff found one unused intermediate;
the derivative now consumes that intermediate. Final Ruff: all checks
passed. The bounded evidence oracle does not claim to replace the
functional-analytic existence or source-applicability arguments.
