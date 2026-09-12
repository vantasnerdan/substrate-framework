# Bounded review correction for 0041

The independent `0041` review identified domain, wording, and convention
overstatements without changing the exact threshold, flux, ray, or abstract
block formulas. This correction makes the minimum scoped repairs:

1. `threshold-reduction.md` now states that the steady branch and an assumed
   local conservative KdV form fix `beta` by compatibility. The
   time-dependent full-Euler adjoint projection, pairing, and complement
   normalization remain open.
2. `ray-exclusion.md` retains strict axial no-return, the integrable one-pass
   strain bound, and zero high-frequency ray exponent, while withholding a
   literal scattering-matrix limit until a polynomially conjugated or
   weighted integrability estimate is proved.
3. `block-scattering-lemma.md` now declares the closed-generator/Kato
   resolvent realization and two-sided inhomogeneous supersmooth hypotheses.
   The P-channel supersmooth estimate is an open scalar resolvent step. The
   document no longer infers exclusion of full-Euler high-frequency poles
   from the principal-ray exponent. The Euler factorization, high-frequency
   operator exclusion, physical adjoint conversion, and scaled coupling
   estimate stay explicit construction targets.
4. Both block documents use `A_P,0=+sigma partial_y L_*` in the declared
   right-moving coordinate `y=X-sigma T`. This is the sign obtained by
   linearizing the displayed KdV equation about `A_*(y)`.

No equation-level oracle changed, so no unchanged verifier was rerun.

Pre-correction SHA-256:

- `threshold-reduction.md`: `0cfe43b2138ce910940189c0ecc56648aca6625ae841baf6150ffe6f93edf900`
- `ray-exclusion.md`: `f7fa76057492c242e2e211aff6fd7036654d0af7bbb8a9d544e303a559f8fc56`
- `block-scattering-lemma.md`: `e68ff5c3213c412e82fc78afc5861481c709c2771a76a5ee262e1890f7c00195`
- `scattering-block-reduction.md`: `994e160d8c785f338b0502e1a8f63dc32c597b2259fc8e6e44402499ec0e73ba`

Post-correction hashes are appended after the scoped checks.

Post-correction SHA-256:

- `threshold-reduction.md`: `5a2b619df9ec69f3f9865a6efcad11c9090d6b3ce0bab543051fa1c07df13a23`
- `ray-exclusion.md`: `fe3a0aa0b2e1c2070ad6848b9d5cbe05de49e188e4f71919271e04254819fc52`
- `block-scattering-lemma.md`: `ccc65503594b67ad88c0aa9641f67995f3288b48bd5b9d5ac4c4b050e7100d90`
- `scattering-block-reduction.md`: `35fbb57dc23435fcf7158249c4041a86ca5173d7227c365f796b8120cbb9f5ac`

The scoped `git diff --check` completed with exit zero.
