# Independent closed-core Kelvin construction review

Reviewer root; frozen base v0.181.0 /4bae709. The prior straight-lift
suggestion is disclosed in README.md. Root did not implement the compact
return, Piola map or global pressure transfer reviewed here. Central
registration and schema activation (269 accepted claims) preceded this
source-body review. No accepted statement is changed by this record.

## Strongest supported positive result

0222 establishes its stated actual field-changing Kelvin/Euler response
on the SAME0211 smooth stationary ring with a literal nonzero closed
Beltrami core. Its positive fixed material tag is moved by that same
Kelvin displacement. The physical phase and energy have positive margins;
the initial spin/phase inertia is normalized exactly by a positive tag
fraction. Actual finite-time angle, G and spin approach the common
positive oscillator normalization in the stated ordered limits, with
the full pressure and the lowest ring harmonic retained.

The statement is accepted as a research construction without correction.
It does not claim a spatial continuum or a global spectral pole, and
neither is needed for this finite-time statement. The parent remains
active through its actual common-K, density and acoustic joining work.

## Load-bearing proof audit

The compact inverse is licensed by an actual scalar source moment, not
by a radial wall. For the m=1 radial Laplacian,
`s^2 Delta_1 f=(s^2 f'-s f)'`; regularity at zero and the canceled source
moment remove the only decaying exterior 1/s coefficient. The return
is outside the observed tag but inside the literal-curl region. It has
signed amplitude, not signed material density. Smooth flat matching
makes both the potential and its divergence completion compact.

The right-handed tube map has positive Jacobian `1+X/R`. Its Piola
pushforward preserves the exact divergence equation, and the integer
carrier kR=-1 is single-valued. The full initial Kelvin velocity is the
whole-space projection of the actual compact force. The literal lambda
relation is used only on its support. Subsequent fields solve Euler/Lin;
the compact inverse moment is not incorrectly required to survive
radial dephasing at later times.

The pressure bound is genuinely global. The Newton kernel integrated
over a fixed-width tube about a radius-R circle costs O(1+log R) per
source point, including the near cross-section and far arc. The squared
amplitude O(R^-2) and tube volume O(R) therefore give
`||div W_app||dotH^-1^2=O(log R/R)`. Derivative versions follow from the
same support and derivative bounds. Relative to the actual O(sqrt R)
field norm, the Leray correction is O(sqrt(log R)/R), including its
noncompact exterior. This is not the high-angular-mode estimate from
an unrelated route.

The actual background derivative bounds are uniform globally: the local
singularity uses smooth vorticity cancellation, and the far Biot-Savart
kernel is integrable against the circular tube's linear volume-growth
bound. The residual after projection retains the transport of its
nonlocal pressure correction. Higher derivative norms can be selected
before applying the fixed-time Euler energy estimate. Duhamel followed
by Lin reconstruction thus gives the stated normalized O(log R/R)
transfer with constants depending on the already fixed tag/core scale.

The complete cotangent, not vorticity phase alone, determines the sign.
In the straight column the axial contribution to `w+Du xi` cancels.
Integrating the remaining phase radially gives
`beta=rho*pi*L*m*integral Z' A^2`, positive for m=-1 and Z'<0. The full
Jacobi energy reduces after radial integration to
`rho*pi*L*lambda^2*m^2*integral s*Omega^2*A^2/2`.
Its Hessian is consequently beta*Omega_0 in the small-core limit, not
twice that value. The return alters the actual positive energy and
phase integrals and is retained in their normalization.

The physical moment tests keep the moving centroid, symmetric shape,
and the extra toroidal G/spin terms. The error estimates are observation
row bounds rather than division by a sinusoid at its zeros. Choosing
the small core first and then R makes the full transfer errors small
relative to the nonzero phase and tag signal. The reversing half-turn
indeed takes the reflection-symmetric axisymmetric ring to its time
reversal, so the stated real initial parity registration is available.
The exact positive initial ratio can then scale a fraction below one;
neither finite-time equality nor a nonzero density limit is silently
inferred from that initial match.

## Evidence and validation boundary

Read the entire 272-line proof, 98-line verifier, first output and receipt.
The 16 existing exact checks expose the full cotangent, energy factor,
compact moment, divergence completion and nonzero curvature terms.
The analytic global estimate is the transfer oracle; its scale checks
are corroborating anchors, not a simulation or formal PDE proof.
Source hashes match the frozen receipt exactly:

- Proof: `9b0ce9340cb58988778944d67fb2787aa14261852dd832259def4e1ab6e40175`.
- Verifier: `6e666fdf21db89973c99b7f35ce576fa6d0c4f4c446aa0acebc1f33e122e9051`.
- First output: `e70e5f7cb56dcdaa42a102fef7033098728f1fdd21297f63e20ff57dc585f07a`.

The banked first execution is 16/16, exit0; Ruff and source diff checks
are unchanged. No duplicate oracle or full-suite replay was needed.
Only this new evidence delta and its immediate mathematical inputs were
reviewed; no canonical API, old claim or generated document changed.

## Decision and continuation

Route verdict: established as stated. Evidence scope: exact compact
Kelvin construction with an analytic controlled finite-time full-pressure
response and positive physical action/current normalization. Review:
accepted at this construction scope. Compatibility: compatible research
extension. Epistemic: active parent. Correction check: not needed.

The next parent construction is the common spatial/current and acoustic
join on the same stationary ensemble; a commit does not close it.
