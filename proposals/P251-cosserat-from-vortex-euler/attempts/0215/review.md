# Independent review: actual constant-curl core on the global Euler ring

Reviewer `/root`; base v0.181.0, checkpoint 61ac2ea. The candidate-
suggestion role is disclosed in README.md. This review concerns the new
implemented seed, complete border, global continuation and physical return
map. The reviewer did not implement those proofs or their verifier.

`route_verdict: established`

`evidence_scope: independently reviewed smooth global steady Euler ring
with a nonzero literal constant-curl inner closed core and persistent
unknotted inner torus, in its actual translating steady frame`

## Strongest supported result

Accept the frozen0211 construction as stated. It establishes a geometric
object on the actual global field, not a local approximate cylinder:
compact smooth toroidal vorticity, a nonzero closed circular core, and an
open inner neighborhood satisfying curl u = sigma lambda u for one
constant lambda. A suitable inner invariant torus has nonzero physical
flux twist and Diophantine boundary rotation. Its boundary persists under
the stated sufficiently small analytic divergence-free perturbations;
the nondegenerate elliptic periodic core persists separately.

The outer variable curl factor and uniform far velocity remain part of
this same global solution. Neither a globally constant-curl compactly
supported vorticity field nor persistence of every irrational core
frequency is asserted. These distinctions are already in the submitted
proof, not additional restrictions imposed by this review.

## Source and equation audit

Read both proof bodies completely, the entire final verifier, both raw
24-check outputs, the exact first-source change, receipt, and the archived
EPS Theorem7.6 and equations7.28--7.32 directly. The archived primary text
digest agrees with its receipt. The unchanged0186/0195 global Green
construction supplies the stated source-to-potential estimate; the new
seed's invertibility is proved here and is not borrowed from the old
plateau profile.

1. **New seed and mass border.** The flat taper gives H=phi-delta/2
   exactly above delta. Consequently f=GG' has the stated affine inner
   equation and Q=f' includes the taper product derivative. The bound on
   Q is uniform in delta, whereas higher smooth norms are used only after
   fixing delta. At the taper entry the Bessel radial derivative is
   strictly negative. Integrating the radial ODE across the O(delta)
   layer establishes the exterior logarithmic coefficient and its strict
   nonzero limit. No sampled eigenvalue substitutes for this argument.

2. **Complete inverse.** The positive translation solution V=-phi'
   gives the all-angular ground-state identity with its boundary term.
   Regularity and exterior decay remove that term for the kernel test.
   Higher angular orders have the strict m^2-1 correction. The mass row
   removes the radial logarithmic solution, and the center pairing is
   exactly -Gamma. The translation cokernel pairing also fixes the linear
   harmonic drift in the bordered equation; the remaining constant is
   fixed by the radial Green equation. This is the full mass/center border,
   not an inverse on a hand-selected radial subspace.

3. **Global field and frame.** The effective source R^2 f/r^2 gives the
   actual force-free stream equation. The proved new inverse allows the
   existing global kernel continuation after delta is fixed. The outer
   maximum principle excludes spurious positive source components. Flat
   current at the source edge gives smooth velocity/vorticity and zero
   swirl on the axis. The far field is -U e_z in the steady frame, and
   E-U I_z is retained. Direct cylindrical curl differentiation verifies
   the literal inner curl identity for both signs without a boost.

4. **Physical return map.** The toroidal advance integrates the actual
   nonuniform transit and r metric. Coarea gives J_R'(h) rho_R(h)=-R.
   Thus the twist is differentiated against section flux, not area. The
   Bessel expansion yields the strictly positive reduced central twist
   lambda^2/(8A). Fixing a positive-radius annulus before R makes contour
   differentiation legitimate. The actual core return ratio has derivative
   lambda/2+O(log R/R), allowing a nonresonant radius choice without any
   arithmetic conjecture about Bessel roots.

5. **Persistence hypotheses.** In the inner affine region the elliptic
   equation and coefficients are analytic, so the velocity and return map
   are analytic there. The EPS theorem requires an analytic invariant
   Diophantine circle, nonzero normal torsion and the same preserved
   positive measure. The proof verifies these and explicitly applies the
   nearby-flux Moser identification before the same-measure theorem,
   exactly as in the primary source. A slightly larger tube leaves room
   for the perturbed invariant boundary. Nonunit core multipliers give
   the separate periodic-core IFT; strict ellipticity is open for its
   area-preserving transverse map. No Diophantine-core robustness is
   inferred.

The symbolic oracle checks the derivative, source, Bessel, full angular,
cylindrical curl and physical flux identities used by this analytic proof.
Its coarea check is an algebraic normalization anchor, not an independent
existence proof. The all-m strengthening correctly replaces three sample
orders with an exact nonnegative-integer polynomial. Reuse both recorded
successful executions; rerunning unchanged algebra would add no evidence.

## Frozen boundary and continuation

- affine-inner-seed.md SHA256:
  `a477727fa9a5d90becf2db8dea7b62882eea19ac6ca56f5ef421b50d3ed3dc52`
- same-ring-beltrami-torus.md SHA256:
  `8c2fa0c980c740e4d6acd55da324b66ac4a93c9a2862f0c9b2a816cf553291a6`
- verify.py SHA256:
  `79d055d7c86d9872ca0f693895ca5becc5f336a9caf283c9a69fca0b8fba6aea`
- strengthened.stdout SHA256:
  `f508240b991e9972416ee6abaa23db3a29ca2a6ebc76b212285cd7dea8731c2d`
- primary EPS text SHA256:
  `854f09f74d8806d5dd694c6f72c85eb8a63c3b37d443d8fc331b94dd41f549a1`

No scientific correction is requested. This positive geometry supplies
the same-field home and structural profile identities for0213's actual
Kelvin residual/current calculation. Its optical action and the coupled
continuum are distinct active parent obligations, not acceptance criteria
newly added to this geometric theorem. No unrelated replay or full
validation was needed for this review. Promotion will materialize the
reusable geometric definitions without changing the proven statement.
