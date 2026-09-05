# C-CST-011 — individual scientific review

Reviewer: `/root/construction_review`  
Date: 2026-09-05  
Transaction: the new 0147 construction, as proposed in 0157/README.md and
corrected by 0157/claim-correction.md.  
Verdict: **established as corrected; accepted for scientific promotion**.
Extraction, tests, registry promotion and release remain the parent's work.

## Independence and frozen boundary

This is one substantive pass and one correction check. I did not author or
implement 0147. I authored the unchanged 0145 twist dependency; that dependency
is an explicitly declared input, not a claim independently reviewed here.
This review assesses its new use in 0147. Accepted C-CST-008/009/010 and the
base release v0.175.0 are unchanged. No new empirical comparator or numerical
selection is used.

Reviewed raw artifacts and SHA-256:

| Artifact under attempts/0147 | SHA-256 |
| --- | --- |
| finite-packet.md | `5ea3aa86247c099e928b360274a0be6a78b76160aa4df34b35f1864f517ecd19` |
| toroidal-transfer.md | `be8f6d74a1e9d314c32e03e5e7cdde27340bbd8a8ea48585126e93f33f5a96b8` |
| verify.py | `a0267f2defa986819dadf682fbcfd89a7e58cbf06e92416361bb60bad6eeb498` |

The proof, verifier and captured repaired execution were read directly. The
recorded 24 exact checks pass. This review does not convert those algebraic
checks into a numerical PDE convergence claim; the analytic construction and
estimates are the load-bearing oracle.

## Strongest supported statement

The corrected proposed statement is supported: for a prescribed finite
optical window, sufficiently high but finite carrier produces an actual
finite-action, Kelvin-prepared Euler packet on a smooth stationary
constant-curl background with a robust unknotted material vortex torus.
The packet, its full velocity/pressure tails, the transported tag and the
optical observation concern the same tube. The registered physical
sheet-quadrupole angle has the claimed positive carrier curvature, and one
fixed nonnegative material tag supplies literal centroid angular momentum
matching the specified action row through two carrier derivatives, with
the stated finite-window error. The moving two-phase action retains its
time and parameter connections and has positive mass.

For the stated n=8 packet, with a=1/(1+c^2), the leading result is

    p_*^2 partial_p^2 gamma^2
      = a^2 sqrt(2) Omega^2 delta_*/3
        + O_T(delta_*^2 polylog),

with all geometric/approximation choices ordered after the finite packet
and its positive margins. The periodic background and its registered tube
family admit the stated stationary isotropic finite-energy-density law.
The localized packet is conditioned on a selected tube, as the correction
now explicitly says.

This is a meaningful microscopic positive result, not an autonomous
Cosserat closure, an exact monochromatic eigenmode, or a claim that the
literal tag spin is separately conserved. Those distinctions do not weaken
any conclusion actually constructed in 0147.

## Load-bearing checks

1. **Finite action and physical clock.** Plancherel uses the actual axial
   packet: beta_packet is the integral of squared spectral amplitude times
   fiber beta, rather than a fiber action multiplied by an assigned length.
   The Gaussian tag/packet convolution gives the factor a in the carrier
   derivative and hence a squared in the curvature. Real envelope
   derivatives remain in the mass and connection. The sign of the physical
   clock is not substituted for the different intrinsic-frequency
   curvature.

2. **Fixed physical tag.** Both material deformation contributions enter
   literal spin. A radial-only tag has the exhibited second-jet envelope
   mismatch; the two even axial controls repair it. The radial polynomial
   and exponential rows, the Laguerre coprimality argument and the two
   twenty-row blocks supply the needed rank. Reference-phase constraints
   and the J=7 remainder survive the two scaled carrier derivatives. The
   finite-dimensional implicit-function step is applied to the actual
   normalized observation, not to an independently chosen frequency.

3. **Same torus and full pressure.** The common-circle CK construction
   matches the Lundquist normal jets by positive angular quadrature and a
   finite Vandermonde correction. It retains a nonzero twist margin on the
   actual circle. The weighted high-angular-frequency pressure estimate
   controls the full Green operator and its far tails; it does not insert a
   radial wall or replace pressure by a local projector. In particular it
   addresses the otherwise fatal combination of growing global background
   norm and a merely algebraic local comparison error. Packet periodization
   uses integer angular modes with continuously varying spectral weights,
   not noninteger periodic modes.

4. **Actual histories and ordered transfer.** The arbitrary finite-order
   Kelvin preparation imported from 0142 allows the order to cover the
   finite derivative/conditioning losses before the small parameter is
   chosen. The full finite-time Euler comparison, material transport,
   centroid/frame corrections and nonlinear-amplitude-last ordering are
   retained. The later periodic approximation has the uniform global
   coefficient bounds needed by its local-plus-tail comparison. The new
   construction does not infer a uniform acoustic-time or all-wave-number
   theorem from these finite optical windows.

No further load-bearing equation, control-rank or pressure-transfer defect
was found in this pass.

## Two findings and the completed correction check

The original proposed wording admitted two stronger readings not supplied
by the proof:

- The eta=1/2 option could be read as a constructed counterpropagating
  standing pair. The evidence constructs a single finite spectral band;
  eta=1/2 is a permitted predeclared matching normalization. Construction
  of a second branch is a separate task.
- Phase translation modulo a periodic cell could be read as making the
  single localized packet stationary on all of R^3. It establishes the
  background/registered-family law, not that localized-packet law.

Both minimum corrections are now explicit in claim-correction.md, SHA-256
`a95a580cd83102aa5034b6f35357f8e3dec90d81fd99b0cbe4e14c63bcbed11d`.
I read and hash-checked that file. It resolves both findings without
changing equations, evidence, the quantitative microscopic claim, or the
user's objective. This is the single correction check; no scientific rerun
is indicated by this wording-only delta.

## Disposition and remaining integration

Verification: exact analytic construction with symbolic verification.  
Review: accepted, with the two incorporated corrections.  
Compatibility: additive microscopic extension; unchanged accepted inputs.  
Epistemic state: proposed until the parent completes canonical promotion.

Reusable extraction and its independent integral/derivative tests are
assigned to 0159. The parent still needs that implementation receipt,
impact-bounded consumer replay, and synchronized registry, release,
generated documentation and memory before claiming promotion complete.
This review supplies the individual scientific acceptance, not those
operational receipts.

The complete coupled-continuum objective remains active. Its acoustic and
joint-response constructions are not silently inherited from this optical
claim, and they are not hidden debt within the corrected C-CST-011 scope.

Signed: `/root/construction_review`.
