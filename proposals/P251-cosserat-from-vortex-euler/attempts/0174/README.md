# 0174 — same-EPS optical common-vector response

Base v0.177.0/dfe495c; parent objective and conditional continuum scope
unchanged. C-CST-011 establishes a finite-packet axial carrier response.
0172 requires its actual full common macroscopic K action before averaging
whole-field orientations, not independently rotated copies of K.

Target: construct the physical angle/current phase and full pressure/KKS
response through second derivatives in one common K, on that same finite
EPS cell, with transverse derivatives either computed or quantitatively
bounded below the positive axial curvature margin after averaging. The
choice of packet preparation and microscopic phase correctors is explicit;
neither phase-gauge invariance nor an axial projection is assumed to provide
the needed response. Averaging solution-column actions precedes elimination.

Candidates: locally transported transverse phase compensation with full
pressure tails; direct tilted-carrier Euler/Lin expansion; high-carrier
weighted-pressure/multipole localization of the cell Bloch response. Select
by exact physical observable, common-K consistency, retained microscopic
equations, error scale versus positive curvature, and natural EPS fit.
No empirical comparator. Stay analytic while these derivatives are tractable;
load small-ratio instructions before any numerical spectral design.

The ordinary local reduction gives fourth orientation moments1/5 and3/5
only after this common-K preparation and response has been licensed. The
construction may retain gradient mass and time-dependent current; those
terms are not an autonomous Cosserat law by naming them moduli.

Primary sources are frozen C-CST-011/0147 and finite-cell0153, the exact
constant-curl affine response0164 and actual0166 packet calculus where
its hypotheses transfer. Worker owns0174 only; root owns central metadata.
One native run per frozen verifier with captured output; blast-radius-only
validation. A route result activates the next parent obligation in-run.
