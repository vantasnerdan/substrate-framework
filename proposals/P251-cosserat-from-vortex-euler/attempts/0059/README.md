# 0059 — stationary random-Beltrami common-rotation density

Parent P251 / issue #198; owner `/root/construction_review`; this directory
only. Complementary to 0057's one-action material assembly, not a review of
its eventual conclusion. Source: 2006.15033, already archived by 0057; read
that source without editing or duplicating it.

Frozen positive task: on the actual stationary finite-variance Gaussian
Beltrami ensemble, derive finite origin-independent KKS and energy-Hessian
densities for a common physical rotation paired with stationary local compact
reaction fields. Include hard-core good-patch selection, exact projection
off the three Euclidean translation moments, full Leray projection, and
Euler Lie–Poisson/Jacobi identities. Determine rather than presume whether
the common rotation is physical or merely changes realization labels.

Candidate A: a co-rotated random field and its selected compact generator
family, with density defined by stationary averaging and centered local
moment cancellations. Candidate B: a marked local material-frame rotation
and conditional law, retaining its physical spatial registration rather
than quotienting by isotropy. Structural selection is an actual physical
angular observable, nonzero finite KKS density, well-defined positive action,
one material mass normalization, and explicit ensemble/coherence premises.
An isotropic distribution by itself is not evidence of a free physical rotor.

Oracle: source theorems on stationarity, finite variance, support and
ergodicity; exact angular/translation moment identities, stationary
integration by parts, finite-patch action bounds and symplectic reduction.
No numerical eigenvalue, empirical comparator, fitted coefficient or assumed
favorable sign. The source introduction was inspected to identify this
candidate; no comparator value is used. Contract frozen before its density
or action calculation.

Status: active. All relative and thermodynamic limits must identify what is
held fixed, transformed, or independently varied.

Append-only replacement candidate (requested by `/root`): the original
positive object needs only U and Phi, not an independently cyclic global
rotation B. Use a compact physical relative-angle pair q,s and set
Phi=beta+q, beta=curl U/2. Construct a KKS-isotropic fixed right inverse of
all eleven incompressible affine moments (three translations and eight
trace-free linear fields), so the same Euler action retains no uncomputed
affine/momentum pairing. Its positive relative-rate kinetic term is carried
through the already derived 0053 low-gradient normal form, with boundary
currents retained. This avoids inferring an independent physical rotor from
the rotation invariance of a probability law.

## Completed route receipt

`stationary-relative-angle.md` contains the exact derivation. The original
global-B density route is **blocked at physical material-rotor identification**:
its finite origin-independent KKS pairing does not separate the actual orbital
motion of material centers from an ensemble-label rotation. This is a
route-specific missing construction, not a parent-object no-go.

The replacement route is **established as stated**: a physical compact
relative-angle/cage sector with all eleven incompressible affine KKS moments
removed by an explicit isotropic right inverse, unchanged core jet, positive
full stationary Euler action density, and full reaction-space elimination.
Its coefficients are `J=D* P^-1 D` and `K=H_QQ-N* P^-1 N`; the operator
contains all stationary Leray interactions. Independent time-reversed fluid
reactions cancel the odd gyro after variation. Isotropic full-action averaging
gives `J=j I`, `K=4 alpha I`; `Phi=curl U/2+q` supplies the physical relative-rate
kinetic matrix rather than an independently postulated common-body mass.

Evidence scope: exact conditional stationary/marked-ensemble microscopic
action construction and slow-affine input, not completed parent claim
promotion. Source results are read from 0057's archived 2006.15033; full
material pressure/current and spatial-gradient assembly are the consuming
0057/0062 work. The source Gaussian law is not silently declared rotationally
invariant: rotation/parity averages and independent reaction fields are stated.
No fitted modulus, empirical comparator, numerical eigenvalue threshold,
isolated-cell inverse, or all-wavelength invariance requirement is used.

`verify.py` imports the canonical `hermitian_schur_jet` and `CheckLedger`.
Its first execution, preserved in `stdout.txt`, passed **22/22 checks**, exit
zero, in 1.821 seconds. The exact checks corroborate the analytic proof and
expose wrong response isotropy, scalar-potential sign, prematurely tied
time-reversal reactions, and isolated inverse substitutions. Ruff and scoped
diff whitespace checks pass. These algebraic checks do not stand in for the
source support theorem or the functional-analytic stationary bounds proved
in the text. No numerical small-ratio design is present.

Next parent achievement: consume the same full action and force-moment
regularity in the material/gradient assembly, carry its complete reaction
Schur jets through the normal form, and review the resulting joint claim.
