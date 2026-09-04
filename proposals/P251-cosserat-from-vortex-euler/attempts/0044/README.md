# Attempt 0044: slow Bloch angle/shape action of the stationary Beltrami tube

## Frozen delegated contract

Parent objective remains P251's complete conditional smooth-Euler continuum,
including stationary EPS compatibility and the original slowly varying affine
ensemble. This child owns only new attempt 0044. Attempt 0040 and its
independent review are frozen inputs; no shared module is edited here.

Positive deliverable: derive the full Hermitian energy Hessian and
anti-Hermitian KKS form of the actual stationary tube's angle/shape sector
under exactly divergence-free slow Bloch variations. Compute optical
dispersion through second order in macroscopic wave number along each of
three axes, retaining the wave-number dependence of both forms. Identify
positive curvature coefficients or the named candidate failure, then repair
the actual representation if needed. No all-wave-number invariant Euler
closure requirement is added to the original affine objective.

Inputs: the explicit stationary Beltrami tube and physical angle/axial-return
generators of 0040; the fixed noncircular ratio `a=2b`; Fourier Leray
projection and the exact Euler orbit Hessian and KKS form. Slow generators
are defined by `xi_i(k)=curl_k A_i`, where the periodic Coulomb potential
has `A_i=curl xi_i/|wave|^2` mode by mode. This retains the return terms
required by incompressibility. No observed coefficient or target dispersion
is used to choose geometry or shape.

Exact finite Fourier algebra is the oracle. First derive full rational
matrices, verify Hermiticity and the zero-wave-number limit, then expand the
Hamiltonian characteristic equation. Chiral first derivatives are retained;
any later parity/isotropy average is a separate stated ensemble operation.
No numerical soft-mode calculation, optimizer, or floating-point threshold
is designed. The method-repair alternative is direct real-space derivative
and independent energy/KKS rederivation, not numerical threshold tuning.

## Active record

The complete Bloch Hessian/KKS calculation is preserved in `bloch_sector.py`
and `repaired-stdout.txt`. All three axis matrices are exact rational
functions. Their first derivatives vanish by computation; no chiral term
was removed by averaging. Both forms recover 0040 exactly at zero wave
number. The macroscopic return and shifted Leray projection are retained.

The oriented cell has negative axial optical curvature after the conjugate
shape and wave-number-dependent KKS form are eliminated. This candidate
result activated the original isotropic affine ensemble continuation.
`effective_action.py` derives the action coefficients before averaging and
shows that its isotropic periodic-domain curl/divergence action is strictly
positive. The explicit boundary divergence distinguishes this statement from
pointwise positivity on arbitrary affine gradient matrices. No geometry
parameter was retuned.

`bloch-action.md` contains the raw matrices, derivative inertia, collective
field normalization, negative oriented coefficient, exact isotropic tensor
average, positive curl/divergence representation, and physical units. It
also names the local core-axis angle observable precisely. The raw matrices
were supplied to the parent's common-angle/body construction; normalized
speeds are not substituted unchanged into that different coupled sector.

## Verification and route receipt

The first execution in `stdout.txt`/`stderr.txt` stopped on a SymPy structural
equality issue: `25/96` was represented with complex-unit factors. The
independently simplified difference was exactly zero. Replacing structural
equality by simplification of the matrix difference produced 27/27 checks
in `repaired-stdout.txt`; `repaired-stderr.txt` is empty. No equation or
numerical tolerance changed. The subsequent action continuation passed
15/15 first-run checks in `action-stdout.txt`; `action-stderr.txt` is empty.
An independent real-space Bloch-return and core-jet calculation passes 8/8
checks in `jet-stdout.txt`, with empty `jet-stderr.txt`; it records the
second-order difference between raw Bloch amplitude and literal relative
section angle, and proves the normalized coefficient invariant under that
field map. Ruff passes all three scripts. No numerical spectrum or energy
fitting is used.

Route verdicts: the all-parallel-cell positive axial-curvature route is
refuted at its explicitly negative coefficient; the original isotropic,
periodic affine action route is established at its stated internal-sector
scope, subject to independent review. This does not automatically complete
the parent common-angle, stationary EPS, or promotion obligations.

Reproduce from the worktree with `PYTHONPATH=src` and the repository venv
Python, running `bloch_sector.py` and `effective_action.py`. Both import the
frozen Fourier operations from 0040 rather than altering or copying them.
Also run `verify_core_jet.py` for the independent real-space observable check.
The sphere moment contraction reuses the importable homogenization module.
