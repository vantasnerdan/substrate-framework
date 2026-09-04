# Attempt 0036: finite-core and smooth continuation of the hexagon optical action

## Frozen delegated contract

Positive deliverable: construct the finite-core and then smooth Euler
continuation of attempt 0032's positive optical angle action, with exact source
hypotheses and spectral/variational persistence proved at the scope actually
used. Parent objective remains the complete exact smooth-Euler micropolar
continuum with translational action and stationary EPS compatibility.

Owner: delegated `orientation_construction`; authorized write surface: this new
directory only. Parent owns proposal/shared APIs and attempts 0035/0037. Base
commit: `3626fbf`; accepted release: v0.171.0. Sources: proposal P251, attempt
0032 exact six-vortex action, and primary theorem sources García
arXiv:2010.07361, Long–Wang–Zeng arXiv:1809.06425, and Radu–Stevenson
arXiv:2511.18592. These external papers are route inputs rather than accepted
framework claims.

The fixed microscopic conventions are constant density, incompressible Euler,
fixed circulation, isovortical variations, and positive point-vortex optical
Hessian at the regular hexagon. The optical mode is the alternating centroid
radial/tangential mode in a rotating frame. Existence of a steady patch,
smoothness of vorticity, and persistence of that spectral mode are distinct
obligations. No assumption of one will be used as proof of another.

The candidate route is dependency ordered: (1) audit the precise finite-patch
polygon theorem, (2) establish its constrained centroid/shape Hamiltonian and
spectral continuation, (3) audit and apply a smooth desingularization theorem,
(4) transport the optical action or identify its exact missing construction.
Alternative source routes are generic concentrated-vortex reduction and smooth
semilinear Euler desingularization. Structural fit, exact hypothesis coverage,
and preservation of the physical symplectic form control selection. No
empirical comparator, target frequency, or fitted coefficient is an input.

The oracle is exact theorem-hypothesis matching plus symbolic derivation of the
full variational/Poisson objects. Numerical work is not yet licensed: the
analytic representation and residual proposition will be recorded before any
solver or spectrum approximation is designed. Source excerpts, computations,
and route outcomes are appended below with their exact scope. The parent
campaign remains active throughout.

## Active result record

The source-only route was replaced after exact hypothesis matching exposed
its missing multiple-component and spectral-transfer steps. This is an
append-only candidate expansion motivated by those concrete mechanisms:
construct a smooth radial core with a full bordered inverse, solve the exact
nonlinear polygon equation directly, then pull back the physical Euler
Hamiltonian action to the declared affine molecular variations.

The resulting analytic artifacts are:

- `source-audit.md`: source hypotheses, scope and saved PDF digests.
- `radial-core-gap.md`: explicit separation of nontranslation radial-core
  modes, including the weighted and ordinary centroid gauges.
- `smooth-polygon-construction.md`: exact C-infinity, compact-vorticity,
  rotating polygon construction using a smooth Lane--Emden profile and two
  bordered implicit-function steps, with circulation fixed.
- `finite-core-angle-action.md`: exact isovortical affine trial-family Euler
  action, positive finite-core optical inertia and stiffness derived from
  interaction integrals, and the distinction from full-PDE spectral closure.

The strongest current result is an exact smooth finite-core polygon plus its
Euler-derived affine angle action. This closes the finite-core existence and
restricted microscopic action routes analytically, subject to independent
review of the stated proof. It does not assert the full continuum objective,
an invariant optical subspace of unrestricted Euler, or stationary EPS
compatibility. No finite-wavelength equality outside the original slow-varying
affine coarse-graining objective is added as a new requirement.

## Verification receipt

First executions of `radial_core_gap.py` and `smooth_core_construction.py`
passed 18/18 and 13/13 exact symbolic checks respectively. Outputs are
preserved in `gap-stdout.txt`, `gap-stderr.txt`, `construction-stdout.txt`,
and `construction-stderr.txt`. Ruff passes both scripts. These receipts check
the analytic identities and normalizations; the functional-analysis proof is
in the accompanying documents, not replaced by the tally. No numerical
spectrum, soft eigenvalue estimate, quadrature, or floating-point threshold
was used. The positivity error is a stated analytic cross-kernel bound.

Replay from the worktree root:

```
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/attempts/0036/radial_core_gap.py
PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P251-cosserat-from-vortex-euler/attempts/0036/smooth_core_construction.py
ruff check proposals/P251-cosserat-from-vortex-euler/attempts/0036
```

Independent review is requested for the exact bordered inverse and the
finite-core orbit-form restriction before any claim promotion. The next
parent construction uses this microscopic action with the explicitly stated
affine ensemble and stationary three-dimensional compatibility; this child
does not declare that parent construction complete.
