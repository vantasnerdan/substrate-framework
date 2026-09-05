# 0053 — isotropic action and second-gradient kinetic normal form

Main owns this fixed algebraic theorem and its importable API. It continues
N3/N4 from the complete orbit reductions in 0048–0052, without assuming that
those microscopic constructions already supply every joint coefficient.
Parent objective and base release remain those of 0035 (v0.171.0).

Freeze before implementation: average a scalar cell angle about a unit axis
and its symmetric positive gradient matrix over simultaneous spatial
rotations. Derive the isotropic fourth-order tensor by invariant contractions,
then check it independently on exact spherical moments. Derive the action
normal form of a time-even, reflection-paired isotropic translation/spin
quadratic pencil through spatial order two. Keep all kinetic gradient and
translation/spin cross terms until the field transformation is displayed.
The transformation preserves uniform rigid translations and rotations but
does change nonuniform coarse observables; its physical field map remains
part of the result.

Inputs are the complete reduced Euler-action coefficients, not phenomenological
values. The theorem licenses exact coefficients of the declared slow-gradient
expansion, not exact equality at arbitrary wave number, a global assembly of
Euler cells, or deletion of boundary currents. The common positive spin
inertia is a prerequisite; the structure-free branch is taken before any
division by it. No empirical comparator, numerical remainder, tolerance,
fitted coefficient, or solver is involved. Small-ratio prescriptions are
satisfied by exact symbolic identities and explicit positive input bounds.

Route: derive invariant contractions and pull back both kinetic and potential
matrices by the same near-identity map. Independent symbolic variation,
wrong-sign/omitted-gradient mutations, and direct tensor moment checks are
the strongest practical oracle. Canonical implementation extends only the
new unpromoted euler_orbit.py; direct consumers are its new test file and
this attempt. Existing accepted modules/claims remain unchanged.

Result: established as stated. All nine euler_orbit tests pass in the saved
first-run receipt; five new cases cover this extension and four preserve the
0049 boundary. Independent review found no load-bearing correction. GitNexus
reported UNKNOWN because the new symbols are absent from its old index;
direct source search found tests and the new 0052 caller only. Ruff and diff
checks pass. This licenses the coefficient map, not the missing assembly.
