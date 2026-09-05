# 0162 — complete physical configuration/rate action chart

Root owns a new additive helper in euler_phase.py, its direct tests and
this new attempt. Base v0.176.0 at dbf0c04; C-CST-008..011 and all frozen
claims/attempts remain unchanged. This supports0158's actual joint phase
construction without replacing physical observations by canonical names.

The positive object is an exact n-coordinate physical configuration/rate
chart of a given actual finite Euler phase family, including the complete
moving symplectic form and Hamiltonian, measured momentum, configuration
Poisson bracket and the obstruction to an ordinary second-order action.
The caller supplies actual observation rows and their derivatives; the
helper does not construct an Euler mode, enforce a missing moment, or
assume physical configurations commute merely because their initial rows do.

One exact algebraic route suffices: derive the physical rate row from the
actual generator, invert the complete observation chart, and use the
existing moving_phase_pullback for every connection. Compare direct
Euler-Lagrange variation in two coupled physical coordinates and expose
a noncommuting-coordinate counterexample. Scalar reduction is checked
against the existing physical_scalar_chart without changing that API.
No empirical comparator or numerical eigenvalue/sign design enters.

New reuse surface: physical_configuration_chart and its result object.
Existing public definitions are unchanged; source search and GitNexus
impact guide tests of the existing phase/scalar-chart module consumers
plus the new helper. Follow the user's explicit2026-09-05 blast-radius
validation instruction: targeted tests and relevant structural checks,
reusing0160's full2597-test receipt rather than running another full suite.

Passing this child supplies a reusable actual-observation/action identity,
not the physical moment construction, homogeneous continuum closure or
parent completion.0158 retains ownership of the actual Euler construction.

Status: preregistered; central validation precedes implementation.
