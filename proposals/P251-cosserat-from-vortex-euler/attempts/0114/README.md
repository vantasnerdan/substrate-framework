# 0114 — physical-frame packet action and its sign

Owner `/root/smooth_core_review`, this directory only. Parent P251 / issue
#198 and the stronger same-EPS-cell objective remain active. This is new
constructive frontier after 0112, not a reinterpretation of accepted claims.
Physics-erdos-loop governs the exact analytic construction. No numerical
eigenvalue, comparator, or sign selection from a Floquet log is used.

Target: the actual same-EPS packet KKS form and quadratic Euler action in a
fixed physical comoving frame, its energy/Krein sign, and its complete
physical angle/moment map. Preserve all carrier and frame-generator terms.
An autonomous canonical oscillator is not automatically a positive
physical-spin rotor; retain that distinction in the result.

Frame frozen before sign selection: the geometric Frenet normal/binormal
frame of the oriented EPS core (nonzero curvature is a source hypothesis),
followed by the zero-winding periodic uniform-twist registration obtained
by subtracting the time integral of its actual torsion from its own time
average. No additional integer frame winding is allowed in this route.
The carrier is the transported phase fixed in 0112, not an unrecorded
rotating factor. Floquet continuation is anchored to the unwrapped
geometric total torsion in this frame, not to a post-selected principal
logarithm.

Registered physical candidates: the two nonzero curl-eigenvalue signs
allowed by EPS Theorem 6.8 at fixed sufficiently thin geometry, and the
oppositely oriented/mirrored geometric torsion candidates. Selection uses
the derived fixed-frame KKS/Hamiltonian sign, while preserving the same
topological and finite-time Euler control obligations. The sign of a
coordinate frequency alone is not a selection oracle.

Main owns the generic exact moving-frame pullback in 0115; this attempt
consumes that algebra when available and derives the actual packet
coefficients, physical sign choice and observations. Earlier attempts and
accepted claim files remain immutable and outside this mutation boundary.

Result: `physical-frame-action.md` establishes a positive nondegenerate
same-EPS angular packet pullback for the physically selected curl-sign
candidate, with exact finite-interval autonomous normalization close to
the frozen frame and all physical moment/current rows retained. It also
constructs a particle-monodromy-invariant positive covariance that gives
nonconstant actual material spin, not merely nonzero squared spin.

`verify.py` passes 22/22 exact checks; Ruff passes. The first execution's
factored-versus-expanded structural-equality checker defect is preserved
in `first-run.txt`, and the corrected execution plus added nonconstant-spin
continuation is preserved in `corrected-run.txt`. There was no numerical
solver, empirical comparator, or numerical Krein-sign selection.
