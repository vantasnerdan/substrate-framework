# 0112 — independent centroid motion and actual EPS-core dynamics

Owner `/root/smooth_core_review`, this directory only. Parent issue #198 /
P251, original positive objective, and accepted release v0.171.0 remain
fixed. Uses the physics-erdos-loop analytic/source workflow. No empirical
comparison, discretized spectrum, or numerical selection is involved.

First exact target: extend 0109 by the spatial-translation/Galilean symmetry
of the actual Euler solution, proving independent physical centroid initial
position/velocity and optical orientation at zero macroscopic wave number.
The nondecaying mean sector and the decaying optical Kelvin leaf remain
distinct. This is a symmetry proof, not a new independent-review request.

Next constructive target: an actual optical/material-moment sector inside
the SAME invariant EPS tube, not a remote packet with a prescribed tied
coordinate. Registered routes: (A) use the EPS local Beltrami construction
to realize a controlled near-Lundquist segment in that tube; (B) propagate
a WKB packet along its actual elliptic periodic core with the exact Euler
amplitude equations and pressure/complement remainder; (C) use a directly
compatible axisymmetric toroidal seed before robust knotted continuation.
Selection criteria are actual source hypotheses, Kelvin compatibility,
physical material angle and moment maps, finite-time full Euler control,
and whether any extra scale hierarchy can actually be constructed.

Source inventory: the archived primary 1210.6271v2 (EPS tube existence,
thin-tube Beltrami estimates, elliptic periodic orbit, global approximation)
and the previously reviewed source imports. Source-supported streamline
ellipticity will not be equated with Euler perturbation spectral stability.
The already accepted or frozen parent conditional claims and reviewer
files are outside this attempt's mutation boundary.

Result: `galilean-chart.md` establishes the exact independent mean chart.
`eps-core-floquet.md` executes the same-core route using the primary
source's free small eigenvalue, an exact unit-determinant Euler amplitude
return, quantitative trace margin, and arbitrary finite-order compact
Kelvin/pressure recursion. It includes actual material moment rows and
a constructive positive covariance yielding nonzero spin. `verify.py`
passes all 21 exact checks on its first execution, saved in `first-run.txt`;
Ruff passes. This is a finite-time same-core Floquet/packet construction,
not a claim that the remaining continuum objective is complete.
