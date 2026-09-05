# 0136 — stationary force-free realization of a shaped finite vortex core

Owner /root. Baseline a93e78e, accepted v0.175.0 and008..010 unchanged.
Positive target: construct actual stationary Euler cores whose vorticity
gradient changes the intrinsic optical response, while preserving a
force-free structural connection to stationary vortex tubes. The
uniform-rotation sign identity does not cover these shaped cores.

Candidates: a generalized Beltrami column with a freely shaped swirl
and a derived large axial velocity; a large-major-radius toroidal
semilinear continuation; and a strong-Beltrami perturbative realization.
The comparison distinguishes exact stationary Euler, generalized versus
constant-factor Beltrami structure, physical perturbation dynamics, and
global stationary tube existence. No nonconstant curl factor is silently
inserted into the accepted constant-factor APIs or EPS theorem.

Primary inventory before opening new bodies: Gallay–Smets,
Spectral stability of inviscid columnar vortices, arXiv1805.05064;
Enciso–Poyato–Soler, Stability results, almost global generalized
Beltrami fields and applications to vortex structures in the Euler
equations, arXiv1605.06626 / CMP360(2018),197–269. Search abstracts are
only route-generation evidence. Nearest source is0120's exact thin
force-free torus and0135's finite-core optical response continuation.

Oracle: analytic Euler/force-free ODEs, physical pressure, exact second
variation and available existence hypotheses. No numerical eigenvalue,
stability threshold or comparator is selected. A positive stationary
profile does not itself license an optical eigenmode or global EPS
embedding; those remain explicit downstream constructions.

Failure-derived global representation repair: a general stationary
planar Euler field with globally defined Bernoulli B=p/rho+|v|²/2
and B bounded above has the exact generalized force-free lift
u=(v,sqrt(2(C-B))) for C>sup B. This uses the actual Bernoulli function,
not a globally single-valued vorticity function of streamfunction.
The planar Euler subsystem and its fixed-orbit action remain exact.
This construction is extracted into an additive, explicitly unpromoted
`euler_forcefree` module with direct tests; no existing canonical API or
accepted claim changes. GitNexus reports zero incoming edges for nearby
stationarize_planar, but rg finds its tests; that graph omission is not
treated as absence of consumers. The new module has no existing callers.
