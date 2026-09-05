# 0062 — complete second-gradient momentum Schur jets

Main owns this exact canonical-algebra extraction for 0057/0059. Freeze a
Hermitian quadratic action with C² coefficient jets P(k), N(k), H(k), where
P(0)>0 is the full reaction-momentum block and the Taylor convention is
P=P0+k P1+k² P2+o(k²). Derive the inverse and reduced Schur jets before
isotropic averaging or kinetic normalization. In particular, retain both
the P1 P0^-1 P1 inverse term and all mixed N1 terms.

The microscopic stationary C² proof and coercivity are supplied separately
by the full-fluid construction; this API does not manufacture those inputs
or replace a nonlocal operator by independent cells. The matrix formula
also states the bounded-operator identity when that construction licenses
it. Exact symbolic matrix inversion provides an independent oracle, with
a noncommuting complex-Hermitian example and omitted-term mutation.
No numerical remainder, tolerance or comparator is used.

Impact: additive helper in euler_orbit.py and its tests. The freshly indexed
pr199-completion graph locates the module; direct source search still
supplies more complete caller coverage than the graph's incomplete test
edges. Existing reduction/normal-form functions stay unchanged.
