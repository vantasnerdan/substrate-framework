# Source and authority audit

The affine weighted theorem uses the standard smooth whole-space Euler local
existence theorem only for the unweighted `H^(s+1)` solution. The new result is
the a-priori dyadic weighted vorticity estimate and the zero-mean
Biot--Savart reconstruction displayed in `derivation.md`; it is not imported
from a named weighted theorem or inferred from unweighted continuity.

An earlier velocity formulation tried to propagate `s` weighted velocity
derivatives directly from the pressure stress. It lost one derivative because
`T=u tensor u` with `s` weighted derivatives gives `grad p` only `s-1`
weighted derivatives at the displayed far-field order. That route is
withdrawn. The final proof transports `q=curl(u-U_*)` through order `s-1` and
uses the order-minus-one Biot--Savart reconstruction, so no pressure derivative
enters the top-order estimate.

The exact pressure convention is

    p=(-Delta)^(-1) partial_i partial_j(u_i u_j),

with kinematic pressure and `G=1/(4*pi*|x|)`. The two derivatives eliminate
pressure monopole and dipole orders. The leading stress multipole and its
gradient are rederived by the attempt verifier and reusable API. They are an
independent low-order consistency check, not the analytic engine of the
weighted propagation theorem.

P253/0072 is used only at author scope for the smooth source-free completion,
harmonic transform, fixed-frame field, and exact cross-energy atoms. Its later
independent review controls these inputs before any parent promotion.

Roman Shvydkoy, *Homogeneous solutions to the 3D Euler system*,
arXiv:1510.03378v1, supplies the exact stationary sphere system at
`alpha=2` and only the axisymmetric classification. It does not establish a
general nonaxisymmetric no-go and is not used to prove weighted local
propagation. Cached source and version remain those pinned in `0072`.

No empirical comparator, production numerical result, KKS theorem, particle
force, or quantum premise is imported.
