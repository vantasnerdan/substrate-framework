# 0263 analytic receipt

Owner: `herdr geometry-review pane w3:p4`.  Central activation records
`schema.exit=0` and 270 accepted claims.  All new derivation is confined to
0263; no numerical production, source/central edit, commit, or long script was
made.

## Strongest supported statement

For the exact radial logarithmic collar, the one-sided full-operator trace
problem is established.  The proof retains `d^2 v_dd` and constructs two
explicit supersolutions.  They give, uniformly for `|m|>=2`,

\[
 |v_m(0)|\le Ce^{-c|m|}|v_m(\delta_1)|
 +C\langle m\rangle^{-2/3}
       \sup_{0<d<\delta_1}|(\mathscr L_mv_m)(d)|.
\]

Thus the inner-data Poisson operator is exponentially smoothing and the
source-to-edge trace gains `2/3` of an angular derivative.  The estimate is a
full second-order result, not the first-order toy kernel from 0261.  The two
low blocks are finite: `m=1` is the translation sector and `m=0` is controlled
by the 0261 radial spectral gap.

The exact Itô and Stratonovich generators also check.  The normal process has
quadratic diffusion and strictly inward drift, so `T=0` is entrance-only and
`T=delta_1` is the exit.  Its Feynman--Kac formula selects the actual outer
trace from the source and inner Dirichlet datum.  A uniform exit-time and
exponential-moment bound gives a `C^0` representation for nearby smooth
Hanzawa/profile coefficients even when the zeroth-order term has no sign.
The explosive scalar branch maps to an `e^{1/d}` physical field and is outside
the form domain; the bounded branch maps to a flat physical field.  Hence the
probabilistic and physical selections agree at the radial base.

The correct usable high-regularity statement has a collar shrink
`delta_0<delta_1`.  The positive half-collar is invariant under the diffusion,
so a stopped resolvent on a smooth two-sided coefficient extension restricts
to the physical one-sided solution and gives an actual smooth extension
across `T=0`.  Interior hypoelliptic estimates on the shrunken collar therefore
need only a low norm of the remote inner datum.  This proves the corrected
`k_0` form rather than assuming an arbitrary extension.  In the global
physical solve, interior elliptic regularity controls that datum, and conormal
gluing is automatic because collar and interior are restrictions of one form
solution.

## Route verdicts

`candidate_V_radial_route_verdict: established as stated`.

`candidate_P_representation_verdict: established as stated for the exact
generator, entrance classification, selected trace, exit moments, killing
sign, and coefficient-uniform C0 bound`.

`nonradial_tame_perturbation_verdict: established as stated in the intrinsic
projective-limit scale`.

The constructed extension converts the boundary point into an interior
Hörmander point without changing the positive-side process.  Uniform local
estimates control the trace and all intrinsic derivatives, and differentiation
of the extended equation gives smooth tame parameter dependence.  The
coefficient perturbation is no-loss on this intrinsic graph scale, so the
base right inverse perturbs by a Neumann/resolvent argument while the finite
translation projection persists.

`finite_R_candidate_I_verdict: established as stated`.

The 0261 range equation, Hadamard division, simple 0256/0261 balance border,
and scalar implicit-function theorem now give an actual solution for all
sufficiently large finite `R`.  Its meridional support is one axis-separated
disk, its revolution is a compact solid torus, and its velocity, vorticity,
and full pressure have smooth zero extensions.  Smooth convergence on the
retained Bessel core supplies the 0258 `C^4` premise, so the same field has
nonzero normalized flux-action twist on a fixed regular core annulus.

## Parent and downstream boundary

The actual finite-`R` compact Euler carrier and retained twist are established
on the declared dependency closure.  Disjoint compact-support copies give the
stationary positive-density supplier with no cross field or pressure term.
For lattice spacing `L` exceeding twice the support radius and independent
occupancy probability `p_occ`, the exact number density and measured core
fraction are `p_occ/L^3` and `p_occ V_core/L^3`; both are positive and require
no amplitude or pressure renormalization.

The same-field 0250 response, physical action/current, actual coupled
histories, their common normalization, parent completion, promotion, and
exhaustion remain outside 0263.

`evidence_scope: exact radial full-normal trace smoothing; exact nearby
entrance-diffusion representation; constructed invariant two-sided resolvent;
uniform intrinsic trace/right-inverse family; physical gluing; finite-R
compact Euler carrier; retained 0258 twist; and exact disjoint-copy stationary
density normalization`.

`next missing construction: consume this stationary compact positive-density
supplier in the independent same-field 0250 response/action/current and
coupled-history join`.
