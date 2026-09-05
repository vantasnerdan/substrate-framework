# 0037 — exact smooth Euler displacement and stationary construction

Contract: 0035; parent N2–N7 stays active. Fixed exact subtheorem, with no
empirical comparator or production numerical remainder. Inputs are smooth
constant-density Euler, Euclidean metric, and stated material/boundary data.

## Full variational derivation before reduction

Let g(t,a) be a volume-preserving material map, u0 a smooth stationary Euler
solution, and xi(t,x)=delta g(t,g0^-1(x)) its displacement. The Eulerian velocity
perturbation is v=Dt xi-(xi.grad)u0, Dt=partial_t+u0.grad. Differentiating Euler
without dropping any block gives

    Dt v + (v.grad)u0 = Dt² xi - (xi.grad)((u0.grad)u0).

Thus its equation is rho Dt² xi+Hess(p0)xi+grad(delta p)=0, div xi=0.
This holds about every smooth EPS field on its actual domain; it does not
require the Rankine discontinuity or point-vortex approximation.

The same result follows from the second variation of the full constrained
material action integral[rho |g_t|²/2+p0(det Dg-1)]. For the local deformation
I+epsilon Dxi, det''=(div xi)²-tr((Dxi)²). On div xi=0,

    p0 tr((Dxi)²) = xi.Hess(p0).xi
                    + div[p0 (xi.grad)xi - xi (xi.grad p0)].

Periodic or compatible zero-boundary-flux variations therefore yield

    S2 = 1/2 integral[rho |Dt xi|² - xi.Hess(p0).xi] dx dt.

The positive rho |xi_t|²/2 is material inertia, not circulation-centroid
inertia. The mixed rho xi_t.(u0.grad xi) is retained; dropping it changes the
equations. The static quadratic form is Hess(p0) minus the background transport
norm, so positive time inertia alone does not establish a stable spin mode.
Pressure boundary flux remains for a moving subparcel, providing reaction.

## Planar stationary continuation

For a rotating Euler profile v=J grad psi, angular rate Omega, define the
changed physical flow w=v-Omega Jx. Direct momentum substitution gives

    p_w = p_v - 2 rho Omega psi + rho Omega² |x|²/2.

This is a genuine stationary Euler solution, with vorticity shifted by
-2 Omega. The far field/background changes; it is not the original flow
silently relabeled stationary. The construction extends a smooth rotating
planar core to a stationary planar core but does not yield knotted 3-D EPS
tubes. In 3-D, curl(e_z cross w)=-partial_z w for div w=0: the Coriolis term
generically cannot be absorbed in pressure. The planar extrusion is not a
nonzero constant-lambda Beltrami field.

## Transfer and verification

C-COL-001 and collective_coordinates.py establish conditional pullback
discipline, not an Euler profile or closure; only that discipline transfers.
The full action above now supplies the correct field-space input for a future
collective pullback. Existing vortex_dynamics.py is a singular 2-D model and
does not supply smooth EPS data. GitNexus reports new API names UNKNOWN; direct
search finds no previous definitions/callers. Impact is additive new module,
its new tests and this attempt, no accepted claims or releases.

Oracle: general symbolic transport identity and determinant second variation,
independently checked on a nonconstant-pressure smooth Euler field; wrong
pressure-Hessian and omitted gyroscopic terms must be exposed. Exact symbolic
calculus has no mesh, eigenpair, zero-mode floor, or precision uncertainty.
No continuum spectral sign is inferred from these checks.

Status: construction in progress. Next achievement: use this exact action to
derive the reaction-bearing collective projection and its closure conditions.

Completed child receipt: verify.py exits zero with 11/11 exact checks;
tests/test_euler_displacement.py adds five passing API tests. First-run logs
are retained here. Route verdict: established for the full displacement
identity, Jacobi action and stated planar stationarization. These are inputs
to the continuing joint-action construction, not N4 completion.
