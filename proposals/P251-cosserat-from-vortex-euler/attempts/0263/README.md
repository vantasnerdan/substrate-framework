# 0263 — one-sided characteristic trace estimate and gluing

Owner: `herdr geometry-review pane w3:p4`, continuing the independent
0253/0255/0261 finite-radius geometry route.  Write scope is attempt 0263 only.
Base checkpoint is `754558d`; the parent is the unchanged issue200/P251
objective.  No prior attempt, source module, central proposal entry, registry,
release, generated file, or commit is changed here.

This contract is frozen before opening any new external source body or carrying
out new estimate/probabilistic analysis.  The coordinator must add the central
0263 entry and a successful schema receipt before substantive work begins.
Only the already opened 0261 equations and source licenses are inherited at
their recorded scope.

## Four frozen levels

The parent objective remains one stationary smooth constant-density Euler
ensemble supporting actual coupled histories, inherited physical action and
current, and compact Euclidean nonzero finite-core tubes at usable positive
measured density.  The periodic supplier, 0250 response integration, and other
root constructions remain independent.  This attempt cannot reduce the parent
to a local regularity or stationary geometry statement.

The active obligation is the precise missing lemma in 0261: a one-sided,
full-operator estimate that constructs the physical moderate solution,
determines its actual free-boundary trace from the interior data and source,
sums all angular modes with finite derivative loss, persists under the
nonradial profile/Hanzawa perturbations, and glues to the interior physical
spectral inverse.

This concrete attempt compares:

1. **V — exact Volterra/mode construction.**  Retain the full normal
   second derivative, prove a uniform Fourier trace estimate for the exact
   radial collar operator, and perturb it in an intrinsic scale.
2. **P — inward-diffusion exit construction.**  Treat the full drift/sum-of-
   squares operator as a killed diffusion with the outer edge an entrance
   boundary and the inner collar boundary the exit surface; derive the trace
   from the exit distribution and occupation kernel, then prove its finite-loss
   smooth bounds and perturbation stability.

The candidates may be combined if the probabilistic representation proves the
one-sided trace bound while the Volterra formula supplies sharp radial
normalization and physical-domain identification.  A route verdict applies
only to the route analyzed.  No failure implies exhaustion.

## Frozen base operator and physical domain

On the exact straight logarithmic collar let

\[
 d=a-\rho,\qquad 0<d<\delta_1<a,                    \tag{1}
\]

and retain the complete 0261 Fourier operator

\[
 \mathscr L_m v=d^2v''+
 \left(2-4d-\frac{d^2}{a-d}\right)v'
 +\frac{d^2(1-m^2)}{(a-d)^2}v.                       \tag{2}
\]

Its exact integrating factor and moderate equation are

\[
 \mu(d)=(a-d)d^{-4}e^{-2/d},                          \tag{3}
\]

\[
 (\mu v_m')'+\mu\frac{1-m^2}{(a-d)^2}v_m
       =\frac{\mu}{d^2}F_m,                          \tag{4}
\]

\[
 \mu(d)v_m'(d)=\int_0^d\mu(s)\left[
  \frac{F_m(s)}{s^2}
 +\frac{m^2-1}{(a-s)^2}v_m(s)\right]ds.             \tag{5}
\]

The discarded integration constant gives
`v_m' proportional to e^{2/d}d^4/(a-d)` and is excluded by the
physical quadratic-form domain under the exact conjugacy

\[
 D\mathcal F(T_0)v=-\frac{T_0^4}{\Phi}
 A_0\left(\frac{\Phi}{T_0^2}v\right).              \tag{6}
\]

The attempt must prove both directions of this domain identification; it may
not define “moderate” only by deleting the unwanted solution after the fact.
The natural finite-radius measure remains
`(1+epsilon x)^{-1} dx dz`.

In physical collar coordinates, a nearby pulled-back linearization has

\[
 P v=\sum_{j=1}^2X_j^2v+X_0v+c v,                    \tag{7}
\]

\[
 X_j=T\partial_j,qquad
 X_0=(2-5T)\nabla T\mathbin\cdot\nabla
 -\frac{\epsilon T^2}{1+\epsilon x}\partial_x.    \tag{8}
\]

On `T=0`, the eikonal gives `|grad T|^2>0` and
`[X_0,X_j]=2|grad T|^2 partial_j`.  These inherited brackets license
testing an intrinsic estimate; they do not license an arbitrary two-sided
extension or a boundary inverse.

## Candidate V — full-mode Volterra estimate

Candidate V must keep (2)'s `d^2v''` term.  It must:

1. construct, for each Fourier mode, the unique physical moderate solution
   from source plus Cauchy data at `d=delta_1`;
2. derive constants uniform in `m`, including the transition scale
   `d about |m|^{-2/3}`, rather than summing modewise existence statements;
3. prove the selected boundary trace is a finite-loss or smoothing operator
   of the declared angular Sobolev/Hölder order;
4. control the `m=0` block using the 0261 radial spectral gap, the `m=1`
   block with its translation projection/gauge, and `m>=2` without replacing
   the full equation by the first-order toy model;
5. perturb the estimate to the nonradial coefficients induced by the Hanzawa
   graph and nearby profile, with a quantitative smallness norm; and
6. show that the collar trace and conormal derivative match the uniformly
   elliptic interior solution with no extra boundary condition.

The exact radial formula (5) is a starting identity, not the requested
`m`-uniform theorem.

## Candidate P — inward-diffusion exit representation

Interpret the second-order and drift terms in (7) as the backward generator of
a diffusion in the one-sided collar.  At the straight edge the normal model
has diffusion coefficient proportional to `d` and inward drift tending to
`+2 partial_d`; tangential noise has size `d/(a-d)`.  The outer boundary
`d=0` is therefore a candidate entrance/characteristic boundary, while
`d=delta_1` is the declared exit surface carrying the interior matching data.

Candidate P must derive, rather than assume:

1. the exact stochastic differential operator, including Itô/Stratonovich
   drift corrections, the zeroth-order term, and the sign convention relating
   its generator to (7);
2. existence up to the inner exit time, nonattainment or the correct boundary
   classification at `d=0`, and finiteness of the required exit/occupation
   expectations;
3. a Feynman--Kac/Dynkin formula with any killing shift stated explicitly, so
   negative physical eigenvalues or a sign-indefinite zeroth-order coefficient
   are not hidden behind a maximum-principle slogan;
4. that the exit law and occupation kernel select the actual boundary trace
   from the source and `d=delta_1` data, rather than prescribing arbitrary
   Cauchy data at `d=0`;
5. finite-loss smooth estimates for this trace.  Hörmander smoothness of a
   transition density alone is insufficient without bounds uniform as the
   starting point approaches the characteristic edge and uniform in the
   Hanzawa/profile perturbation;
6. equivalence of the probabilistic solution with the physical moderate
   quadratic-form solution in (6), followed by collar/interior gluing.

The toy first-order smoothing kernel from 0261 is a limiting check, not a
permitted replacement for the diffusion containing the full normal second
derivative.

## Corrected usable estimate

The `k_0`-only interior-data term in 0261 equation (34) is not frozen as a
valid standalone estimate.  Local high-regularity control ordinarily requires
high-order matching data unless an interior elliptic estimate is composed with
it.  Freeze two collar widths

\[
 0<\delta_0<\delta_1,                                \tag{9}
\]

and first target the honest local estimate

\[
 \begin{aligned}
 &\|v\|_{\mathcal X^{s+2}(0,\delta_0)}
 +\|\gamma_0v\|_{\mathcal B^{s+\sigma}(\partial D)}
 \\
 &\qquad\le C_s\Big(
   \|Pv\|_{\mathcal Y^{s+r}(0,\delta_1)}
  +\|\gamma_{\delta_1}v\|_{H^{s+r_1}}
  +\|\gamma_{\delta_1}\partial_dv\|_{H^{s+r_2}}
  +\|v\|_{L^2(0,\delta_1)}\Big),                   \tag{10}
 \end{aligned}
\]

where `mathcal X,Y,B`, `sigma`, and the finite losses
`r,r_1,r_2` are derived from the full operator.  They are not selected after
seeing a desired conclusion.

The campaign-useful global estimate may replace the high-order inner traces by
a low norm only **after** uniformly elliptic interior regularity on a still
larger region has bounded them by the global source plus the physical
`L^2` norm.  The collar shrink `delta_0<delta_1` and that interior
composition must be explicit.  Thus an acceptable final form is

\[
 \|v\|_{\mathcal X^{s+2}(0,\delta_0)}
 +\|\gamma_0v\|_{\mathcal B^{s+\sigma}}
 \le C_s\left(\|Pv\|_{\mathcal Y^{s+r}(D)}
                 +\|v\|_{L^2(D)}\right),            \tag{11}
\]

only when accompanied by the interior estimate and range projection that make
(11) true.

## Perturbation, gluing, and nonlinear consumer

Use the 0261 Hanzawa chart with an unknown boundary graph.  The cosine graph
coefficient is the radial-center gauge, evenness removes the vertical
translation, and every other trace mode remains selected by the collar/interior
solve.  Candidate success requires:

1. an estimate (10), and then (11), uniform on a stated neighborhood of the
   balanced radially nondegenerate profile;
2. a global right inverse modulo the one even translation cokernel in the
   physical weighted domain;
3. smooth/tame dependence on the profile coordinate, boundary graph, and
   iterate, with a fixed finite derivative loss;
4. compatibility with 0261's desingularized scalar equation
   `S(epsilon,beta)=epsilon S_tilde(epsilon,beta)` and
   `S_tilde(0,beta) proportional to Q(Phi_beta)`.

Passing these rows licenses the range inverse consumed by the 0261 nonlinear
construction.  It does not by itself solve the scalar balance equation,
construct the finite-`R` Euler field, or transfer the 0258 twist; those
conclusions follow only after their already stated downstream hypotheses are
checked on the same solution.

## Frozen source and verification policy

After activation, inspect primary sources for:

- boundary classification and killed-diffusion exit representations for
  degenerate generators with inward drift;
- smooth exit densities or Poisson kernels under a drift Hörmander condition;
- one-sided maximum principles and intrinsic boundary estimates for
  nonnegative-characteristic operators; and
- perturbation/gluing of the resulting trace operator.

Every theorem is transferred with its domain, coefficient regularity,
boundary classification, killing/sign, and uniformity hypotheses.  A result
about an interior transition density is not a boundary Poisson estimate.

The oracle is analytic: exact generator/sign calculation, a proved
one-sided representation, an `m`-uniform estimate with collar shrink, domain
equivalence, and perturbative gluing.  No numerical remainder is frozen, so
the small-ratio numerical skill does not bind and no numerical production is
licensed.  A short exact symbolic check may be used after activation only for
a concrete drift or conjugacy uncertainty.

## Typed ladder and maximum verdict

| Step | Positive intent | Pass licenses | Does not license | Failure scope | Unlocks |
| --- | --- | --- | --- | --- | --- |
| 1 | One-sided physical trace construction | Moderate solution and actual edge trace from source/interior data | Uniform mode summation or perturbation | Exact radial Volterra or killed-diffusion representation | Estimate (10) |
| 2 | Finite-loss collar estimate | Full-mode (10) with explicit losses and collar shrink | Global inverse without interior/range gluing | Declared intrinsic scale | Global estimate (11) |
| 3 | Perturbed physical right inverse | Uniform Hanzawa/profile inverse modulo translation | Finite-`R` solution before scalar/nonlinear closure | Declared perturbation neighborhood | 0261 Lyapunov--Schmidt consumer |

`requires`: activated 0261 exact operator, physical conjugacy, radial spectral
gap, boundary graph, and scalar reduction; new primary sources only after
activation.

`pass_licenses`: at maximum, the actual smooth/tame one-sided right-inverse
family needed by the 0261 finite-radius construction.

`does_not_license`: arbitrary two-sided extension; prescribed zero trace;
deleting `d^2v_{dd}`; using a smooth transition density without boundary
bounds; finite-`R` Euler existence; 0258 twist transfer; compact density;
0250 response/action/current; parent completion; promotion; or exhaustion.

`maximum_verdict`: an analytic, perturbatively stable, one-sided
full-operator estimate and physical collar/interior right inverse satisfying
(10)--(11).

`failure_scope`: the exact radial Volterra summation route or the declared
killed-diffusion/maximum-principle representation and their perturbation
spaces.  Failure of either activates the other or a new boundary calculus; it
does not refute the finite-radius obligation.

`route_verdict: pending central/schema activation and substantive estimate
construction`.

`evidence_scope: preregistered full normal operator, Volterra and probabilistic
exit candidates, corrected collar-shrink estimate, physical-domain
equivalence, perturbation/gluing, and parent boundary only`.

Next coordination achievement: `w3:p1` adds the central 0263 manifest entry
and a schema-valid receipt.  Until then this pane opens no new source body and
performs no substantive Volterra, probability, maximum-principle, or gluing
analysis.

## Append-only completion status

Central activation subsequently passed with `schema.exit=0` and 270 accepted
claims.  The substantive result is recorded in `trace-estimate.md`,
`perturbation-gluing.md`, `source-transfer.md`, and `receipt.md`.
The activation instruction explicitly authorized continuing to downstream
finite-`R` closure if the registered inverse was earned; the final
Candidate-I verdict below consumes that authorized continuation and does not
retroactively alter the blinded inverse criteria above.

`candidate_V_radial_route_verdict: established as stated`.

`candidate_P_representation_verdict: established as stated for trace
selection, boundary classification, exit moments, signs, and nearby-coefficient
C0 stability`.

`nonradial_tame_perturbation_verdict: established as stated by the invariant
two-sided stopped-resolvent construction and interior intrinsic estimates`.

`finite_R_candidate_I_verdict: established as stated; the remaining parent
dependency is the independent same-field response/action/current and
coupled-history join, including normalization of that response on the common
density ensemble`.

The same-field stationary density normalization is equation (11) of
`perturbation-gluing.md`: disjoint compact copies give positive number and
core-volume densities without cross terms.  The remaining normalization is
the independent response/action/current normalization on that ensemble.
