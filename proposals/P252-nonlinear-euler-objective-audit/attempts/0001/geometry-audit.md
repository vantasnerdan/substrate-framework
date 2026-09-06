# Independent geometry audit: compact finite-radius supplier and nonlinear route

Reviewer boundary: fresh non-author audit of the geometry conjunct of
`C-CST-018`, reconstructed from attempts 0253, 0255, 0256, 0258, 0261, 0263,
and 0265 before reading the favorable 0268 review.  The accepted registry is
the pinned authority under audit, not evidence for the theorem.  This audit is
analytic; it consumes no numerical force, soft-mode, splitting, or stability
verdict, so the small-ratio numerical prescriptions do not bind.

## Verdict

`route_verdict: blocked with the missing construction named`.

`evidence_scope: the straight compacton, constrained radial nondegeneracy,
exact radial full-operator trace smoothing, conditional finite-ring twist, and
conditional disjoint periodic assembly survive; the nonlinear nonradial
free-boundary inverse and hence the claimed finite-R Bessel-core compact ring
are not established by the recorded proof`.

This is a proof gap, not a refutation of compactly supported steady Euler rings
and not a counterexample to the proposed finite-R branch.  It directly reaches
the first existential conjunct of `C-CST-018`: the registry says that there
exists an exact periodic assembly of finite-radius rings
(`governance/claims.yaml:15260-15268`) and identifies 0263 as its finite-R
existence proof (`governance/claims.yaml:15389-15393`).  The cited source chain
does not finish that proof.  If the missing right-inverse family is taken as a
declared hypothesis, the strongest valid result is the useful conditional
theorem stated below; the current evidence does not discharge that hypothesis.

## What the source chain does establish

1. **Exact straight compacton and the correct finite-radius remainder.**
   The inverse-designed radial profile gives a smooth, compact-velocity
   straight solution with an unchanged nonzero constant-curl Bessel core and a
   logarithmically flat edge (`0253/construction.md:22-127`).  Direct
   substitution shows that bending the circular levels produces the nonzero
   curvature residual (30), so a deformed free boundary is genuinely needed
   (`0253/construction.md:133-175`).  The defining-function equation
   (33)--(35) is an exact and useful regularization of the non-Lipschitz label
   law, while the same source explicitly says that it does not prove
   solvability (`0253/construction.md:177-231`).

2. **Translation compatibility and constrained radial nondegeneracy.**
   The first translation solvability row is the exact balance
   \(Q(\Phi)=0\), independently supported by the finite-radius virial identity
   (`0256/translation-balance.md:13-48,95-133`).  A fixed-core/fixed-edge
   profile family attains both signs and has a simple balanced member
   (`0256/translation-balance.md:146-171`).  In 0261, the physical radial
   endpoint is limit point with compact resolvent, and the three-weight Bessel
   buffer calculation makes the constraint/eigenvalue determinant nonzero
   (`0261/radial-construction.md:52-73,121-172`).  The finite-dimensional IFT
   then removes a possible radial zero while preserving the whole core, edge,
   support, and balance, and retains a separate transverse balance coordinate
   (`0261/radial-construction.md:174-210`).  Together with the ground-state
   factorization for \(m=1\) and strict positivity for \(m\ge2\)
   (`0255/first-order-construction.md:163-184`), this supports the claimed
   straight-operator kernel statement.  I found no contradiction in this
   radial part.

3. **A real all-mode radial edge estimate.**
   The exact Fourier operator retains \(d^2v''\).  The occupation barrier gives
   an \(O(|m|^{-2/3})\) source-to-edge bound and the inner-data barrier gives
   exponential smoothing (`0263/trace-estimate.md:93-198`).  Consequently the
   displayed estimate (23) for each high mode, with the finite low blocks
   separated, is supported (`0263/trace-estimate.md:200-229`).  The exact
   integrating factor also identifies the explosive branch and excludes it
   from the physical form domain (`0263/trace-estimate.md:236-252`).  These are
   substantive positive results and should be preserved in any correction.

4. **Conditional twist and assembly.**
   Given a finite-R solution with the specified local \(C^4\) convergence, the
   contour/area formulas transfer the nonzero Bessel-core flux-action twist
   without an additional tuning parameter (`0258/core-twist.md:44-92`).  Given
   a compact smooth template whose pressure is flat and exterior-constant,
   periodized disjoint copies solve stationary Euler exactly because all cross
   terms vanish pointwise (`0265/assembly-and-chart.md:8-29`).  The action-angle
   chart, finite-R frequency overlap, and cohomological denominator estimates
   are likewise valid conditional consequences of that template
   (`0265/assembly-and-chart.md:31-85`).  None of these implications supplies
   the missing template; 0265 declares that dependency at lines 1--6.

The strongest geometry theorem presently supported is therefore:

> There is an admissible balanced straight compacton with a fixed nonzero
> constant-curl Bessel core, logarithmically flat physical edge, no radial zero
> mode, only the two translation zero modes in the full straight
> linearization, and an exact high-mode one-sided trace estimate.  If the exact
> bordered nonradial Hanzawa linearization admits the uniform smooth tame right
> inverse claimed in 0263, then Lyapunov--Schmidt/Nash--Moser continuation gives
> the finite-radius compact ring; that ring conditionally retains nonzero
> flux-action twist and admits exact disjoint periodic assembly at positive
> density.

## The load-bearing gap

The gap begins precisely where 0261 left it.  That attempt says the missing
lemma is a one-sided full-operator Poisson/trace estimate uniform for nearby
coefficients, with declared intrinsic spaces, finite loss, physical-domain
equivalence, and collar/interior gluing (`0261/inverse-analysis.md:329-363`).
It correctly warns that without this lemma, Nash--Moser and compact Euclidean
existence are unsupported (`0261/inverse-analysis.md:365-382`).

Attempt 0263 proves enough to reduce that problem, but does not provide the
promised theorem:

* Its frozen contract requires definitions of \(\mathcal X,\mathcal Y,
  \mathcal B\), the gain/loss indices, and estimates uniform in the Hanzawa and
  profile perturbations (`0263/README.md:166-210`).  It separately requires a
  global right inverse, smooth/tame parameter dependence, and compatibility
  with the desingularized scalar equation (`0263/README.md:212-227`).
* The result instead describes the spaces only schematically, by the intrinsic
  derivatives they are intended to control, and states the interior estimate
  (26)--(27).  The only explanation of ordinary derivative loss is
  that the bracket conversion “costs finitely many derivatives at every
  order”; neither that loss nor the boundary trace scale is computed
  (`0263/trace-estimate.md:274-320`).  Thus there is no checkable domain and
  codomain on which the later inverse acts.
* The cited sources do not fill this hole.  The campaign's own transfer receipt
  says that they supply interior hypoellipticity and local derivative estimates,
  but not a killed characteristic-boundary Poisson estimate or Hanzawa-parameter
  tame bounds (`0263/source-transfer.md:39-54`).  It calls the invariant
  extension and uniform tame composition “new analytic work”
  (`0263/source-transfer.md:56-63`).  A stopped extension can plausibly turn
  each fixed one-sided solution into an interior distributional solution; local
  hypoellipticity alone does not prove the quantitative, coefficient-uniform
  boundary and tame estimates required by the frozen contract.
* The full perturbation is then compressed into
  \(R_q=(I+R_0K_q)^{-1}R_0\) and one tame product inequality
  (`0263/perturbation-gluing.md:60-82`).  Here \(R_0\) is only a projected
  inverse modulo translation.  No fixed bordered domain/codomain is defined;
  the varying kernel and cokernel projections and the graph trace variables
  are asserted to persist but are not inserted into the displayed operator.
  Small spectral clusters can persist under perturbation; an exact kernel need
  not persist at a nonsolution iterate or after the radial translation symmetry
  is broken at finite \(\epsilon\).  What is needed is the actual bordered
  inverse on the moving spectral subspace, with its projection derivatives and
  estimates, rather than the wording “kernel and cokernel persist.”
  The estimate needed to invert in a base norm and bootstrap through every
  higher intrinsic seminorm is not proved.
* Finally, the proof states that the range solution \(U(\epsilon,\beta)\)
  exists and immediately applies Hadamard division and a scalar IFT
  (`0263/perturbation-gluing.md:90-124`).  It never states a Nash--Moser theorem
  with the actual spaces, smoothing operators, neighborhood, nonlinear map,
  full derivative (including boundary graph and border variables), or a smooth
  tame right inverse satisfying that theorem's hypotheses.  This is the step
  that would turn the radial estimate into an actual nonlinear free-boundary
  solution, so it cannot be deferred as routine bookkeeping.

The importable implementation confirms the evidence boundary: its module
docstring explicitly says that it derives the conditional Grad--Shafranov
algebra and “does not solve the finite-R inverse or certify existence”
(`src/substrate_framework/euler_compact_ring.py:1-7`).  Its tests exercise
Euler residual and normalization identities, not the missing existence
theorem.

This reading also matches the actual external theorem scopes.  Feehan--Gong--
Song treats Feynman--Kac uniqueness for degenerate partial-boundary problems;
Bramanti--Zhu proves **local** Schauder/\(L^p\) estimates for operators
structured on Hörmander fields with drift; and Hamilton's inverse theorem
assumes a smooth tame inverse/right-inverse family.  None of those statements,
as cited, supplies this problem's moving characteristic-boundary Poisson map or
the bordered family that Hamilton requires.  This is an applicability finding,
not a demand for a new axiom or a duplicate oracle.

There is no counterexample here to equations (23), (26), or the candidate
branch.  The contradiction is between the proof's declared inputs and its
verdict: 0263's own source audit excludes an imported theorem for the
moving-boundary tame estimate, while `perturbation-gluing.md:126-133` declares
that estimate and finite-R existence established without supplying it.

## Reconciliation with review 0268

I opened 0268 only after reaching the conclusion above.  The prior review
repeats that the displayed intrinsic estimates and tame product bound provide
the uniform entrance resolvent (`0268/review.md:31-37`) and then says the
perturbative right inverse and scalar reduction produce the compact Euler field
(`0268/review.md:39-52`).  It does not define the missing spaces or losses,
write the bordered operator, derive a uniform parameter estimate, or check a
specific nonlinear inverse theorem.  It therefore is not an independent
closure of the gap.  Reviewer count and the accepted status do not change this
equation-level finding.

## Minimum repair for the existing finite-R route

The shortest positive repair is a finite-regularity Banach inverse theorem;
the proof need not begin with a full projective-limit Nash--Moser apparatus if
the asserted no-loss estimate is true.

1. Fix one reference disk and an explicit Hanzawa map.  Define weighted
   \(C^{k,\alpha}\) or intrinsic Sobolev spaces for the interior defining
   function, the free-boundary graph, source, edge trace, center gauge, and
   translation cokernel.  State how \(u=(\Phi/T^2)v\) embeds into the physical
   quadratic-form domain.
2. Write one bordered linear operator containing the PDE range equation, all
   free graph modes, the center gauge, the translation projection, and the
   profile border.  Prove that it is an isomorphism at \(\epsilon=0\).  The
   existing Fourier barriers and physical spectral gap are the principal
   ingredients, but the proof must include the collar/interior partition,
   trace and conormal compatibility, and uniform estimates for every mode.
3. Prove a neighborhood estimate for the exact pulled-back operator with one
   declared derivative count and coefficient norm.  If it is no-loss at a
   sufficiently large finite \(k\), use the ordinary Banach IFT and bootstrap.
   If loss is unavoidable, define a tame grading and smoothing operators and
   prove the smooth tame right-inverse family before invoking the named
   Nash--Moser theorem.
4. Only then solve the range equation, perform the already credible Hadamard
   division and scalar balance IFT, and check positivity, one-disk topology,
   axis separation, flat zero extension, pressure substitution, and quantitative
   \(C^4\) convergence on the fixed twist annulus.

This repair would restore the geometry conjunct without weakening the accepted
statement.  It also has an exposing failure test: the theorem must remain valid
after a small nonradial high-frequency boundary perturbation in the declared
coefficient norm; failure of the inverse bound or growth beyond the stated
loss identifies exactly where the present qualitative extension argument is
insufficient.

## Executable alternatives and the full nonlinear rotational-Euler objective

The geometry gap is upstream of the accepted prepared theorem, but closing it
still would not yield the user's full nonlinear objective.  `C-CST-018`
explicitly excludes an unrestricted Euler invariant manifold and nonlinear
finite-amplitude stability (`governance/claims.yaml:15370-15380`).  The current
construction prescribes linear Euler/Lin histories for each retained initial
amplitude; it does not prove that finite objective rotations and strains evolve
on a closed nonlinear material state space.

Three positive continuations remain concrete:

1. **Repair this ring and construct a nonlinear modulated manifold.**  After
   the bordered finite-R theorem is complete, parameterize finite translations
   and rotations by actual volume-preserving diffeomorphisms, with the rotation
   variable in \(SO(3)\), and substitute the resulting modulated compact-ring
   field into the full Euler action.  Compute the complete nonlinear residual,
   including pressure and intercell correctors.  Solve its normal component and
   prove a finite-window error theorem uniform over a fixed finite-amplitude
   neighborhood.  The reduced action must be frame-indifferent and its
   small-rotation Hessian must recover the accepted prepared coefficients.
   This is the direct route from the present geometry to a full nonlinear
   rotational law.
2. **Change the stationary supplier to a published compact toroidal shell.**
   The repository already records the Gavrilov/Constantin--La--Vicol smooth
   compact-shell existence route and Baldi's exact action-angle chart
   (`0251/source-transfer.md:17-30`).  This avoids the unproved custom
   Bessel-core free-boundary theorem, at the cost of redoing the same-field
   field-changing response, gain, and coefficient calculation on the annular
   shell.  It is especially attractive if the nonlinear target needs finite
   rotations and twist but not the exact Bessel central core.
3. **Retain exact projected dynamics with memory before claiming locality.**
   Project full Euler onto the finite rotation/strain observables and derive the
   exact normal-mode remainder or memory kernel.  A local nonlinear Cosserat or
   Vikulin law is earned only if that kernel is shown to reduce in a controlled
   scale-separated limit.  If it does not, the exact nonlocal equation is still
   a stronger positive result than assigning a local finite-rotation law from
   the linear prepared Hessian.

The first route is the minimum repair that preserves all of P251's chosen
geometry.  The second is the fastest independent geometry supplier.  The third
prevents the nonlinear objective from being silently replaced by another
prepared linear closure.

## Addendum: exact Bramanti--Zhu transfer and the spectral projection wording

I checked the primary [Bramanti--Zhu paper](https://msp.org/apde/2013/6-8/apde-v6-n8-p01-p.pdf),
not only 0261's source summary.  This sharpens the finding above.

The operator is locally of the right *algebraic* kind after the positive-side
problem has already been extended smoothly through \(T=0\).  Bramanti--Zhu
consider

\[
  L=\sum_{i,j=1}^q a_{ij}X_iX_j+a_0X_0
\]

with smooth real vector fields satisfying their weighted Hörmander condition,
uniformly positive \((a_{ij})\), and \(a_0\) bounded above and away from zero
(paper, Assumption (H), pp. 1797--1798).  For 0263's proposed extension,
\(a_{ij}=\delta_{ij}\), \(a_0=1\), and
\(X_1=T\partial_x,X_2=T\partial_z\).  At \(T=0\), \(X_0\) is a nonzero normal
field and \([X_0,X_j]=2|\nabla T|^2\partial_j\); on a compact coefficient
neighborhood with a uniform eikonal margin, this supplies the required finite
weighted bracket span.  A bounded zeroth-order \(c_T\) is not in the displayed
Bramanti--Zhu operator, but it is a standard local lower-order perturbation and
is not the load-bearing objection.

There is one literal hypothesis mismatch: the paper formulates its setting as
\(q+1<n\), while the meridional problem uses \(q=2\) diffusion fields on
\(n=2\) variables.  This means Theorem 2.1 is not directly applicable exactly
as cited.  I do **not** treat that dimensional inequality as a refutation: an
explicit auxiliary Rothschild--Stein-type lift or a theorem stated without the
redundant-field restriction may remove it.  The campaign neither constructs
such a lift nor cites that replacement, so direct applicability remains to be
written down.

More decisively, Theorem 2.1 is only a local second-intrinsic-order a priori
estimate on \(\Omega'\Subset\Omega\):

\[
 \|u\|_{C_X^{2,\alpha}(\Omega')}
 \le C\bigl(\|Lu\|_{C_X^\alpha(\Omega)}+\|u\|_{L^\infty(\Omega)}\bigr).
\]

That is exactly useful for a fixed two-sided extension and supports 0263's
base interior estimate.  It supplies no boundary trace or Poisson operator,
no solvability or global Fredholm inverse, and no parameter differentiability.
Most concretely, the authors say that extending the estimate to weighted order
\(k+2\) in terms of \(k\) derivatives of \(Lu\) is nontrivial and that they do
not address it (paper, p. 1797, lines immediately before Section 2).  Thus this
source cannot justify 0263's claims of “all intrinsic commutations,”
differentiation to every parameter order, or a smooth tame projective-limit
inverse.  Those may be provable by a new commutator induction for these special
smooth fields, but that proof is absent.  This makes the unsupported transfer
mathematical rather than merely a missing citation or a request for more
exposition.

The constant in Theorem 2.1 depends on the vector fields, the nested domains,
the ellipticity margin, and coefficient Hölder norms.  Uniformity over a
Hanzawa/profile family is plausible only after placing that family in one
fixed two-sided domain and proving uniform bracket, coefficient, and nesting
bounds.  Smooth or tame dependence on the parameter does not follow from
uniform constants alone.

Finally, `0263/perturbation-gluing.md:78-82` should be read in spectral-cluster
language.  At the straight solution the translation eigenspace is exact.  At a
general nonsolution Nash--Moser iterate, and for the curvature-broken radial
translation at finite \(\epsilon\), the corresponding small spectral cluster
and Riesz projection may persist while an exact kernel and cokernel do not.
That is repairable: form one bordered operator with the moving Riesz
projection, graph center condition, and profile balance coordinate, then prove
its inverse and the parameter derivatives.  Merely asserting persistence of
the “kernel and cokernel” does not construct the right inverse used in (4).
