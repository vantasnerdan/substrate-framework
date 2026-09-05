# Uniform measured third-K jets after sparse-tag normalization

This is the parent-diagonal continuation of the established 0250 fixed-ring
history result. It does not alter the finite-window gain construction or claim
an exact-in-K action normalizer.

## Normalized source and measured functional

Fix one sparse tag rung `N=N_j`, with

    chi=chi_0[1+epsilon sum_k c_k cos(N_k I)],   c_j=N_j^-j.

The exact 0260/0250 lift has, for every fixed derivative order,

    ||xi_N||_s+||v_N||_s <= C_s N^(s+1),

while its resonant physical angle, current and hybrid gains have measured
scale

    g_N=N c_N g_0(omega)[1+O(N^-1)],   |g_0|>=g_*>0

on the selected compact frequency band. The unit physical-output source is
therefore `(xi_hat_N,v_hat_N)=(N c_N)^-1(xi_N,v_N)`. Its raw graph/source
cost remains explicit:
`||xi_hat_N||_s+||v_hat_N||_s=O(N^(j+s))`.

For each retained physical output row, let `O_N(epsilon K,t)` be the exact
tag-plus-ambient observation from the defining material Fourier integral,
including moved centroids, absolute material velocities, pressure projection
and continuous ambient. Define

    J_N^(3)(K,t)=partial_epsilon^3 O_N(epsilon K,t)|_{epsilon=0}.

The bound below is on the measured jet and its finite-`K` parameter dependence,
not on a raw displacement norm.

## Sparse-sequence quantifier order

The tag is built inductively, rather than by choosing one exponent `j` after
the whole infinite sum is already fixed. At stage `j`, freeze the previous
carriers `N_1,...,N_{j-1}`, their measured constants, the finite
history/derivative list and the tolerance budgets for all earlier rungs. Then
choose `q_j>=j`, choose `N_j` larger than all prior scales and
resonance-separation thresholds, and set `c_j=N_j^-j`. Finally choose future
carriers recursively so fast that, for every already frozen rung `ell<=j`, the
future-tag contribution to each measured row is at most
`2^(-k-j)` times that row's remaining tolerance at stage `k`. At each new
stage only finitely many earlier inequalities must be met, so increasing
`N_k` makes them simultaneous; the same choice can enforce
`N_k^(m-k)->0` for every previously requested derivative order `m`, giving
one positive `C-infinity` tag.

The exponent test `q_j>=j` controls only the power of `N_j`. Its constant
`C_{q_j,r,T}` depends on the frozen sparse tag prefix, previously selected
carriers and the chosen finite derivative list; it is not asserted uniform in
`j`. It is included among the stage-`j` constants before `N_j` is increased.
The resonant constant in (1), by contrast, is a bound on a fixed compact
action band, fixed observation kernel and fixed finite gain-inverse class; it
does not grow with the carrier label. At stage `j` choose `N_j` large enough
that its normalized nonresonant tail is at most `2^-j` (and all earlier
tolerances hold). Thus the constructed sequence has one uniform resonant
bound plus a summable tail, even though the integration-by-parts constants
used to reach that bound are stage dependent.

## Resonant terms

Every degree-three coefficient is a finite sum of integrals

    integral chi(I) a_N(I,theta,t)
      x_{i_1}x_{i_2}x_{i_3} J(I)dI dtheta,

and first variations with one factor among `xi_N`, `D_t xi_N`, `v_N`, or a
derivative of `chi`. Evaluate these in material/Lagrangian coordinates before
estimating them: apparent products such as `chi' xi_N` are paired with the
moved-domain/Jacobian term and reduce to the single transported material
coefficient, so they do not create an artificial `N^2 c_N` factor. On the
fixed compact ring all spatial and background coefficients are bounded. The
lift has only one carrier-sized factor: `xi_N` and `D_t xi_N` contribute
`O(N)`, while exact pressure and lower material terms contribute `O(1)`. The
resonant tag coefficient contributes `c_N`; its differentiated Fourier
coefficient is already included in that same single material factor. Hence,
for every finite
time/derivative list `r`,

    |J_{N,res}^(3)|_{C^r([-T,T])} <= C_{3,r,T} N c_N.       (1)

This includes the complete 0241 second/third material-momentum tensors and
centroid phase. The ambient term is handled by the oscillatory estimate next,
not by inserting `||xi_N||_s`.

## Oscillatory tails

For a smooth compact action amplitude, integration by parts gives

    integral a(I,theta,t)e^{iNI}dI
      =(iN)^-q integral partial_I^q a(I,theta,t)e^{iNI}dI. (2)

There is no boundary term. The same identity applies to every nonresonant tag
harmonic. The sparse C-infinity tag makes all differentiated Fourier
coefficients summable. After removing the single resonant conjugate `N_j`
term, angular Fourier selection followed by (2) gives, for every finite `q,r`,

    |J_{N,tail}^(3)|_{C^r([-T,T])} <= C_{3,q,r,T} N^-q.     (3)

The whole-space Leray/pressure contribution has the same finite-order
oscillatory representation; its noncompact tail is included in the accepted
finite moment constants, not estimated by a growing displacement norm.

Combining the unit-gain normalization with (1)--(3) gives

    |J_hat_N^(3)|_{C^r}
      <= C_{3,r,T}+C_{3,q,r,T} N^(j-q-1).                 (4)

Choosing `q_j>=j` makes the normalized complete measured third jets bounded
for that rung, with the tail tending to zero as its already-frozen `N_j` is
increased. The bound is uniform in the finite `K`-ball as follows.

## Uniform third jets on a `K`-ball

Let `K_0` be the fixed long-wave neighborhood on which the finite physical
gain inverse and its parameter-dependent source amplitudes are `C^3` with
bounded derivatives. Direct differentiation of the defining material Fourier
integral gives, for `|alpha|=3`,

    partial_K^alpha integral chi U exp(-iK.r)dx
      =integral chi U (-i r)^alpha exp(-iK.r)dx.          (5)

The ring radius bounds `|r|`, while the complete material variation uses the
same Lagrangian formula as above; the continuous ambient and whole-space
Leray projector contribute their accepted finite third moments. The projector
is a fixed spatial operator, independent of macroscopic `K`, so it introduces
no `K` singularity. Differentiating the actual `C^3` gain-inverse source map
adds only its bounded first-through-third parameter derivatives. Therefore,
for every finite `T,r`,

    sup_{j:N_j>=N_*} sup_|K|<=K_0
      |D_K^3 O_hat_{N_j}(K,t)|_{C^r([-T,T])}
      <= C_{3,r,T,K_0}.                                   (6)

The same bound applies to the reflected whole-law average. This is a bound on
the full third derivative throughout the ball, not merely on
`D_K^3 O_hat_N(0,t)`. Taylor's integral remainder consequently gives

    |O_hat_N(K)-O_hat_N(0)-D O_hat_N(0)K
      -1/2 D^2O_hat_N(0)[K,K]|
      <= C_{3,r,T,K_0}|K|^3/6.                            (7)

## Parent K/N diagonal

Let `C_3` be the finite bound in (6), after the finite observation/gain
inverse and requested time derivatives have been selected. Choose the
macroscopic wave number from this measured constant first:

    C_3 |K| <= epsilon/3.                                  (8)

The third Taylor remainder in (7) is then at most `epsilon |K|^2/3`; no
relation `K=N^-D` is required. With this fixed `K`, choose the sparse rung and
the full-pressure parametrix order `M` so that the exact-evolution-minus-
parametrix error, after unit-gain normalization, obeys

    E_{M,r,T} N^(j-M-1) <= epsilon |K|^2/3.                (9)

Taking `M>j+1` and then increasing the selected carrier makes (9) hold. The
full WKB parametrix itself differs from its leading clock/observation by a
separate normalized subleading term

    |O_hat_N^param-O_hat_N^lead| <= C_clock,r,T/N.          (10)

This is not the exact-evolution-minus-parametrix error. After `K` has been
chosen from (8), increase the carrier so that
`C_clock/N <= epsilon |K|^2/3`; then choose `M` for (9). The finite
band-synthesis coefficient sum and derivative graph cost `O(N^(j+s))` are
constants at this stage and remain in the preparation ledger. Remaining
third-order and initial-map terms receive the last third of the error budget.

Thus the established finite-window result survives the parent diagonal:
measured third jets are uniformly bounded after resonance normalization,
while source/action costs remain explicit. The exact-in-K full-form normalizer
under development may improve this ordering, but it is not needed to compare
an unnormalized `O(1/N)` remainder directly with `K^2`.

`route_verdict: normalized measured third-K jets uniformly bounded and parent
K/N diagonal established at the stated finite-window observation scope`

`remaining scope: exact-in-K full action/Jacobi normalization, final periodic
0145/0147 transfer and independent geometry/density remain separate`
