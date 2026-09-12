# P253/0024: reduced monodromy powers and modulated normal form

## Preregistration and activation boundary

This README alone preregisters attempt 0024. No new primary-source body,
operator derivation, spectral calculation, verifier, or numerical result may
be opened or produced for 0024 until root:

1. appends 0024 to the central P253 obligation and candidate ledger;
2. records ownership, dependencies, and the non-inheritance boundary; and
3. places a repository-schema exit of zero in this directory.

After activation, 0024 owns only append-only files in this directory. Root
owns central proposal, source/API/test, memory, commit, and promotion changes.
The accepted base is release `v0.183.0`, whose manifest pins
`substrate-framework@b39f392f8feee21c93a6052a5f971bc2df862779`.

Target kind: a fixed linear Hamiltonian/Floquet theorem for powers of the
actual P253/0021 reduced two-contour return map, with competing analytic
representations. The target is not a new periodic solution and not a
finite-dimensional center-history surrogate.

No empirical comparator is involved. No production numerical verifier is
registered at freeze. Any later calculation of a soft multiplier, Krein
collision, spectral edge, small action splitting, or near-unit resolvent
first requires the small-ratio-numerics skill and a separately frozen error
budget. Long production checks, if later licensed, run only in the root-
coordinated scripts pane.

## Frozen object and exact claim hierarchy

The candidate input from P253/0021 is the fixed-source-orbit evolution

\[
 U(t,r):X^s\longrightarrow X^s,
 \qquad M=U(\mathcal T,0),                                   \tag{1}
\]

on two independent, non-reversible, individual-area contour tangents. Its
fixed-level and symmetry flag is

\[
 G=\operatorname{span}\{v_{\rm ph},v_z\}
 \subset K=\ker d\mathscr H_c\cap\ker dP\subset X^s,         \tag{2}
\]

and the reduced space and return map are

\[
 \mathcal N_{\rm red}=K/G,
 \qquad M_{\rm red}:\mathcal N_{\rm red}\to\mathcal N_{\rm red}.
                                                                    \tag{3}
\]

The symbols in (2) use the normalized 0015 convention; multiplying
\(\mathscr H_c,P,\Omega\), and the action by the common positive physical
factor \(2\pi\rho_m\) changes neither (2), the Hamiltonian vector field, nor
the reduced dynamics.

The exact primary success target is

\[
 \boxed{\sup_{n\ge0}\|M_{\rm red}^n\|_{\mathcal N^s_{\rm red}
             \to\mathcal N^s_{\rm red}}<\infty}              \tag{4}
\]

in a declared Hilbert/Sobolev topology equivalent to the 0021 contour norm.
The following conclusions are strictly ordered and cannot inherit from one
another:

1. one-period boundedness of \(M_{\rm red}\) — already the 0021 candidate
   result, but insufficient for (4);
2. spectral containment in the unit circle — insufficient for (4) in
   infinite dimension and insufficient at a Jordan point;
3. power boundedness (4), earned by an equivalent positive invariant metric,
   exact unitary similarity, or complete spectral/Riesz control;
4. a quantified polynomial normal form when (4) fails only through named
   finite nilpotent shears;
5. nonlinear orbital stability — outside 0024 unless a separate nonlinear
   iteration and modulation theorem is activated.

The full, non-reduced monodromy is expected to retain phase/translation and
conserved-level companions. In an 0021 splitting it has the flag-compatible
form

\[
 M=\begin{pmatrix}
 M_{\rm red}&0&B\\
 C&I_2&D\\
 0&0&I_2
 \end{pmatrix}.                                               \tag{5}
\]

The companion coordinates are physical nearby-energy/impulse directions,
not gauge. A full-map estimate must modulate them or retain their shear; (4)
does not silently assert \(\sup_n\|M^n\|<\infty\).

## Frozen source and dependency inventory

Only the following already exposed artifacts may be used before activation:

| Input | Exact role | Scope limit |
|---|---|---|
| García–Hassainia–Hmidi, `arXiv:2603.21644v1`, cached PDF SHA-256 `7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6` | Actual axisymmetric no-swirl translating-frame-periodic two-patch family, exact contour equation, and restricted periodic operator estimates | No full non-reversible Floquet spectrum, monodromy power bound, nonlinear stability, or patch-family derivative across the Cantor set |
| P253/0015 corrected construction, current SHA-256 `90eec03cd419ca92e6796e6e4d3839b6c6d2152ffa2b04ef49048eacd3a5bff9` | Actual labeled leaf, normalized/physical Hamiltonian and KKS conventions, local critical-loop action, Jacobi equation, and symmetry modes | Raw action is indefinite; no restoring norm or bounded patch monodromy powers |
| P253/0021 construction, SHA-256 `d1c27273e34903e76126233076fdcb87cbe41f47ffb37bb8c8700b7087c301aa` | Candidate full two-label IVP propagator, symplectic monodromy, exact invariant flag, companions, and bounded reduced one-period map | Proposal evidence pending its own review; no spectral, power-bounded, or nonlinear conclusion |
| P253/0021 source receipt, SHA-256 `d09dd90b2f3c70ca953464bcca583dc17d306d73405c9fc169b825d9da54e143` | Page-level source domain and operator inventory | Does not broaden the source restrictions |
| Accepted `C-FLO-001` | Exact finite-dimensional criterion: unit-circle eigenvalues require semisimplicity for bounded powers | No unbounded-operator domain, Riesz basis, essential spectrum, contour flow, or infinite-dimensional sufficiency |

No additional external source is selected at freeze. If activation shows that
an infinite-dimensional symplectic/Krein, Riesz-basis, exponential-dichotomy,
or evolution theorem is necessary, record its title, version, URL, hash,
exact operator domain, and proposed implication before opening its body.

### Exact García source scope to preserve

For fixed \(0<a<b\), \(\kappa>0\), sufficiently small \(\epsilon\), and
\(\lambda\) in the source Borel Cantor set \(\mathcal C_\epsilon\), the source
constructs one sufficiently smooth finite-core Euler patch pair. It is
periodic only after a common axial translation, and its two labels obey a
half-period relation in the existence chart.

The source periodic linearized operator acts on a space imposing
zero area, reversibility, half-period label reduction, and first-mode
conditions. On the higher-mode block its audited normal form is

\[
\begin{aligned}
 \Phi^{-1}\mathcal L_\perp\Phi
 ={}&\epsilon^2L\omega\partial_\phi
   +c_1\partial_\theta+c_2\mathcal H
   +\epsilon^2c_3\partial_\theta\Lambda_2\\
 &+\epsilon^2L^{1/2}\Pi^\perp\partial_\theta
      R_{4,-3}\Pi^\perp+E_{4,n},                            \tag{6}
\end{aligned}
\]

with \(L=|\log\epsilon|\), \(c_2=1/2+O(\epsilon^2)\), smoothing
\(R_{4,-3}\), finite-step error \(E_{4,n}\), and diagonal normal frequencies

\[
 \mu_{j,2}=c_1j+c_2\operatorname{sign}j+\epsilon^2c_3d_j,
 \qquad d_j=-\frac{j}{4(j^2-1)|j|},\quad |j|\ge2.            \tag{7}
\]

Equations (6)-(7) are estimates for a periodic boundary-value/invariance
operator on the source-restricted space. They are not yet an IVP
conjugation, do not include two independent labels or all first modes, and do
not prove (4). A periodic inverse can exclude a restricted periodic kernel;
it cannot control powers of a return map.

## Common operator domain and required reductions

After activation, every route works first on the actual 0021 fixed model

\[
 X^s=H^s_0(\mathbb T)\oplus H^s_0(\mathbb T),
 \qquad A_*(t):X^{s+1}\to X^s,                              \tag{8}
\]

or on an explicitly proved equivalent tame scale. Both labels are
independent; all first modes remain until the physical center, phase, and
axial modulation is performed. Tangential relabeling alone is gauge.

The quotient (3) must be implemented by bounded phase and axial-centering
conditions. Conserved-level companions \(w_E,w_P\) are retained with

\[
 d\mathscr H_c(w_E)=1,\ dP(w_E)=0,
 \qquad d\mathscr H_c(w_P)=0,\ dP(w_P)=1.                   \tag{9}
\]

Any source conjugation must be lifted through this complete decomposition.
It may not delete \(B,C,D\) in (5), impose reversibility, identify the two
labels, or call a first mode gauge merely because the source removed it for
solvability.

## Route A — full IVP reducibility and positive Krein metric

Lift the source transformations to periodic strongly continuous operators
\(\Phi(t)\) on the full two-label IVP and prove

\[
 \Phi(t)^{-1}(\partial_t-A_*(t))\Phi(t)
   =\partial_t-\bigl(D_*+R_*(t)\bigr),                       \tag{10}
\]

with exact domains, derivative losses, symplecticity, bounds on
\(\Phi^{\pm1}\), and all center/shape cross blocks stated. The diagonal
part must be skew-adjoint in a positive equivalent metric, not merely have
imaginary formal frequencies.

Route A earns (4) only by one of the following complete closures:

- exact reducibility \(R_*=0\) and a bounded periodic similarity to a unitary
  group;
- a positive periodic metric \(G(t)\) satisfying
  \(\dot G+A_*^*G+GA_*=0\), with
  \(mI\le G(t)\le MI\); or
- a full Krein/Riesz theorem proving that every residual coupling is
  spectrally stable, semisimple at the unit circle, and has uniformly bounded
  spectral projections.

A small nonzero one-period remainder alone gives at best
\(\|M_{\rm red}^n\|\lesssim(1+\delta)^n\) and therefore cannot establish
(4). Source Cantor nonresonance denominators must be translated into the
actual evolution domain rather than copied from the periodic inverse.

**Route-A success:** an explicit similarity or positive invariant metric and
a finite constant \(C\) with \(\|M_{\rm red}^n\|\le C\) for every integer
\(n\ge0\).

## Route B — direct reduced resolvent and Krein-signature construction

Without importing the source boundary-value conjugation as an IVP theorem,
start from the symplectic \(M_{\rm red}\) constructed in 0021. Determine its
essential spectrum from the self-interaction principal operator, and retain
the smooth cross interaction as a compact operator only after proving the
relevant Sobolev compactness.

For each isolated spectral cluster near the unit circle, construct the Riesz
projection

\[
 \Pi_\Gamma=\frac1{2\pi i}\oint_\Gamma(z-M_{\rm red})^{-1}\,dz, \tag{11}
\]

prove a resolvent bound uniform over clusters and high angular modes, compute
the KKS/Krein signature on \(\operatorname{ran}\Pi_\Gamma\), and exclude both
opposite-sign collisions and nontrivial Jordan chains. Completeness and
uniform conditioning of the spectral subspaces must be proved; a list of
unit-modulus multipliers is not enough in infinite dimension.

**Route-B success:** a uniformly conditioned Riesz decomposition, unit-circle
spectral location with definite signatures, and semisimplicity sufficient to
sum the modal bounds and obtain (4).

Route B is materially different from Route A: it analyzes the already
defined return operator and its resolvent rather than converting the source
periodic normal-form machinery into an evolution conjugation.

## Route C — exact modulated polynomial normal form

If A or B exposes a genuine unit-multiplier Jordan obstruction rather than a
missing estimate, retain it and change the target representation. Construct
bounded spectral/modulation projections and a similarity

\[
 S^{-1}M S=U\oplus\bigoplus_{k=1}^N(I+N_k),
 \qquad U^*U=I,\qquad N_k^{m_k+1}=0,                         \tag{12}
\]

on the reduced-plus-modulated physical space. Then derive exactly

\[
 \|M^n\|\le C(1+n^{m_*}),
 \qquad m_*=\max_k m_k,                                     \tag{13}
\]

and identify which terms are phase twist, axial drift, energy/impulse
companion shear, or genuine shape-sector growth. The source filament twist
may motivate a block but cannot be transferred to the finite-core patch
without an actual differentiable family or a directly constructed
generalized Jacobi vector.

**Route-C success:** the strongest exact polynomial-growth normal form on the
actual full/modulated patch tangent, with every nilpotent block physically
identified. This is a constructive fallback if (4) is refuted by a named
Jordan mechanism; it is not relabeled power boundedness.

## Frozen route comparison and continuation ladder

Selection criteria, fixed before activation, are:

1. exact coverage of the actual finite-core pair and full non-reversible,
   independent-label perturbation class;
2. correct operator domains, derivative losses, and periodic similarities;
3. preservation of KKS symplecticity, phase/translation reduction, and
   energy/impulse companions;
4. proof of infinite-dimensional completeness/conditioning rather than a
   finite sampled spectrum;
5. explicit dependence on \(\epsilon,\lambda,\kappa,s\) and source Cantor
   constants;
6. a statement about all integer powers, not one period or a finite window;
7. the strongest route-scoped result with no inheritance to nonlinear or
   particle claims.

Route A is preferred if the source transformations extend to the full IVP
with a positive metric. Route B is preferred if direct compact-perturbation
and resolvent estimates expose the spectrum more cleanly. If either route
fails, execute method repair and then the other representation in the same
attempt. A proved Jordan obstruction activates Route C. A missing estimate is
recorded as that route's missing construction and does not become a Jordan
obstruction or obligation-level no-go.

Each tested route receives exactly one verdict: `established`, `refuted` with
the exact mechanism, or `blocked` with the missing construction. The attempt
must preserve the strongest positive theorem and name its next executable
dependency.

## Evidence and non-inheritance contract

After activation, required artifacts are:

- source/interface receipt with exact pages, versions, hashes, operator
  domains, tame losses, Cantor quantifiers, and what does not transfer;
- full two-label IVP conjugation or direct monodromy resolvent construction;
- exact phase/translation and energy/impulse reduction, including block (5);
- a positive-metric/unitary-similarity proof, complete Krein/Riesz proof, or
  exact polynomial normal form;
- symbolic checks for tractable conjugation, symplectic, projection, and
  Jordan-power identities;
- one route verdict per tested route and the next materially different
  construction.

Maximum 0024 verdict is reduced linear power boundedness (4), or the exact
polynomial-growth statement (13) if a named obstruction refutes (4). Neither
licenses nonlinear orbital stability, an open Euler neighborhood for all
time, general three-dimensional or swirl stability, helicity, intrinsic or
quantum spin, electron or neutrino identification, parent completion, or a
global no-go.

## Activation handoff

Root should centrally register 0024 and place the activation-schema receipts
in this directory. Substantive work begins only after
`activation-schema.exit` records zero.
