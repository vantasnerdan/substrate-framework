# P253/0021: full contour monodromy or controlled filament-to-patch Floquet transfer

## Preregistration and activation boundary

This README alone preregisters attempt 0021. No new primary-source body,
derivation, verifier, spectrum, or numerical result may be opened or produced
for 0021 until root:

1. appends 0021 to the central P253 candidate/obligation ledger;
2. records the exact attempt ownership and non-inheritance boundary; and
3. places a repository-schema exit of zero in this directory.

After activation, 0021 owns only its append-only attempt artifacts. Root owns
central proposal/API edits and commits. The target is the actual
García–Hassainia–Hmidi labeled periodic patch pair from 0015, not a prescribed
center history, two independently evolved rings, or a finite-dimensional
surrogate.

Target kind: fixed periodic-orbit linear evolution theorem with two competing
representations. Route A constructs the full two-contour propagator and
reduced monodromy directly. Route B may transfer the exact source filament
return map to the patch only through a proved center/shape decomposition and
operator-norm remainder bound.

No empirical comparator is involved. No numerical verifier is frozen here.
Any later soft Floquet multiplier, spectral edge, or action splitting first
requires the small-ratio-numerics skill and a separate analytic-closure
receipt.

## Frozen inputs before any new body

Only already audited artifacts may be read before activation:

| Input | Frozen role | What it does not provide |
|---|---|---|
| García–Hassainia–Hmidi, arXiv:2603.21644v1, cached PDF SHA-256 7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6 | Actual axisymmetric no-swirl two-patch periodic family, source contour equation, periodic Sobolev chart, mode reductions, and Nash–Moser inverse estimates already audited in 0015 | No full non-reversible initial-value propagator, full period map, Floquet stability, or nonlinear orbital stability |
| P253/0015 source receipt and construction | Exact translating frame, labeled product leaf, \(E-cP\), corrected KKS action, full Euler Jacobi equation, symmetry modes, and source-scope audit | Its period map was conditional on a bounded evolution realization |
| P253/0015 action-sign correction | Shared sign convention \(i_{X_H}\Omega=dH\), \(-d\Theta=\Omega\) | No stability conclusion |
| P253/0014 and the corrected 0005/0011 action audit | Independent sign check and warning that Hamiltonian Jordan shear is not restoring coercivity | Different carrier/sector; no theorem about the García pair |

No new external source is selected at freeze. If activation reveals that a
well-posedness, Hamiltonian PDE, or Floquet theorem not already in this
inventory is necessary, record its title, version, exact proposed role, and
scope limit before opening its body.

## Actual source estimates and their exact domain

The following are frozen from the completed 0015 source audit, not newly
inferred from a source body.

The source uses periodic phase \(\phi\) and contour angle \(\theta\), with
\(3<s_0\le s\le s_{\rm up}\), \(1<\tau\le3/2\), and Lipschitz parameter
spaces over an interval \(\mathcal O\). Its nonlinear existence chart is

\[
 X_s=\operatorname{Lip}_\nu
 \bigl(\mathcal O,H^s_{\circ,\mathrm{even}}(\mathbb T^2)\bigr),
 \qquad
 Y_{s-1}=\operatorname{Lip}_\nu
 \bigl(\mathcal O,H^{s-1}_{\circ,\mathrm{odd}}(\mathbb T^2)\bigr).
                                                               \tag{1}
\]

Here the source spaces impose zero spatial average, reversibility, a first
sine-mode condition, and a global first-cosine normalization. The normal
space \(X_s^\perp\) additionally removes the first cosine mode, so its spatial
Fourier modes satisfy \(|j|\ge2\). Label 2 is not independent in this chart:
it is label 1 shifted by half a period.

At the renormalized state \(g\), Corollary 7.2 gives a periodic invariance
operator

\[
 \mathcal L(\epsilon,g)=\partial_gF(\epsilon,g):
 X_s\longrightarrow Y_{s-1},                                  \tag{2}
\]

whose leading terms are

\[
 \epsilon^2|\log\epsilon|\,\omega(\lambda)\partial_\phi
 +\frac12\mathcal H-\partial_\theta,                           \tag{3}
\]

plus the source's variable transport, heterogeneous Hilbert, shifting,
first-mode modulation, logarithmic integral, and smoothing terms. These
terms must all be retained after activation; equation (3) is not a replacement
operator.

For the normal block, Proposition 8.3 and Sections 9.1–9.4 construct
periodic, reversibility-preserving changes of variables \(\Phi\) satisfying

\[
 \|\Phi^{\pm1}h\|^{\operatorname{Lip},\nu}_s
 \lesssim
 \|h\|^{\operatorname{Lip},\nu}_s
 +\|g\|^{\operatorname{Lip},\nu}_{s+3}
   \|h\|^{\operatorname{Lip},\nu}_{s_0},                       \tag{4}
\]

and, on the stated Cantor restrictions,

\[
\begin{aligned}
 \Phi^{-1}\mathcal L_\perp\Phi
 ={}&\epsilon^2|\log\epsilon|\,\omega\partial_\phi
     +c_1\partial_\theta+c_2\mathcal H
     +\epsilon^2c_3\partial_\theta\Lambda_2\\
 &+\epsilon^2|\log\epsilon|^{1/2}
   \Pi^\perp\partial_\theta R_{4,-3}\Pi^\perp+E_{4,n}.          \tag{5}
\end{aligned}
\]

The audited estimates include

\[
 c_2=\frac12+O(\epsilon^2),\qquad
 \|\partial_\theta^4R_{4,-3}h\|_s
 \lesssim\|h\|_s+\|g\|_{s+2\tau+6}\|h\|_{s_0},                 \tag{6}
\]

and

\[
\|E_{4,n}h\|_{s_0}
\lesssim
\epsilon^5N_0^{\mu_2}N_{n+1}^{-\mu_2}\|h\|_{s_0+2}
+\epsilon N_n^{s_0-s}\|g\|_{s+2\tau+6}\|h\|_{s_0},            \tag{7}
\]

with \(\mu_2\ge4\tau+3\). The diagonal normal frequencies are

\[
 \mu_{j,2}=c_1j+c_2\operatorname{sign}j
 +\epsilon^2c_3d_j,\qquad
 d_j=-\frac{j}{4(j^2-1)|j|},\qquad |j|\ge2.                   \tag{8}
\]

The source approximate inverse loses angular derivatives and carries the
small-divisor factor \(\nu^{-1}\); the first-cosine block has an inverse of
size \(O((\epsilon^2|\log\epsilon|)^{-1})\) in its restricted periodic
space. Exact loss exponents, Cantor sets, and finite-step errors are to be
copied from the pinned pages after activation before any transfer claim.

These are estimates for a periodic boundary-value/invariance operator on
(1), with half-period label reduction and reversibility. They are not
estimates for an initial-value evolution family \(U(t,s)\). Even an exact
inverse of \(\partial_t-A(t)\) on periodic histories would at most exclude
multiplier \(1\) on its stated restricted domain; the source supplies an
approximate tame inverse during Nash–Moser, not control of
\(\sup_n\|\mathcal M^n\|\).

## Corrected Hamiltonian/action convention

In the translating frame,

\[
 \mathscr H_c=E(q_1+q_2)-c_\epsilon P(q_1+q_2),\qquad
 \dot q=X_{\mathscr H_c}(q).
\]

For \(\delta_\chi q_j=\{q_j,\chi_j\}\), freeze

\[
 \Omega_q(\delta_\chi q,\delta_\psi q)
 =-\sum_{j=1}^2\int q_j\{\chi_j,\psi_j\}\,d\rho\,dz,
 \qquad i_{X_{\mathscr H_c}}\Omega=d\mathscr H_c.              \tag{9}
\]

A local potential obeys

\[
 -d\Theta=\Omega,
\]

so the correctly signed local action is

\[
 \mathscr A[q]=\int_0^{\mathcal T}
 \left[\Theta_q(\dot q)-\mathscr H_c(q)\right]dt,              \tag{10}
\]

or \(-\int_\Sigma\Omega-\int\mathscr H_cdt\) in symplectic-area
form. In a linear symplectic chart its quadratic part is

\[
 \boxed{\mathscr A^{(2)}[\eta]
 =\int_0^{\mathcal T}
 \left[-\frac12\Omega(\eta,\dot\eta)
       -\frac12D^2\mathscr H_c(q_*)[\eta,\eta]\right]dt.}       \tag{11}
\]

Its stationarity gives the actual Jacobi equation. The sign in (11) must be
checked independently against the direct contour linearization; flipping
only the kinetic term is not an admissible repair.

## Full physical contour domain to construct

The source existence chart is too narrow for orbital perturbations. After
activation, build a local contour chart around the two actual source
boundaries with:

- two independent labeled normal graphs;
- exact zero-area constraint for each patch;
- all physical first Fourier modes retained;
- tangential reparametrizations quotient out;
- no reversibility restriction;
- no imposed relation \(h_2(t)=h_1(t+\mathcal T/2)\).

One candidate chart separates centers and shapes:

\[
 \mathcal X^s_{\rm full}
 =\mathbb R^4_{\rm ctr}\oplus
   H^s_{\ge2,0}(\mathbb T)\oplus H^s_{\ge2,0}(\mathbb T),      \tag{12}
\]

where the two shape factors have zero area and no first modes. The
center/shape coordinate map must be derived from the source radial-square
boundary variable and proved locally invertible; (12) is not assumed merely
because Fourier projections exist.

The intended evolution domain is

\[
 C([0,\mathcal T],\mathcal X^s_{\rm full})
 \cap C^1([0,\mathcal T],\mathcal X^{s-1}_{\rm full}),         \tag{13}
\]

or the exact tame scale dictated by derivative loss. Linearizing the complete
two-contour equation must yield

\[
 \partial_t\eta=A_*(t)\eta,\qquad
 A_*(t+\mathcal T)=A_*(t),\qquad
 U(t,s):\mathcal X^s_{\rm full}\to\mathcal X^{s-\sigma}_{\rm full},
                                                               \tag{14}
\]

with a proved loss \(\sigma\), composition law, and one-period bound. Only
then define

\[
 \mathcal M=U(\mathcal T,0).                                  \tag{15}
\]

## Exact symmetry, momentum, and period-twist reduction

The two Hamiltonian symmetry fields are

\[
 v_{\rm ph}=X_{\mathscr H_c}(q_*)=\dot q_*,
 \qquad
 v_z=X_P(q_*)=-\partial_zq_*.                                 \tag{16}
\]

They obey \(i_{v_{\rm ph}}\Omega=d\mathscr H_c\) and
\(i_{v_z}\Omega=dP\). The regular reduced tangent is therefore

\[
 \boxed{
 \mathcal N_{\rm red}
 =\frac{\ker d\mathscr H_c(q_*)\cap\ker dP(q_*)}
        {\operatorname{span}\{v_{\rm ph},v_z\}}.}              \tag{17}
\]

A Poincaré phase condition and an axial-centering condition choose a
representative of this quotient. They must not delete physical \(|j|\ge2\)
shape modes.

The quotient alone is insufficient until the conserved-level companions are
accounted for. Construct \(w_E,w_P\) so that the \(2\times2\) matrix

\[
 \begin{pmatrix}
 d\mathscr H_c(w_E)&d\mathscr H_c(w_P)\\
 dP(w_E)&dP(w_P)
 \end{pmatrix}                                                \tag{18}
\]

is nonsingular and normalize it to the identity. These companions carry
changes of conserved energy and impulse; they are not discarded as gauge.
For a relative-periodic family satisfying

\[
 \tau_{-a(\lambda)}
 \Phi^{\mathcal T(\lambda)}(q_\lambda)=q_\lambda,
\]

with \(\tau_aq(\rho,z)=q(\rho,z-a)\), differentiation gives the
convention-fixed generalized relation

\[
 (\mathcal M_{\rm rel}-I)\partial_\lambda q_\lambda
 =-\mathcal T'(\lambda)v_{\rm ph}
   +a'(\lambda)v_z.                                           \tag{19}
\]

The source filament theorem has \(\mathcal T'(\lambda)\ne0\), so its phase
Jordan shear is real. Variation of the common radial momentum parameter
\(\kappa\), or an independently constructed leaf direction, is needed for
the impulse companion. Because the full patch family exists only on a Cantor
set, differentiability of \(q_\lambda\) and (19) for the patch cannot be
inferred from source existence; it must be proved or kept as a filament-only
input.

For fixed-level spectral analysis use (17). For orbital estimates allowing
nearby conserved levels, modulate phase, axial shift, and the two companion
coordinates and prove their drift remains tangent to a nearby orbit tube.
A phase quotient by itself cannot remove the nonsemisimple period-twist
block.

## Route A — actual full-contour propagator and monodromy

Route A derives (12)–(15) directly from source contour equation (2.10).
The analytic ladder is:

1. Derive the exact two-label normal linearization, including self and cross
   kernels, center/shape coupling, speed terms, and the area/tangential gauge.
2. Derive the same operator from the quadratic action (11); equality fixes
   all KKS, time, translation, and kernel signs.
3. Prove tame initial-value well-posedness through one full period on (13),
   with every derivative loss and dependence on \(\epsilon,\nu,s\) explicit.
4. Prove symplectic transport
   \(U(t,0)^*\Omega_{q_*(t)}U(t,0)=\Omega_{q_*(0)}\), then define
   the bounded monodromy (15).
5. Construct the four-dimensional phase/translation/companion block and the
   reduced map on (17). State algebraic and geometric multiplicities at
   multiplier \(1\); do not silently assume semisimplicity.
6. Prove the strongest available reduced conclusion: bounded one-period map,
   spectral location with Krein signatures, power boundedness, or a
   quantified growth bound. These are distinct verdict strengths.

Route-A success requires an actual bounded evolution/return operator on the
full non-reversible two-label domain. Repackaging the source periodic inverse
as \(\mathcal M\), or proving only that a restricted periodic kernel is zero,
does not pass.

## Route B — controlled filament-to-patch Floquet transfer

Route B begins only from the locally invertible center/shape chart in (12).
Write the full one-period map in blocks,

\[
 \mathcal M_{\rm patch}=
 \begin{pmatrix}
 M_{\rm ctr}&M_{\rm ctr,sh}\\
 M_{\rm sh,ctr}&M_{\rm sh}
 \end{pmatrix}.                                                \tag{20}
\]

The comparison object is the exact source filament return map
\(M_{\rm fil}\), including its period-twist Jordan chain, together with the
shape evolution obtained from the full time-dependent contour generator.
The required transfer estimate has the form

\[
 \left\|
 \Pi_{\rm red}\left[
 \mathcal M_{\rm patch}
 -(M_{\rm fil}\oplus M_{\rm sh}^{(0)})
 \right]\Pi_{\rm red}
 \right\|_{\mathcal X^s\to\mathcal X^{s-\sigma}}
 \le C\epsilon^\alpha|\log\epsilon|^\beta,                    \tag{21}
\]

with derived \(\alpha>0,\beta,\sigma,C\), not a formal
small-core statement. Cross blocks in (20) must be estimated, not set to
zero.

Because \(M_{\rm fil}\) is nonsemisimple at multiplier \(1\), ordinary
eigenvalue perturbation or a Bauer–Fike estimate on the unreduced block is
invalid as a restoring argument. First extract (19), fix the conserved-level
companions, and compare the reduced block. The source frequencies (8) can
support a transfer only after the periodic conjugations in (5) are upgraded
to an initial-value transformation and the smoothing remainder is controlled
over repeated periods.

Route-B success requires (21) plus a separation/Krein or resolvent estimate
strong enough to survive its right-hand side. If it proves only a one-period
approximation, the maximum verdict is finite-time Floquet transfer, not
power-bounded or nonlinear stability.

## Frozen route comparison

Selection criteria are fixed before activation:

1. exact coverage of the actual labeled patch pair and arbitrary nearby
   source-regular contour data;
2. correct KKS/action signs and exact self/cross interaction;
3. no reversibility, half-period label, or first-mode restriction hidden in
   the physical perturbation class;
4. explicit phase/translation quotient and energy/impulse companions;
5. quantified derivative loss and \(\epsilon,\nu\) dependence;
6. a bounded initial-value operator rather than a periodic solvability
   surrogate;
7. honest separation of one-period, spectral, power-bounded, and nonlinear
   conclusions.

Route A is preferred if the source contour calculus yields a tame propagator
directly. Route B is preferred if the center/shape scaling supplies a sharper
operator bound while retaining all cross blocks. Failure of either route
activates a representation change to a reduced resolvent/Krein estimate for
the same monodromy; it does not close the restoring obligation.

## Verdict and evidence contract after activation

Each tested route receives exactly one verdict:

- established only for the complete statement named by that route;
- refuted only with an exact ill-posedness, spectral, or transfer obstruction
  and its mechanism;
- blocked with the missing propagator, reduction, or norm estimate named.

Required artifacts after activation are:

- source-estimate receipt with exact pages, hypotheses, domains, losses, and
  Cantor quantifiers;
- full two-contour chart and tangent/gauge proof;
- direct and action-derived linearizations with an independent sign check;
- propagator and monodromy construction or the exact obstruction to it;
- phase/translation/energy/impulse block and reduced-space definition;
- controlled transfer bound if Route B is used;
- exact symbolic checks for tractable sign, block, and Jordan identities;
- strongest route verdict and next executable representation.

No production numerics are licensed by this README. Maximum 0021 verdict
before nonlinear estimates is an exact full-contour monodromy theorem or a
quantified reduced filament-to-patch Floquet transfer. It cannot license
nonlinear orbital stability unless a separate nonlinear iteration closes,
nor full three-dimensional/swirl stability, a restoring frequency without a
semisimple physical mode and kinetic normalization, helicity, intrinsic or
quantum spin, electron/neutrino identification, parent completion, or a
global no-go.

## Activation handoff

Root should centrally register 0021 and place the schema receipts here.
Substantive source or derivation work begins only after
activation-schema.exit records zero.
