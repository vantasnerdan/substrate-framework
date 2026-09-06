# P253/0015 — actual periodic-pair variational and Floquet construction

## Result in one sentence

The García–Hassainia–Hmidi pair is an exact critical loop of the local
Euler/Kirillov–Kostant–Souriau (KKS) action on its labeled patch leaf, but it
cannot be a maximizer of the autonomous translating energy–impulse
functional: the complete first variation pairs with the nonzero
translating-frame velocity. The raw loop action is also strongly
indefinite. What remains is a reduced evolutionary Floquet estimate on the
full, non-reversible patch tangent space; the source's periodic
Nash–Moser inverse does not supply that estimate.

## 1. Exact source family and conventions

Write \(L=|\log\epsilon|\), use

\[
 \rho=\frac{r^2}{2},\qquad x=(\rho,z),\qquad dx=d\rho\,dz,
\]

and let

\[
 \{a,b\}=a_\rho b_z-a_zb_\rho,
 \qquad \nabla^\perp a=(-a_z,a_\rho).
\]

The source equation is

\[
 \partial_tq+\{\Psi[q],q\}=0,
 \qquad
 \Psi[q](x)=\frac1{\sqrt2\,\pi}\int_{\Pi}G(x,y)q(y)\,dy.       \tag{1}
\]

For fixed \(0<a<b\), \(\kappa>0\), sufficiently small \(\epsilon\), and
\(\lambda\in\mathcal C_\epsilon\subset(a,b)\), source Theorem 1.1 gives

\[
 q^{\rm lab}_*(t)=q_{1,*}^{\rm lab}(t)+q_{2,*}^{\rm lab}(t),
 \qquad
 q_{j,*}^{\rm lab}(t)=\epsilon^{-2}{\bf1}_{D_j^{\rm lab}(t)}. \tag{2}
\]

The slow time is \(\tau=Lt\). Corollary 3.1 splits the common axial
center into \(U_\epsilon\tau+W_\epsilon(\tau)\), with
\(W_\epsilon\) periodic of period
\(T=T(\epsilon,\lambda,\kappa)\). Thus the physical translating speed and
period are

\[
 c_\epsilon=L U_\epsilon,
 \qquad \mathcal T_\epsilon=\frac{T}{L}.                       \tag{3}
\]

Define the uniformly translating-frame labels by

\[
 q_{j,*}(t,\rho,z)
 =q_{j,*}^{\rm lab}(t,\rho,z+c_\epsilon t).                    \tag{4}
\]

The source periodicity and its half-period label relation imply
\(q_{j,*}(t+\mathcal T_\epsilon)=q_{j,*}(t)\). The pair is labeled:
the half-period exchange does not identify the two factors of the leaf.

The source shape variable obeys \(w_j^2=1+2\epsilon f_j\), with
\(\int_0^{2\pi}f_j(\tau,\theta)d\theta=0\). The ellipse map

\[
 Z_j=(2p_{j,1})^{1/4}\cos\theta
      +i(2p_{j,1})^{-1/4}\sin\theta
\]

has determinant one. Hence the normalized cross-section area is exactly

\[
 \frac12\int_0^{2\pi}w_j^2\det(Z_j,\partial_\theta Z_j)d\theta
 =\pi,
\]

the physical \((\rho,z)\)-area is \(\pi\epsilon^2\), and each label has
circulation

\[
 \int_\Pi q_{j,*}\,dx=\pi.                                    \tag{5}
\]

For small \(\epsilon\), the source center separation is
\(O(L^{-1/2})\), whereas each core is \(O(\epsilon)\); the two material
domains are therefore disjoint. Common area-preserving transport retains
their amplitudes, areas, simple connectivity, disjointness, and labels.

### Nonstationarity is an exact source fact

Proposition 3.4 gives

\[
 p_{1,1}(\tau+T/2)=p_{2,1}(\tau),
\]

and at \(\tau=0\)

\[
 p_{1,1}(0)-p_{2,1}(0)
 =r_\epsilon(2\kappa)^{1/4}\lambda\ne0.                        \tag{6}
\]

Consequently the first labeled patch changes radial location after half a
period even after the common axial drift is removed. The translating-frame
pair is a nonconstant periodic orbit, not a relative equilibrium.

## 2. Actual invariant labeled leaf and conserved functional

Let \(D_{j,0}\) be the translating-frame source domains at \(t=0\). The
smooth material coadjoint leaf is

\[
 \mathcal O_j=
 \left\{\epsilon^{-2}{\bf1}_{\eta_j(D_{j,0})}:
 \eta_j\in\operatorname{Diff}_{c,\rm area}(\Pi),\
 \eta_j(\Pi)=\Pi\right\},
 \qquad
 \mathcal O_{\rm pair}=\mathcal O_1\times\mathcal O_2.         \tag{7}
\]

Only a sufficiently small open subset with disjoint supports is needed for
the independent variations below. Its equimeasurable weak closure can be
used for compactness questions, but it is not the smooth KKS leaf and is not
used to manufacture a first variation.

The common Euler flow, followed by the common axial translation (4), is an
area-preserving diffeomorphism. Therefore

\[
 (q_1(t),q_2(t))\in\mathcal O_{\rm pair}
\]

for every source time. Independent elements of the product leaf are not
claimed to evolve by independent velocities: both labels obey the same
total stream function.

With \(q=q_1+q_2\), define

\[
 E(q)=\frac12\int_\Pi q\Psi[q]\,dx
 =\frac1{2\sqrt2\,\pi}\iint_{\Pi^2}G(x,y)q(x)q(y)\,dxdy,
 \qquad
 P(q)=\int_\Pi\rho q(x)\,dx.                                  \tag{8}
\]

The last expression is Cao's
\(P=\frac12\int r^2\zeta\,r\,dr\,dz\) after \(\rho=r^2/2\).
Symmetry of \(G\) gives \(\delta E/\delta q=\Psi[q]\). Axial
translation invariance gives conservation of \(P\), and Euler transport
conserves \(E\). Thus the actual translating Hamiltonian

\[
 \boxed{\mathscr H_c(q_1,q_2)=E(q_1+q_2)-c_\epsilon P(q_1+q_2)} \tag{9}
\]

is conserved and contains both self energies and the full cross energy.
No individual impulse or self/cross energy is asserted to be conserved.

Writing

\[
 h[q]=\frac{\delta\mathscr H_c}{\delta q_j}
     =\Psi[q]-c_\epsilon\rho,                                  \tag{10}
\]

the two translating-frame label equations are exactly

\[
 \boxed{\partial_tq_j=\{q_j,h[q]\}},\qquad j=1,2.              \tag{11}
\]

Indeed, (11) is equivalent to
\(\partial_tq_j+\{\Psi-c_\epsilon\rho,q_j\}=0\), the result of
differentiating (4).

## 3. Route A: complete first variation and exact obstruction

Let \(\chi_j\in C_c^\infty(\Pi)\) generate an area-preserving flow with
velocity \(\nabla^\perp\chi_j\). The corresponding labelwise isovortical
tangent is

\[
 \delta_{\chi_j}q_j=\{q_j,\chi_j\}.                            \tag{12}
\]

Integration by parts in the Poisson bracket, with no boundary contribution,
gives the complete first variation

\[
\begin{aligned}
 D\mathscr H_c(q)[\delta_\chi q]
 &=\sum_{j=1}^2\int_\Pi h[q]\{q_j,\chi_j\}\,dx\\
 &=\sum_{j=1}^2\int_\Pi q_j\{\chi_j,h[q]\}\,dx\\
 &=-\sum_{j=1}^2\left\langle\chi_j,\{q_j,h[q]\}\right\rangle.
                                                                    \tag{13}
\end{aligned}
\]

At a source phase, equation (11) turns this into

\[
 \boxed{
 D\mathscr H_c(q_*(t))[\delta_\chi q_*]
 =-\sum_{j=1}^2\langle\chi_j,\partial_tq_{j,*}(t)\rangle.}     \tag{14}
\]

Equation (6), together with uniqueness of the autonomous contour evolution,
shows that the full state velocity never vanishes: if it vanished at one
phase, that phase would be an equilibrium and the entire orbit would be
constant. Thus at every phase at least one distribution
\(\partial_tq_{j,*}\) is nonzero. There is a compactly supported smooth
\(\chi_j\), localized near that patch, for which the pairing in (14) is
nonzero. In boundary language,

\[
 \langle\chi_j,\partial_tq_{j,*}\rangle
 =\epsilon^{-2}\int_{\partial D_{j,*}}\chi_j V_{n,j}\,ds,       \tag{15}
\]

so the named noncritical direction is any area-preserving generator whose
boundary trace has nonzero pairing with the actual nonzero normal velocity.

A differentiable local maximum on the smooth leaf must annihilate every
leaf tangent. Equation (14) supplies an admissible tangent on which it does
not. Therefore neither a single phase nor the entire phase/translation
orbit can be a local or global maximizing set of the autonomous functional
(9). Conservation merely gives the one tangent identity

\[
 D\mathscr H_c(q_*)[\partial_tq_*]=0;
\]

it does not give \(D\mathscr H_c(q_*)=0\).

**Route A verdict — refuted.** The exact mechanism is noncriticality on
the actual labeled rearrangement leaf, witnessed by (14). Cao's steady
one-ring maximizer theorem cannot transfer across this obstruction.

## 4. Route B1: the correct local KKS loop action

For tangents (12), define the product-leaf KKS form

\[
 \Omega_q(\delta_\chi q,\delta_\psi q)
 =-\sum_{j=1}^2\int_\Pi q_j\{\chi_j,\psi_j\}\,dx.               \tag{16}
\]

The quotient by generators that stabilize \(q_j\) makes (16) well defined
and nondegenerate. Equations (13) and (16) give

\[
 \Omega_q(X_{\mathscr H_c}(q),\delta_\chi q)
 =D\mathscr H_c(q)[\delta_\chi q],
 \qquad
 X_{\mathscr H_c}(q)_j=\{q_j,h[q]\}.                           \tag{17}
\]

A global primitive of \(\Omega\) need not exist. For a loop \(q\) in a
small coadjoint chart and a homotopy
\(Q:[0,1]\times S^1_{\mathcal T}\to\mathcal O_{\rm pair}\) from a fixed
reference loop to \(q\), use a local potential satisfying
\(-d\Theta=\Omega\) and define the relative/local action

\[
\boxed{
 \mathscr A_{c,\mathcal T}[q;Q]
 =-\int_{[0,1]\times S^1_{\mathcal T}}Q^*\Omega
  -\int_0^{\mathcal T}\mathscr H_c(q(t))\,dt.}                 \tag{18}
\]

Changing the filling can change the value by a symplectic period, but not
its local first variation. For every periodic tangent \(\eta\), orientation
of the filling can be chosen so that

\[
 \delta\mathscr A_{c,\mathcal T}[q](\eta)
 =\int_0^{\mathcal T}
   \left[-\Omega_q(\eta,\dot q)-D\mathscr H_c(q)[\eta]\right]dt.
                                                                    \tag{19}
\]

Here (17) is precisely \(i_{X_{\mathscr H_c}}\Omega=d\mathscr H_c\);
therefore (19) vanishes for all periodic \(\eta\) exactly when
\(\dot q=X_{\mathscr H_c}(q)\). Equations (3)-(4) make the published pair a
closed labeled loop and (11) is precisely this Hamilton equation.

**Route B1 verdict — established.** The actual published pair, rather than
a frozen center-history surrogate, is a critical loop of the local KKS
action (18), modulo the usual filling-period ambiguity.

## 5. Route B2: why the raw action is not a restoring functional

In a local symplectic trivialization along the loop, the action Hessian has
the standard first-order form

\[
 Q(\eta)=\int_0^{\mathcal T}
 \left[-\Omega_0(\eta,\dot\eta)-\langle S(t)\eta,\eta\rangle\right]dt,
                                                                    \tag{20}
\]

where \(S(t)\) contains the spatial second variation and connection terms.
Choose a fixed smooth spatial tangent two-plane with
\(\Omega_0(e,f)=1\), and set, for
\(n=2\pi k/\mathcal T\),

\[
 \eta_\pm(t)=e\cos(nt)\pm f\sin(nt).
\]

Then

\[
 \int_0^{\mathcal T}\Omega_0(\eta_\pm,\dot\eta_\pm)dt
 =\pm n\mathcal T.                                               \tag{21}
\]

Thus the kinetic contributions to (20) are \(\mp n\mathcal T\).
The \(S(t)\) contribution stays \(O(1)\) as \(k\to\infty\), while these
have both signs and grow like \(k\). The plane may be chosen orthogonal to
any fixed finite collection of phase, translation, area, energy, and impulse
conditions. Therefore the raw Hamilton action is strongly indefinite even
after those finite-dimensional neutral constraints are removed. It cannot
be a positive Lyapunov deficit or a direct action maximum.

**Route B2 verdict — refuted.** A sign-definite Hessian/coercive maximum of
the raw loop action is excluded by the high-time-frequency symplectic term.
This does not refute a reduced Floquet estimate or a different positive
modulated norm.

## 6. Route B3: exact Jacobi equation and neutral periodic modes

Let \(\eta=(\eta_1,\eta_2)\) be a labelwise leaf tangent and put

\[
 \Phi_\eta(x)=\frac1{\sqrt2\,\pi}
               \int_\Pi G(x,y)(\eta_1+\eta_2)(y)\,dy.
\]

Linearizing (11) gives the full translating-frame Jacobi equation

\[
 \boxed{
 \partial_t\eta_j
 =\{\eta_j,\Psi_*-c_\epsilon\rho\}
  +\{q_{j,*},\Phi_\eta\},\qquad j=1,2.}                         \tag{22}
\]

It retains the cross interaction through \(\Phi_\eta\). Its tangent space
has \(\int\eta_jdx=0\), and the patch version consists of boundary-supported
normal variations generated by (12); tangential contour changes are gauge.

Two exact periodic Jacobi fields are

\[
 \eta_{\rm ph}=\partial_tq_*,
 \qquad
 \eta_{z}=-\partial_zq_*.                                      \tag{23}
\]

The first follows by time differentiation of the autonomous equation; the
second follows from axial-translation equivariance. They are independent:
the phase evolution changes the relative radial/axial configuration in
(6), whereas a common axial translation changes neither radial center.

If (22) is realized as a differentiable evolution family
\(U(t,s)\) on a source-admissible full contour tangent space, its period map
is

\[
 \mathcal M=U(\mathcal T_\epsilon,0).                           \tag{24}
\]

For smooth Jacobi solutions the Hamiltonian calculation gives the exact
identity

\[
 \Omega_{q_*(t)}(U(t,0)\eta,U(t,0)\xi)
 =\Omega_{q_*(0)}(\eta,\xi),                                   \tag{25}
\]

so any bounded realization of \(\mathcal M\) is symplectic. Equations
(23) give, conditionally on that realization,

\[
 \mathcal M\eta_{\rm ph}(0)=\eta_{\rm ph}(0),
 \qquad
 \mathcal M\eta_z(0)=\eta_z(0).                               \tag{26}
\]

Thus the exact full problem has at least the phase and translation periodic
Jacobi modes. Equations (22)-(26) are structural identities; they do not
place the remaining multipliers on the unit circle and do not exclude
Jordan chains.

The source's operator \(\partial_gF(\epsilon,g_\infty)\) is the periodic
boundary-value version of a reduced piece of (22), but its domain has already
imposed all of the following:

1. zero spatial average (fixed patch area);
2. the half-period label relation \(f_2(t)=f_1(t+T/2)\);
3. reversibility/evenness in the construction space;
4. speed modulation eliminating the first sine mode;
5. separate treatment or projection of the first cosine mode.

Sections 8-9 diagonalize the normal part only up to smoothing and iteration
errors and construct approximate tame right inverses on nonresonant Cantor
sets. This proves the solvability needed for the existence scheme. It is
not a bounded full initial-value propagator, does not cover arbitrary nearby
patch perturbations, and does not prove Floquet spectral or nonlinear
stability.

**Route B3 verdict — established.** The exact full Jacobi equation, its
KKS conservation law, and its two symmetry-generated periodic modes are
constructed. The existence of a bounded full-source-space period map and a
restoring estimate are deliberately not included in this verdict.

## 7. Exact source-filament monodromy: the missing amplitude companion

The reduced filament system is one degree of freedom:

\[
 \dot x_1=-\partial_{x_2}H_\epsilon,
 \qquad
 \dot x_2=\partial_{x_1}H_\epsilon.
\]

Let \(x(t;\lambda)\) be its periodic orbit and
\(T(\lambda)=T(\epsilon,\lambda,\kappa)\). With

\[
 v_{\rm ph}=\dot x(0;\lambda),
 \qquad v_\lambda=\partial_\lambda x(0;\lambda),
\]

differentiate the exact periodicity identity
\(x(T(\lambda);\lambda)=x(0;\lambda)\). If \(M_{\rm fil}\) is the
linearized return map, then

\[
 \boxed{(M_{\rm fil}-I)v_\lambda
       =-T'(\lambda)v_{\rm ph},\qquad
       (M_{\rm fil}-I)v_{\rm ph}=0.}                            \tag{27}
\]

Source Proposition 3.2 proves
\(\inf_{\lambda\in[a,b]}T_0'(\lambda,\kappa)>0\), and Proposition 3.3 gives

\[
 T'(\lambda)=T_0'(\lambda,\kappa)+O(L^{-1/2}).                 \tag{28}
\]

Hence \(T'(\lambda)>0\) for sufficiently small \(\epsilon\). At the
source initial point, \(v_\lambda=(1,0)\), while
\(v_{\rm ph}=(0,\dot x_2(0))\ne0\), so (27) is a genuine Jordan chain at
the unit multiplier. The source filament orbit is parabolic in the
unconstrained amplitude direction, not a strict restoring oscillator.

**Route B4 verdict — established for the source filament subsystem.** The
exact return map has a nontrivial unit-multiplier Jordan chain generated by
period twist. This identifies an amplitude/energy companion that a
phase-only slice misses. It is not promoted to the full Euler-patch
monodromy: the patch solutions exist only for a Cantor set of \(\lambda\),
and the source does not differentiate a full Euler flow family through an
open \(\lambda\)-interval.

## 8. Remaining restoring construction

The remaining claim is sharply separated from the established action and
Jacobi identities.

**Route B5 verdict — blocked with the missing construction named.** No
source theorem or exact calculation above bounds the full patch period map
on the complement of phase, translation, and amplitude/energy shear.
Therefore neither spectral Floquet stability nor nonlinear orbital
stability of the fixed orbit \(\Gamma\) is established.

The next executable analytic construction is:

1. Reinstate two independent source \(H^s\) boundary normal graphs, retaining
   zero area but dropping the existence proof's reversibility and
   half-period label restrictions. Fix only tangential reparametrization.
2. Prove tame initial-value well-posedness for the full linearized contour
   equation corresponding to (22) through one source period. This earns a
   bounded \(U(t,s)\) and the actual \(\mathcal M\) in (24), rather than an
   inverse of a periodic invariance equation.
3. Perform Marsden-Weinstein/Poincare reduction at fixed total impulse and
   account for the four directions suggested by (23) and (27): phase,
   translation, and their energy/impulse companions. A modulation parameter
   for orbit amplitude is required; simply deleting the phase eigenvector
   leaves the filament Jordan shear.
4. Adapt the source's periodic symplectic conjugations to the initial-value
   generator on every \(|j|\ge2\) shape mode and control the smoothing
   remainder uniformly for arbitrarily many periods. Compute Krein signs or
   an equivalent conserved positive quadratic norm on the reduced space.
5. Only after that reduced linear bound, estimate the nonlinear contour
   remainder in the same \(H^s\) class and close an all-time modulation
   bootstrap for distance to \(\Gamma\).

This route retains the actual cross energy, common velocity, label exchange,
and all patch-shape modes. It introduces no prescribed target history and
uses no no-swirl result to license spin, helicity, an electron, or a
neutrino mechanism.

## 9. Scope ledger

| Claim | Exact status in 0015 | Reason |
|---|---|---|
| Actual pair leaf and \(E-cP\) conservation | established | Equations (7)-(11) |
| Autonomous maximum of \(E-cP\) | refuted | Nonzero admissible first variation (14) |
| KKS action-critical periodic patch loop | established | Equations (16)-(19) |
| Coercivity of raw loop action | refuted | Opposite high-frequency signs (20)-(21) |
| Full Jacobi equation and symmetry modes | established | Equations (22)-(26) |
| Filament amplitude multiplier | established, filament only | Nontrivial Jordan chain (27)-(28) |
| Full patch Floquet restoring estimate | blocked | Full propagator/reduced-space bound absent |
| Nonlinear orbital stability of \(\Gamma\) | open | Requires the five-step construction above |

This is a route-level checkpoint. It is not parent completion and is not an
obligation-level or global no-go.
