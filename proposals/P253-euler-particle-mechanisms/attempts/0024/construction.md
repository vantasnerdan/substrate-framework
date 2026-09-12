# P253/0024 — conditional essential-circle transfer and exact companion power law

## Result in one paragraph

The García–Hassainia–Hmidi normal form transfers exactly to the actual
full-label problem only through its diagonal principal return and its compact
order-minus-three, separated-cross, and finite-center terms. Its additional
finite-step residual \(E_{4,n}\) is controlled on the source-restricted tame
spaces, but the paper proves neither that its lift to the full IVP is compact
nor that those lifts converge to zero. Thus a compact-perturbation theorem for
the actual reduced monodromy is not yet established. If that one missing
compact-transfer condition is supplied, the source Cantor condition forces
the diagonal rotation ratio to be irrational and the essential spectrum is
the whole unit circle, ruling out the registered isolated high-mode Riesz
clusters. Unconditionally, an exact power formula for the actual
phase/translation/energy/impulse flag shows that a power-bounded reduced shape
map can still produce at most quadratic companion/phase growth unless its
resonant shears are modulated away. No route presently proves all-power
boundedness.

## 1. Physical normalization and fixed reduced object

Let \(\rho=r^2/2\) be the canonical radial coordinate and let
\(\rho_m>0\) be constant material density. The 0015 correction fixes

\[
 E_{\rm kin}^{\rm phys}=2\pi\rho_mE,\quad
 I_z^{\rm phys}=2\pi\rho_mP,\quad
 \Omega^{\rm phys}=2\pi\rho_m\Omega,\quad
 \mathscr A^{\rm phys}=2\pi\rho_m\mathscr A.                 \tag{1}
\]

Because the same positive factor multiplies the Hamiltonian and KKS form,

\[
 i_X\Omega^{\rm phys}=d\mathscr H_c^{\rm phys}
 \quad\Longleftrightarrow\quad
 i_X\Omega=d\mathscr H_c.                                   \tag{2}
\]

Thus normalization changes neither the Jacobi propagator nor its return map.
Take the P253/0021 candidate full-contour monodromy \(M\) on the two
independent individual-area tangents and its exact invariant flag

\[
 G=\operatorname{span}\{v_{\rm ph},v_z\}
 \subset K=\ker d\mathscr H_c\cap\ker dP\subset X^s.         \tag{3}
\]

Fix bounded phase and axial-centering functionals and identify

\[
 A:=M_{\rm red}:\mathcal N^s_{\rm red}=K/G
                    \longrightarrow\mathcal N^s_{\rm red}. \tag{4}
\]

All spectral statements below use the complexification of this real Sobolev
space. The physical perturbation class remains two-label and non-reversible;
finite center and constraint coordinates are removed only by (3).

## 2. What the source transformations actually transfer

Write \(L=|\log\epsilon|\), use source phase
\(\phi=\omega Lt\), and hence

\[
 \mathcal T=\frac{2\pi}{\omega L}.                           \tag{5}
\]

On the source higher-mode space, Proposition 8.3 gives

\[
 \Phi^{-1}\mathcal L_\perp\Phi
 =\epsilon^2L\omega\partial_\phi+D_\epsilon
  +\epsilon^2L^{1/2}K(\phi)+E_{4,n},                         \tag{6}
\]

where

\[
 D_\epsilon=c_1\partial_\theta+c_2\mathcal H
                 +\epsilon^2c_3\partial_\theta\Lambda_2,
 \qquad D_\epsilon e^{ij\theta}=i\mu_{j,2}e^{ij\theta},     \tag{7}
\]

\[
 \mu_{j,2}=c_1j+c_2\operatorname{sign}j+\epsilon^2c_3d_j,
 \qquad d_j=-\frac{j}{4(j^2-1)|j|},\quad |j|\ge2.            \tag{8}
\]

Here \(K=\Pi^\perp\partial_\theta R_{4,-3}\Pi^\perp\) belongs to
\(\mathrm{OPS}^{-3}\), and \(E_{4,n}\) is the source finite-step smoothing
residual.
The four source transformations are periodic and instantaneous in \(\phi\):
angular diffeomorphisms, bounded pseudodifferential exponentials, and an
exponential of \(\rho_1(\phi)\mathcal H+\rho_2(\phi)
\partial_\theta\Lambda_2\). Consequently they describe a legitimate
time-dependent IVP conjugation on their stated higher-mode domain. They do
not eliminate the last two terms in (6).

Three restrictions prevent (6) itself from being the full desired
similarity:

1. the source domain imposes reversibility and the half-period label relation;
2. \(\Pi^\perp\) removes the physical first/center modes before reduction;
3. the persistent \(K(\phi)\) and finite-step \(E_{4,n}\) are not proved
   skew-adjoint in a positive metric and are not exactly reduced away.

The singular principal calculation nevertheless transfers labelwise. Apply
the source transport straightening separately to labels 1 and 2, with the
second label's coefficients shifted by half a period. The two self blocks
then have the same one-period diagonal return. The independent-label cross
kernel is smooth because the actual contours remain separated, while the
self remainder in (6) has order \(-3\). First-mode and centering corrections
are finite rank. These pieces are compact after return to a fixed Sobolev
level.

The remaining source term requires separate notation. On the restricted
periodic space, Proposition 8.3 controls \(E_{4,n}\) at the base norm by a
decaying finite-step estimate involving two higher derivatives of the input
and a high-norm coefficient of the approximate state. Neither Proposition
8.3 nor the source existence theorem constructs a common full-label IVP
operator \(\mathcal E_n(t):H^s\to H^s\) that is compact, or a sequence of
full-label conjugations for which \(\mathcal E_n\to0\) in operator norm.
Consequently the honest lifted structural identity is

\[
 \partial_t y=\bigl(-\epsilon^{-2}D_\epsilon
              +\mathcal K_c(t)+\mathcal E_n(t)\bigr)y,       \tag{9}
\]

where the full-label lift of \(\mathcal K_c(t):H^s\to H^{s+\sigma}\),
\(\sigma>0\), is compact, while the existence and compactness of the displayed
\(\mathcal E_n(t)\) on the actual IVP is precisely unproved. Equation (9) is
therefore an interface decomposition, not a completed similarity.

Let

\[
 U_0(t)=\exp(-t\epsilon^{-2}D_\epsilon).
\]

Introduce the explicit compact-transfer hypothesis

\[
 \mathbf{CT}:\quad \mathcal E_n(t)\text{ has a full two-label IVP lift and
 is compact on }H^s\text{ (or vanishes in an operator-norm limit).}          \tag{10}
\]

Under \(\mathbf{CT}\), absorb \(\mathcal E_n\) into \(\mathcal K_c\).
Variation of constants gives

\[
 U(\mathcal T,0)-U_0(\mathcal T)
 =\int_0^{\mathcal T}
   U_0(\mathcal T-t)\mathcal K_c(t)U(t,0)\,dt.               \tag{11}
\]

Every integrand is compact on \(H^s\), norm-continuous in \(t\), and the
compact operators are norm closed. Restoring the bounded periodic
conjugation and the two finite-dimensional reductions therefore yields

\[
 \boxed{A=S\,(U_D+\mathcal C)S^{-1},
 \qquad \mathcal C\text{ compact on }\mathcal N^s_{\rm red},} \tag{12}
\]

conditional on \(\mathbf{CT}\), where \(S^{\pm1}\) are bounded and \(U_D\)
is two label copies of the diagonal unitary return map. The labelwise
principal and compact-cross calculation itself imposes no label relation;
what remains missing is the full-label control of \(E_{4,n}\).

The common physical factor (1) multiplies the symplectic form but has no
effect on compactness, similarity, or spectrum.

## 3. Exact essential spectrum

For \(j\ge2\), (5)-(8) give the diagonal multiplier

\[
 m_j=\exp\left[-\frac{2\pi i}{\epsilon^2L\omega}
       \left(c_1j+c_2+\epsilon^2c_3d_j\right)\right].        \tag{13}
\]

For an actual final source parameter \(\lambda\in\mathcal C_\infty\), fix
\((\ell,j)\) and take the Nash--Moser index \(m\to\infty\) in
\(\lambda\in O_m^1(g_m)\). Since \(g_m\to g_\infty\) and \(|j|\le N_m\)
eventually, the limiting first-transport condition is

\[
 |\epsilon^2L\omega\ell+c_1j|
 \ge \nu |j|^{-\tau},
 \qquad (\ell,j)\in\mathbb Z^2,\quad j\ne0.                 \tag{14}
\]

If \(c_1/(\epsilon^2L\omega)=p/q\) were rational, the choice
\((\ell,j)=(-mp,mq)\), with \(m\) large enough that \(|mq|\ge2\), would
make the left side of (14) zero. Hence

\[
 \alpha:=\frac{c_1}{\epsilon^2L\omega}\notin\mathbb Q.      \tag{15}
\]

The sequence \(\{e^{-2\pi i\alpha j}:j\ge2\}\) is dense in the unit
circle. Multiplication by the fixed phase from \(c_2\) and the vanishing
phase correction \(d_j\to0\) preserve density. A diagonal normal operator's
essential spectrum is the closure of its eigenvalue accumulation set,
including eigenvalues of infinite multiplicity. Thus

\[
 \sigma_{\rm ess}(U_D)=\mathbb S^1.                          \tag{16}
\]

Bounded similarity, finite-dimensional reduction, and compact perturbation
preserve Fredholm essential spectrum, so the conditional identity (12) gives

\[
 \boxed{\mathbf{CT}\Longrightarrow
        \sigma_{\rm ess}(M_{\rm red})=\mathbb S^1.}          \tag{17}
\]

Equation (17) neither places discrete eigenvalues on the circle nor excludes
nontrivial Jordan or nonnormal behavior embedded in the essential spectrum.

## 4. Route A: why partial source reducibility does not control powers

If the persistent terms in (6) are temporarily retained only through the
operator norm of their source-restricted realization
\(\mathcal R_n=\mathcal K_c+\mathcal E_n\), physical-time energy gives

\[
 \|U(\mathcal T,0)\|
 \le \exp\left(\int_0^{\mathcal T}\|\mathcal R_n(t)\|dt\right). \tag{18}
\]

The displayed \(\epsilon^2L^{1/2}K\) term becomes order \(L^{1/2}\) after
division by the physical-time coefficient \(\epsilon^2\), so over (5) its
gross accumulated norm is \(O(L^{-1/2})\). This yields, on exactly the
domain where all source transformations and tame losses close, a bound of
the form

\[
 \|U^n\|\le \kappa(S)
       \exp\bigl(CnL^{-1/2}+n\delta_{4,n}\bigr),             \tag{19}
\]

where \(\delta_{4,n}\) is the integrated finite-step residual in the chosen
tame norm. It provides a finite-window estimate when
\(n(L^{-1/2}+\delta_{4,n})=O(1)\), not a uniform all-power bound.

This limitation is logical, not just a weak estimate. On a symplectic Hilbert
space, append to a diagonal unitary with dense essential spectrum the
two-dimensional Hamiltonian block

\[
 H_\delta=\begin{pmatrix}e^\delta&0\\0&e^{-\delta}\end{pmatrix}.
                                                                    \tag{20}
\]

For every \(\delta>0\), \(H_\delta\) is symplectic, differs from the identity
by an arbitrarily small finite-rank perturbation as \(\delta\downarrow0\),
and preserves the diagonal model's same unit-circle essential spectrum, but
\(\|H_\delta^n\|=e^{n\delta}\). Therefore symplecticity, small compact
remainder, and even the conditional conclusion (17) do not imply power
boundedness.

**Route A verdict — blocked.** The source supplies a periodic bounded
conjugation only on its restricted sector and leaves a smoothing remainder.
The missing construction is either an exact full-label reducibility with
second-Melnikov divisors and zero final remainder, or a uniformly positive
periodic metric \(G(t)\) satisfying
\(\dot G+A_*^*G+GA_*=0\). Equation (20) proves that smallness of the
remainder cannot replace that construction. This verdict is not instability
of the actual patch pair.

## 5. Route B: direct isolated-cluster Riesz analysis

Route B preregistered uniformly separated high-mode contours and bounded
Riesz projections. The exact principal multipliers are dense on
\(\mathbb S^1\) by (15). Under the single compact-transfer condition
\(\mathbf{CT}\), (17) then rules out any family of disjoint contours that
isolates individual high modes or uniformly finite clusters at a positive
uniform distance from the remaining spectrum. A compact perturbation could
also create discrete reciprocal eigenvalue pairs off the circle unless an
additional Krein-definiteness theorem excludes them.

Without \(\mathbf{CT}\), however, the source finite-step residual is not
licensed as a compact perturbation of the actual full-label monodromy. The
dense diagonal model therefore cannot by itself determine the actual
essential spectrum or refute an actual-map Riesz representation.

**Route B verdict — blocked.** The missing construction is a full-label
operator-norm limit of the finite-step conjugations, or a direct proof that
their induced return-map residual is compact on the declared reduced Sobolev
space. Conditional on that construction, the isolated finite-cluster version
of Route B is obstructed by the dense unit circle; a grouped essential-
spectrum or global positive-metric method would remain distinct.

## 6. Route C: exact full-flag powers and the unresolved patch Jordan block

Retain the exact 0021 decomposition

\[
 M=\begin{pmatrix}
 A&0&B\\ C&I&D\\0&0&I
 \end{pmatrix},                                              \tag{21}
\]

where \(A=M_{\rm red}\), the middle block is phase/translation, and the last
block contains normalized energy/impulse companions. Direct induction gives

\[
 M^n=\begin{pmatrix}
 A^n&0&B_n\\ C_n&I&D_n\\0&0&I
 \end{pmatrix},                                              \tag{22}
\]

with

\[
\begin{aligned}
 B_n&=\sum_{k=0}^{n-1}A^kB,\\
 C_n&=\sum_{k=0}^{n-1}CA^k,\\
 D_n&=nD+\sum_{k=0}^{n-2}(n-1-k)CA^kB.                       \tag{23}
\end{aligned}
\]

Consequently, if the still-open reduced estimate
\(\sup_n\|A^n\|\le K_A\) is later proved, then the actual unmodulated full
flag obeys

\[
\begin{aligned}
 \|B_n\|&\le nK_A\|B\|,\\
 \|C_n\|&\le nK_A\|C\|,\\
 \|D_n\|&\le n\|D\|+\tfrac12n(n-1)K_A\|C\|\|B\|.         \tag{24}
\end{aligned}
\]

Thus the companion-to-shape and shape-to-phase couplings can produce
quadratic secular growth even when the reduced shape map is power bounded.
They must be modulated, not deleted. Coboundary components telescope: for
example, if \(B=(I-A)Y\), then

\[
 B_n=(I-A^n)Y,                                               \tag{25}
\]

so only the component of \(B\) not representable as an \((I-A)\)-coboundary
can drive linear accumulation. Because \(1\in\sigma_{\rm ess}(A)\), a bounded
inverse of \(I-A\) is unavailable under \(\mathbf{CT}\); (25) must be proved
on each actual coupling
vector, not inferred from a spectral gap.

The source filament twist gives a genuine size-two Jordan chain only for its
finite-dimensional center system. The finite-core patch family exists on a
Cantor set and supplies no derivative producing the corresponding generalized
Jacobi vector. Therefore no actual patch nilpotent block \(N_k\) has yet been
identified.

**Route C verdict — blocked.** Equations (22)-(25) establish the exact
full-flag power calculus and the maximum quadratic companion growth
conditional on reduced power boundedness. The registered actual polynomial
normal form still requires bounded resonant projections and either a directly
constructed finite-core generalized Jacobi vector or a proof that each
coupling is a coboundary. Filament twist alone cannot fill this dependency.

## 7. Failure-derived next route

The failures of isolated clusters and small-remainder iteration select a
different next representation: construct a global positive symmetrizer over
the whole essential circle, or complete the smoothing KAM reduction with
pairwise divisors

\[
 |\epsilon^2L\omega\ell+\mu_{j,2}-\mu_{k,2}|
 \ge \gamma\,\langle\ell\rangle^{-\tau_1}
               \langle j-k\rangle^{-\tau_2}.                \tag{26}
\]

This is a second-Melnikov condition and is not the source condition (14).
Its next executable tasks are:

1. derive transversality in \(\lambda\) for every divisor in (26), including
   opposite-sign modes where the \(c_2\) terms do not cancel;
2. prove a nonempty/asymptotically large sub-Cantor set inside the actual
   source set rather than assuming a chosen \(\lambda\) satisfies it;
3. solve the smoothing homological equations on both independent labels and
   the finite center/shape cross block with summable tame losses;
4. show the limiting conjugation is bounded, symplectic, periodic, and yields
   a positive constant diagonal metric.

Completion of those steps would give similarity to a unitary and hence the
target \(\sup_n\|A^n\|<\infty\). If a divisor vanishes, its resonant finite
block must instead be retained and classified by its actual Krein form; only
then can Route C assign a polynomial or hyperbolic verdict.

## 8. Strongest exact verdict and exclusions

The strongest new unconditional results are the limiting irrationality (15),
compactness of the order-minus-three, separated-cross, and finite-center
return contributions, and the exact actual full-flag identities (22)-(25).
Under the separately exposed compact-transfer condition \(\mathbf{CT}\), the
additional exact conclusion is (17). The current source and interface
estimates do not decide the actual reduced essential spectrum, off-circle
eigenvalues, embedded Jordan structure, or existence of a positive invariant
metric.

Nothing here licenses power boundedness, nonlinear orbital stability, an
all-time open Euler neighborhood, a global stability or instability theorem,
general three-dimensional or swirl perturbations, helicity, intrinsic or
quantum spin, electron or neutrino identification, parent completion, or a
global no-go.
