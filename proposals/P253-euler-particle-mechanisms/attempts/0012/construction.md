# P253/0012 — robust same-leaf pair and source-localized interval

## Scope and source table

This attempt keeps the axisymmetric no-swirl, two-label Euler system from
0006, but replaces its two exact initial maximizers by a genuine
support-controlled neighborhood on the two individual rearrangement leaves.
It also tests the short-time full-support theorem of Buttà–Marchioro against
an explicitly Euler-scaled Cao polynomial carrier.

| Input | Exact scope used here | Exclusion retained |
|---|---|---|
| Cao et al., Theorems 2.9–2.10 | Compactness of maximizing sequences modulo axial translation and qualitative orbital stability on one rearrangement leaf | No quantitative Hessian gap follows |
| Cao et al., Proposition 1.4 and Theorems 4.2–4.3 | Smooth polynomial-profile traveling ring and one-orbit maximizer set | Fixed carrier is axisymmetric without swirl |
| Cao et al., equation (2.6) | \(G_1(x,y)\le C_G(rr')^2/|x-y|^3\) | Gives an upper interaction bound, not a collision barrier |
| Buttà–Marchioro, Theorem 1.1 | Arbitrary definite-sign concentrated rings at initially disjoint fixed centers; complete support stays in disjoint radius-\(R\) disks for \(0\le t\le T_R\), where \(T_R>0\) is independent of the concentration parameter | \(T_R\) is positive but not all-time |
| Buttà–Marchioro, Assumption 2.1 and Theorem 2.3 | Intensity \(a/|\log\rho|\), \(L^\infty\le M/(\rho^2|\log\rho|)\), and an external field with divergence-free three-dimensional lift and \(O(|\log\rho|^{-1})\) size/Lipschitz constant | The source is a singular family, not the fixed-circulation family of 0006 |

All exact rearrangement, energy, and impulse objects below use Cao's relative
vorticity \(\zeta=\omega_\theta/r\), half-plane measure
\(d\nu=r\,dr\,dz\), and

\[
 E(\zeta)=\frac12\int\zeta G_1\zeta\,d\nu,\qquad
 P(\zeta)=\frac12\int r^2\zeta\,d\nu.
\]

## Route A — a genuine nearby same-leaf neighborhood

### A.1 Data class and exact invariants

Fix one smooth Cao carrier \(\zeta_*\), its traveling speed \(c>0\), and

\[
 \mathcal F(\zeta)=E(\zeta)-cP(\zeta),\qquad
 S=\sup_{\zeta\in\mathcal R(\zeta_*)^w}\mathcal F(\zeta).
\]

Its maximizing orbit is

\[
 \Sigma=\{\zeta_*(\cdot+a e_z):a\in\mathbb R\}
 \subset\mathcal R(\zeta_*).
\]

Let two nonnegative initial labels satisfy

\[
 \xi_j^0\in\mathcal R(\zeta_*),\qquad
 \operatorname{supp}\xi_j^0\subset
 [0,R_{\rm in}]\times I_j,\qquad j=1,2,
\]

where the axial intervals are ordered and have gap

\[
 g_0=\inf I_2-\sup I_1>0.
\]

Both labels have the carrier circulation

\[
 \int\xi_j^0\,d\nu=\kappa,
\qquad
 \|\xi_j^0\|_\infty=\Lambda.
\]

Evolve them by the common autonomous velocity

\[
 \partial_t\xi_j+
 (K[\xi_1]+K[\xi_2])\cdot\nabla\xi_j=0.
\]

The common volume-preserving flow retains
\(\xi_j(t)\in\mathcal R(\zeta_*)\), while injectivity retains the two material
labels. Define the actual conserved total impulse and actual initial cross
energy

\[
 P_T=P(\xi_1^0)+P(\xi_2^0),\qquad
 E_{12}^0=E_{12}(\xi_1^0,\xi_2^0).
\]

Individual impulses exchange:

\[
 \dot P_1=-\iint\xi_1(x)\partial_{z_x}G_1(x,y)\xi_2(y)\,d\nu_xd\nu_y,
 \qquad
 \dot P_2=-\dot P_1,
\]

so \(P_1(t)+P_2(t)=P_T\), without imposing
\(P_j(t)=P_j(0)\).

### A.2 Nonzero initial deficits and their exact evolution

Set

\[
 \Delta_j(t)=S-\mathcal F(\xi_j(t)),\qquad
 \delta_0=\Delta_1(0)+\Delta_2(0).
\]

Because each label stays in the exact rearrangement class,

\[
 \Delta_j(t)\ge0.
\]

Total Euler energy and total impulse are conserved. Expanding total energy as

\[
 E(\xi_1+\xi_2)=E(\xi_1)+E(\xi_2)+E_{12}
\]

and using the same coefficient \(c\) in both deficits gives

\[
 \boxed{
 \Delta_1(t)+\Delta_2(t)
 =\delta_0+E_{12}(t)-E_{12}^0.}
\]

This is the 0006 identity with its actual nonzero initial deficit restored.
Neither \(P_T=2P(\zeta_*)\) nor an exact-maximizer initial condition is used.
In particular,

\[
 E_{12}(t)\ge E_{12}^0-\delta_0
\]

is the exact lower consequence of nonnegative deficits.

### A.3 Cross-energy bound with the actual total impulse

If the two material supports have axial gap \(\ell>0\), the source Green bound
gives

\[
 \begin{aligned}
 E_{12}(t)
 &\le\frac{C_G}{\ell^3}
 \left(\int r^2\xi_1\,d\nu\right)
 \left(\int r^2\xi_2\,d\nu\right)\\
 &=\frac{4C_G}{\ell^3}P_1(t)P_2(t)
 \le\boxed{\frac{C_GP_T^2}{\ell^3}}.
 \end{aligned}
\]

The last inequality is
\((P_1+P_2)^2/4-P_1P_2=(P_1-P_2)^2/4\ge0\).

### A.4 Autonomous support interval

Set \(R=2R_{\rm in}\). While both supports remain in \(r\le R\),
disjoint label transport and equimeasurability give

\[
 \|\boldsymbol\omega\|_\infty
 \le R\Lambda,\qquad
 \boldsymbol\omega=r(\xi_1+\xi_2)e_\theta.
\]

The total three-dimensional vorticity \(L^1\) norm obeys

\[
 \begin{aligned}
 \|\boldsymbol\omega\|_1
 &=2\pi\int r(\xi_1+\xi_2)\,d\nu\\
 &\le2\pi
 \sqrt{\left(\int(\xi_1+\xi_2)d\nu\right)
       \left(\int r^2(\xi_1+\xi_2)d\nu\right)}
 =4\pi\sqrt{\kappa P_T}.
 \end{aligned}
\]

The optimized three-dimensional Biot–Savart split from 0006 therefore gives

\[
 \|u\|_\infty\le U_A,\qquad
 U_A=C_{BS}(4\pi\sqrt{\kappa P_T})^{1/3}
                 (2R_{\rm in}\Lambda)^{2/3},
\qquad
 C_{BS}=\frac{3}{2^{4/3}\pi^{1/3}}.
\]

Define

\[
 \boxed{
 T_A=\min\left\{\frac{R_{\rm in}}{U_A},
                 \frac{g_0}{4U_A}\right\}>0.}
\]

A bootstrap from particle displacement gives, for \(0\le t\le T_A\),

\[
 \operatorname{supp}\xi_j(t)\subset\{r\le2R_{\rm in}\},
\qquad
 \operatorname{gap}_z(\operatorname{supp}\xi_1(t),
                       \operatorname{supp}\xi_2(t))
 \ge\frac{g_0}{2}.
\]

Consequently,

\[
 \boxed{
 0\le\Delta_j(t)\le
 D_A:=\delta_0+\frac{8C_GP_T^2}{g_0^3}-E_{12}^0,
 \quad 0\le t\le T_A.}
\]

The retained \(-E_{12}^0\) is important: replacing it by zero would discard
known initial interaction energy and shrink the certified neighborhood.

### A.5 Qualitative orbital conclusion and a nonempty open neighborhood

For \(R\ge R_{\rm in}\), define the Cao orbital metric

\[
 d_R(\zeta,\Sigma)=\inf_{a\in\mathbb R}\left[
 \|\zeta-\zeta_*(\cdot+a e_z)\|_1+
 \|\zeta-\zeta_*(\cdot+a e_z)\|_2+
 P(|\zeta-\zeta_*(\cdot+a e_z)|)\right]
\]

and the source compactness modulus

\[
 m_R(\eta)=\inf\left\{
 S-\mathcal F(\zeta):
 \begin{array}{l}
 \zeta\in\mathcal R(\zeta_*),\
 \operatorname{supp}\zeta\subset\{r\le R\},\\
 d_R(\zeta,\Sigma)\ge\eta
 \end{array}\right\}.
\]

Cao compactness gives \(m_R(\eta)>0\) for every \(\eta>0\), without assigning
a rate. Therefore

\[
 \boxed{D_A<m_{2R_{\rm in}}(\eta)}
\]

implies

\[
 \boxed{
 d_{2R_{\rm in}}(\xi_j(t),\Sigma)<\eta,\qquad
 j=1,2,\quad 0\le t\le T_A.}
\]

This condition defines a genuine relative neighborhood, not a single
trajectory. To see nonemptiness, first choose \(\bar R\) above the carrier's
radial support and axially translate two exact maximizers far enough apart.
The Green bound makes both \(8C_G(2P(\zeta_*))^2/g_0^3\) and
\(E_{12}^0\) tend to zero as the separation tends to infinity, while
\(m_{2\bar R}(\eta)>0\) is fixed. Then choose strict bounds

\[
 \delta_0<\bar\delta,\quad
 P_T<\bar P,\quad
 E_{12}^0>\underline e,\quad
 g_0>\underline g,\quad
 R_{\rm in}<\bar R
\]

such that

\[
 \boxed{
 \bar\delta+\frac{8C_G\bar P^2}{\underline g^3}
 -\underline e<m_{2\bar R}(\eta).}
\]

Every same-leaf pair satisfying those strict support and invariant bounds has
the uniform velocity and time bounds

\[
 \bar U=C_{BS}(4\pi\sqrt{\kappa\bar P})^{1/3}
 (2\bar R\Lambda)^{2/3},\qquad
 \bar T=\min\left\{\frac{\bar R}{\bar U},
                   \frac{\underline g}{4\bar U}\right\}.
\]

The class is nonempty with nonzero deficits. Let \(\Phi_{j,s}\) be
weighted-volume-preserving diffeomorphisms generated by smooth fields
satisfying

\[
 \partial_r(rv_r)+\partial_z(rv_z)=0,
\]

supported inside the two initial support tubes, and not generating axial
translations. Then

\[
 \xi_{j,s}^0=\zeta_*(\cdot-a_je_z)\circ\Phi_{j,s}^{-1}
\in\mathcal R(\zeta_*).
\]

For small nonzero \(s\), their supports retain the strict bounds,
\(\Delta_j(0)\to0\), and uniqueness of the maximizing orbit implies
\(\Delta_j(0)>0\) whenever \(\xi_{j,s}^0\notin\Sigma\). Continuity of
energy, impulse, and separated-support cross energy retains the displayed
strict neighborhood inequality.

**Route A verdict — established as stated.** There is a nonempty
support-controlled relative-open neighborhood of genuinely deformed,
nonzero-deficit data on the two exact rearrangement leaves. Every member has
full labeled-vorticity orbital control on an autonomous interval derived from
its actual \(P_T\), \(E_{12}^0\), support gap, and deficit sum.

## Route B1 — exact Euler scaling into the Buttà–Marchioro class

### B.1 Scale covariance of the carrier and variational problem

Let \(u_\delta\) be a Cao traveling ring of thickness parameter \(\delta\),
circulation \(\kappa\), speed
\(c_\delta=W\log(1/\delta)\), relative vorticity \(\zeta_\delta\), and
azimuthal vorticity \(\omega_\delta=r\zeta_\delta\). Choose an axial translate
and a point \((0,r_\delta)\) in its support. The Cao diameter and
\(L^\infty\) bounds give constants \(C_s,C_\omega\), independent of small
\(\delta\), such that

\[
 \operatorname{supp}\omega_\delta
 \subset\Sigma((0,r_\delta)\mid C_s\delta),\qquad
 \|\omega_\delta\|_\infty\le C_\omega\delta^{-2},
\qquad r_\delta\to r_*=\frac{\kappa}{4\pi W}.
\]

For \(A,B>0\), Euler scaling is

\[
 u^{A,B}(x,t)=A\,u_\delta(Bx,ABt),
\]

with

\[
 \omega^{A,B}(x,t)=AB\,\omega_\delta(Bx,ABt),\qquad
 \zeta^{A,B}(r,z,t)=AB^2\zeta_\delta(Br,Bz,ABt).
\]

It transforms circulation, energy, impulse, and speed by

\[
 \kappa^{A,B}=\frac AB\kappa,\qquad
 E^{A,B}=\frac{A^2}{B^3}E,\qquad
 P^{A,B}=\frac A{B^3}P,\qquad
 c^{A,B}=Ac_\delta.
\]

Thus

\[
 E^{A,B}-c^{A,B}P^{A,B}
 =\frac{A^2}{B^3}(E-c_\delta P).
\]

The scaling bijects the corresponding rearrangement classes, so the scaled
Cao ring remains a maximizer, its maximizer set remains one axial-translation
orbit, and Cao's qualitative compactness theorem remains applicable. No
Hessian statement is introduced.

### B.2 Match every concentrated-ring hypothesis

Fix a desired physical radius \(r_0>0\) and source intensity \(a>0\). Set

\[
 B_\delta=\frac{r_\delta}{r_0},\qquad
 \rho_\delta=\frac{2C_s\delta}{B_\delta},\qquad
 \frac{A_\delta}{B_\delta}
 =\frac{a}{\kappa|\log\rho_\delta|}.
\]

Then \(\rho_\delta\to0\), and the scaled carrier
\(\bar\zeta_{\rho_\delta}=\zeta_\delta^{A_\delta,B_\delta}\) obeys

\[
 \operatorname{supp}\bar\omega_{\rho_\delta}
 \subset\Sigma((0,r_0)\mid\rho_\delta/2),
\qquad
 \int\bar\omega_{\rho_\delta}\,dr\,dz
 =\frac{a}{|\log\rho_\delta|}.
\]

Its pointwise bound is

\[
 \begin{aligned}
 \|\bar\omega_{\rho_\delta}\|_\infty
 &\le A_\delta B_\delta C_\omega\delta^{-2}\\
 &=\frac{aC_\omega B_\delta^2}
         {\kappa|\log\rho_\delta|\,\delta^2}
 =\frac{4aC_\omega C_s^2}
         {\kappa\rho_\delta^2|\log\rho_\delta|}.
 \end{aligned}
\]

Hence it satisfies the source \(M/(\rho^2|\log\rho|)\) hypothesis with one
\(\delta\)-independent \(M\). Its speed has the required limit:

\[
 A_\delta c_\delta
 =B_\delta\frac{aW}{\kappa}
   \frac{\log(1/\delta)}{|\log\rho_\delta|}
 \longrightarrow
 \frac{r_*}{r_0}\frac{aW}{\kappa}
 =\frac{a}{4\pi r_0}.
\]

This is an exact named Euler-scaled Cao construction, not a generic
concentrated profile.

### B.3 Equal-radius pair and nearby source-admissible data

Take two axial translates of this same scaled carrier with centers
\((z_1,r_0)\), \((z_2,r_0)\), separation \(d=|z_2-z_1|\), and the same
intensity \(a\). Their limiting source trajectories are

\[
 (z_j(t),r_0)=
 \left(z_j+\frac{a}{4\pi r_0}t,r_0\right),
\]

so their limiting-center relative velocity is exactly zero and their center
separation stays \(d\). Choose any \(R\) with

\[
 0<R<\min\{r_0,d/2\}.
\]

Buttà–Marchioro Theorem 1.1 then supplies
\(\rho_R>0\) and \(T_R>0\), independent of
\(\rho=\rho_\delta\), such that for every member of the constructed sequence
with \(\rho_\delta\le\rho_R\)

\[
 \boxed{
 \operatorname{supp}\xi_{j,\rho}(t)
 \subset\Sigma((z_j(t),r_0)\mid R),\qquad
 0\le t\le T_R.}
\]

The conclusion holds not only for the two exact profiles. Any nonnegative
same-leaf rearrangements \(\xi_{j,\rho}^0\) whose supports stay inside the
initial \(\rho\)-disks have the same circulation and the same relative-vorticity
\(L^\infty\) norm. Since \(r\le r_0+\rho\) there, their azimuthal vorticity
obeys the source bound with one enlarged, \(\rho\)-independent constant \(M\).
The small diffeomorphic deformations from Route A fit because the carrier has
the strict \(\rho/2\) support margin.

Therefore the complete supports have the exact coarse bounds

\[
 r\le r_0+R,\qquad
 \operatorname{gap}_z\ge\ell_R:=d-2R>0,
\qquad 0\le t\le T_R.
\]

**Route B1 verdict — established as stated.** Equal-radius,
equal-intensity Euler-scaled Cao rings, including a genuine nearby
same-rearrangement-leaf class, satisfy every hypothesis of the
Buttà–Marchioro full-support theorem and retain a positive support gap for a
concentration-independent time \(T_R\).

## Route B2 — splice the source interval into the deficit modulus

For a fixed admissible \(\rho\), let
\(\bar c_\rho=A_\delta c_\delta\),
\(\bar S_\rho\) be the scaled maximum, and
\(\bar m_{\rho,r_0+R}(\eta)>0\) its Cao compactness modulus. Define the actual
initial quantities

\[
 \delta_{0,\rho}=\Delta_{1,\rho}(0)+\Delta_{2,\rho}(0),\quad
 P_{T,\rho}=P(\xi_{1,\rho}^0)+P(\xi_{2,\rho}^0),\quad
 E_{12,\rho}^0=E_{12}(\xi_{1,\rho}^0,\xi_{2,\rho}^0).
\]

The common coefficient \(\bar c_\rho\), total energy conservation, and actual
total impulse conservation give

\[
 \Delta_{1,\rho}(t)+\Delta_{2,\rho}(t)
 =\delta_{0,\rho}+E_{12,\rho}(t)-E_{12,\rho}^0.
\]

The source support gap gives

\[
 E_{12,\rho}(t)
 \le\frac{C_GP_{T,\rho}^2}{(d-2R)^3}.
\]

Hence the fully explicit source-spliced condition

\[
 \boxed{
 \delta_{0,\rho}
 +\frac{C_GP_{T,\rho}^2}{(d-2R)^3}
 -E_{12,\rho}^0
 <\bar m_{\rho,r_0+R}(\eta)}
\]

implies

\[
 \boxed{
 d_{r_0+R}(\xi_{j,\rho}(t),\bar\Sigma_\rho)<\eta,
 \qquad j=1,2,\quad 0\le t\le T_R.}
\]

For each fixed source-admissible carrier this again defines a nonempty
relative-open neighborhood once the exact pair is given enough axial
separation: the cross-energy upper remainder tends to zero as \(d\to\infty\),
while the fixed-carrier modulus is strictly positive. Increasing \(d\) only
improves the source proof's separated-kernel bounds (with the same coarse
radius \(R\)); the strict inequality then persists under the small
nontranslation diffeomorphisms of Route A.

The source interval is genuinely stronger in the concentration limit than
the direct displacement certificate. Indeed, on the scaled Cao family,

\[
 \kappa_\rho\asymp|\log\rho|^{-1},\qquad
 P_{T,\rho}\asymp|\log\rho|^{-1},\qquad
 \Lambda_\rho\asymp
 \frac{1}{\rho^2|\log\rho|}.
\]

Thus the explicit Route-A velocity majorant satisfies

\[
 U_A\asymp\frac{\rho^{-4/3}}{|\log\rho|},
\qquad
 T_A\asymp\rho^{4/3}|\log\rho|\longrightarrow0,
\]

whereas the source returns \(T_R>0\) independent of \(\rho\).

**Route B2 verdict — established conditionally as stated.** On every fixed
scaled-Cao member satisfying the displayed strict deficit inequality, the
full orbital deformation estimate persists throughout the source's
concentration-independent support interval \(T_R\), with actual
\(\delta_{0,\rho}\), \(P_{T,\rho}\), and \(E_{12,\rho}^0\) retained.

## Route B3 — iterate the source theorem to all time

The theorem cannot be restarted at \(T_R\) with its original small parameter:
its output places the complete support in a fixed radius-\(R\) disk, not in a
new radius-\(\rho\) disk, and supplies no renewed
\(M/(\rho^2|\log\rho|)\) concentrated initial datum. The proof itself chooses
a positive time from a support-growth estimate; it does not state a
semigroup-uniform recurrence of the hypotheses.

**Route B3 verdict — blocked with missing construction.** All-time
equal-radius persistence needs a recurrent thin-support estimate, a
relative-motion invariant tube, or a pair-orbit compactness theorem. Finite
iteration of Theorem 1.1 is not licensed.

## Strongest result and next executable dependency

The strongest 0012 construction is the intersection of Routes A and B:

1. a genuine relative-open same-leaf neighborhood with small positive initial
   deficits and exact invariant bookkeeping;
2. an exact Euler scaling that puts the named Cao polynomial carrier and that
   neighborhood into the equal-radius Buttà–Marchioro class; and
3. full support separation plus qualitative orbital deformation control for
   a positive time \(T_R\) independent of the concentration parameter.

The next executable dependency is the all-time exit: use the exact
García–Hassainia–Hmidi periodic pair as the base orbit, define the pair
rearrangement class and conserved translating-frame functional, and prove
compactness of maximizing sequences modulo the pair's joint phase/translation
symmetries. That route would replace repeated single-ring tubes by one
pair-orbit modulus; it must first verify that the published periodic patches
are maximizers in that class.

Finally, every field here is axisymmetric without swirl:

\[
 u=u_re_r+u_ze_z,\qquad
 \boldsymbol\omega=\omega_\theta e_\theta,\qquad
 u\cdot\boldsymbol\omega=0.
\]

Thus 0012 supplies neither helicity nor a nontrivial axial rotation orbit and
licenses no quantum spin or spin-\(\tfrac12\) statement. It also makes no
full-three-dimensional perturbation, all-time pair, restoring-frequency,
electron, neutrino, parent-completion, or global no-go claim.
