# P253/0006 — autonomous same-family deficit control

## Boundary, variables, and exact source scopes

This attempt continues only P253 R1-C at P2/P3. It fixes one
Cao–Lai–Qin–Zhan–Zou smooth polynomial-profile ring and transports two axial
translates by one autonomous Euler velocity. It does not prescribe a center
history, discard deformation, infer a Hessian gap from uniqueness, or make a
numerical force/stability claim.

Use Cao et al.'s half-plane variables \(x=(r,z)\), weighted measure
\(d\nu=r\,dr\,dz\), Green operator \(G_1\), energy, and impulse

\[
 E(\zeta)=\frac12\int\zeta G_1\zeta\,d\nu,
 \qquad
 P(\zeta)=\frac12\int r^2\zeta\,d\nu .
\]

Buttà–Cavallaro–Marchioro instead write the azimuthal vorticity as
\(\omega=\omega_\theta\). The exact dictionary is

\[
 \zeta=\frac{\omega_\theta}{r},\qquad
 \int_{\Pi}\omega_\theta\,dr\,dz=\int_{\Pi}\zeta\,d\nu .
\]

The source statements used or tested are:

| Source result | Exact useful scope | Transfer in 0006 |
|---|---|---|
| Cao et al., Proposition 1.4 and Theorems 4.2–4.3 | For fixed \(p\ge2\), \(\kappa,W>0\), a smooth compact ring of cross-section \(O(\varepsilon)\), circulation \(\kappa\), speed \(W\log(1/\varepsilon)\), and a one-orbit maximizer set | Names the actual carrier and its maximizing orbit |
| Cao et al., Theorem 2.9 | Every maximizing sequence in the source admissible rearrangement class is strongly \(L^2(d\nu)\)-compact modulo axial translations | Gives a positive qualitative deficit modulus, not a quantitative Hessian gap |
| Cao et al., Lemma 2.11 and Theorem 2.10 | Relative vorticity transported by an admissible common axisymmetric no-swirl velocity retains its exact rearrangement class; the maximizing set is orbitally stable in the stated \(L^1+L^2\) topology (with the additional impulse term when maximizer impulses agree) | Applies label-by-label even though the labels exchange impulse |
| Cao et al., equation (2.6) | \(G_1(x,y)\le C_\delta(rr')^{\delta+1/2}/|x-y|^{2\delta}\), \(0<\delta\le3/2\) | At \(\delta=3/2\), bounds cross energy by the actual exchanged moments |
| Buttà–Cavallaro–Marchioro, published Theorem 1.1 | Definite-sign rings with support in \(\varepsilon\)-disks, intensity \(a_i/|\log\varepsilon|\), \(L^\infty\le M/(\varepsilon^2|\log\varepsilon|)\), and distinct radii \(|r_i-r_j|\ge2D\); for every prescribed finite \(T\), radial support localizes and the vorticity mass concentrates axially near independently translating centers as \(\varepsilon\to0\) | Does not admit two identical axial translates, whose radii agree |
| Buttà–Cavallaro–Marchioro, Theorem 2.2 and Lemma 2.3 | A single concentrated ring remains radially localized and axially mass-concentrated under an external field whose size and Lipschitz constant are \(O(|\log\varepsilon|^{-1})\) | The multi-ring reduction supplies that field only while the full supports remain separated |

The published ZAMP article is controlling. Its Remark 1.1 calls the
distinct-radius condition essential and explains that no axial
support-localization estimate prevents overlap at equal radii. The broader
sentence in arXiv v3 that suggested any separated limiting trajectories is
therefore not imported.

Fix \(p\ge2\), \(\kappa,W>0\), and a sufficiently small, fixed
\(\varepsilon\). Let \(\zeta_\varepsilon\) be the Cao ring, and set

\[
 c_\varepsilon=W\log(1/\varepsilon),\qquad
 \mathcal F_\varepsilon(\zeta)=E(\zeta)-c_\varepsilon P(\zeta).
\]

Let \(S_\varepsilon\) be the maximum of \(\mathcal F_\varepsilon\) over
\(\mathcal R(\zeta_\varepsilon)^w\). Theorems 4.2, 4.3, and 1.6 identify

\[
 \Sigma_\varepsilon=
 \{\zeta_\varepsilon(\cdot+a e_z):a\in\mathbb R\}
 \subset\mathcal R(\zeta_\varepsilon).
\]

The exact hypotheses used below are
\(\zeta_\varepsilon\ge0\), compact support,
\(\zeta_\varepsilon\in L^1\cap L^q\cap L^\infty\) for \(q>2\), finite
weighted impulse, and \(r\zeta_\varepsilon\in L^\infty\). Write

\[
 \int\zeta_\varepsilon\,d\nu=\kappa,\quad
 P(\zeta_\varepsilon)=P_\varepsilon>0,\quad
 \|\zeta_\varepsilon\|_\infty=\Lambda_\varepsilon,
\]

and choose its axial origin so that its support has radial maximum \(R_0\)
and axial half-width \(h_0\).

## Exact cross variations before stability

For two labeled fields define

\[
 E_{12}(\zeta_1,\zeta_2)
 =\iint\zeta_1(x)G_1(x,y)\zeta_2(y)\,d\nu_xd\nu_y.
\]

Self-adjointness gives

\[
 \frac{\delta E_{12}}{\delta\zeta_1}=G_1\zeta_2,
 \qquad
 \frac{\delta E_{12}}{\delta\zeta_2}=G_1\zeta_1.
\]

For frozen translated shapes
\(\zeta_j(r,z)=\widehat\zeta_j(r,z-a_j)\), axial translation invariance gives

\[
 \partial_{a_1}E_{12}
 =\iint\widehat\zeta_1(x)\,\partial_{z_x}G_1
       (r,z+a_1;r',z'+a_2)\widehat\zeta_2(y)\,d\nu_xd\nu_y,
 \qquad
 \partial_{a_2}E_{12}=-\partial_{a_1}E_{12},
\]

and

\[
 \partial_{(a_1-a_2)}^2E_{12}
 =\iint\widehat\zeta_1(x)\,\partial_{z_xz_x}G_1
       (r,z+a_1;r',z'+a_2)\widehat\zeta_2(y)\,d\nu_xd\nu_y.
\]

No sign is assigned to the second derivative. A restoring coefficient would
require its actual sign together with the self-field constrained Hessian and
deformation response.

## Route A — actual constrained Hessian coercivity

Let

\[
 \Psi_\varepsilon=G_1\zeta_\varepsilon-
                 \frac{c_\varepsilon}{2}r^2-\mu_\varepsilon,
 \qquad
 \zeta_\varepsilon=f_\varepsilon(\Psi_\varepsilon),
 \qquad
 f_\varepsilon(s)=\varepsilon^{-2}(s_+)^p.
\]

Take a smooth weighted-volume-preserving axisymmetric flow \(\Phi_s\),
generated by \(v\) with
\(\partial_r(rv_r)+\partial_z(rv_z)=0\), and the isovortical path
\(\zeta_s=\zeta_\varepsilon\circ\Phi_s^{-1}\). Its tangent is
\(\sigma=-v\cdot\nabla\zeta_\varepsilon\). Differentiating the actual
energy–impulse functional along this path gives the constrained quadratic form

\[
 Q_\varepsilon[\sigma]
 =\int\sigma G_1\sigma\,d\nu
  -\int_{\{\Psi_\varepsilon>0\}}
       \frac{\sigma^2}{f_\varepsilon'(\Psi_\varepsilon)}\,d\nu .
\]

For these source-admissible tangents the second integrand equals
\(f_\varepsilon'(\Psi_\varepsilon)(v\cdot\nabla\Psi_\varepsilon)^2\);
the boundary degeneracy is therefore not being treated as a free arbitrary
\(L^2\) direction. Maximality gives \(Q_\varepsilon\le0\) on such paths.
Axial translation has
\(\sigma_T=-\partial_z\zeta_\varepsilon\) and
\(Q_\varepsilon[\sigma_T]=0\).

A quantitative route would need a declared tangent completion and norm plus

\[
 -Q_\varepsilon[\sigma]\ge\gamma_\varepsilon\|\sigma\|_X^2
 \quad\text{for }\sigma\perp\sigma_T,\qquad \gamma_\varepsilon>0.
\]

Cao et al. prove compactness of maximizing sequences and uniqueness of the
maximizing orbit. They do not state this inequality, choose \(X\), or exclude
nontranslation soft sequences at a quantitative rate.

**Route A verdict — blocked with missing construction.** The exact constrained
quadratic form and translation zero mode are identified, but quantitative
coercivity still needs an independent spectral/functional inequality on a
completed source-admissible tangent space. Uniqueness is not used as a Hessian
shortcut.

## Route B1 — qualitative modulus plus exact cross-energy control

### Autonomous labeled Euler system

Set

\[
 \zeta_j(0,r,z)=\zeta_\varepsilon(r,z-a_j),\qquad
 d=|a_1-a_2|>2h_0,\qquad g_0=d-2h_0>0,
\]

and evolve

\[
 \partial_t\zeta_j+
 (K[\zeta_1]+K[\zeta_2])\cdot\nabla\zeta_j=0,\qquad j=1,2.
\]

The common flow transports the initially disjoint labels injectively while
the classical axisymmetric no-swirl solution exists. Lemma 2.11 applies to
each label under the total velocity:

\[
 \zeta_j(t)\in\mathcal R(\zeta_\varepsilon),\qquad
 \int\zeta_j(t)\,d\nu=\kappa,\qquad
 \|\zeta_j(t)\|_{L^s(d\nu)}
 =\|\zeta_\varepsilon\|_{L^s(d\nu)}.
\]

Individual impulse is not conserved. Since
\(K_r[\zeta]=-r^{-1}\partial_zG_1\zeta\), weighted transport and kernel
antisymmetry give

\[
 \dot P_1=-\iint\zeta_1(x)\partial_{z_x}G_1(x,y)\zeta_2(y)\,d\nu_xd\nu_y,
 \qquad
 \dot P_2=-\dot P_1.
\]

Thus the labels exchange impulse while
\(P_1(t)+P_2(t)=2P_\varepsilon\) is exact.

### Nonnegative deficits and exact sum identity

Define

\[
 \Delta_j(t)=S_\varepsilon-\mathcal F_\varepsilon(\zeta_j(t)).
\]

Each label remains in the exact rearrangement class contained in the weak
class over which \(S_\varepsilon\) is the supremum, so

\[
 \Delta_j(t)\ge0.
\]

Also

\[
 E(\zeta_1+\zeta_2)=E(\zeta_1)+E(\zeta_2)+E_{12},
 \qquad P(\zeta_1+\zeta_2)=P_1+P_2.
\]

Total Euler energy and impulse are conserved. Because both deficits use the
same \(c_\varepsilon\), the exchanged individual impulses cancel exactly:

\[
 \boxed{\Delta_1(t)+\Delta_2(t)
 =\Delta_1(0)+\Delta_2(0)+E_{12}(t)-E_{12}(0)}.
\]

Both initial labels are maximizers, hence

\[
 \boxed{\Delta_1(t)+\Delta_2(t)=E_{12}(t)-E_{12}(0)\ge0}.
\]

This retains impulse exchange and full shape deformation; no
collective-coordinate closure is used.

### Cross-energy bound with exchanged moments retained

Cao et al. equation (2.6) at exponent \(3/2\) supplies \(C_G>0\) with

\[
 0<G_1(x,y)\le C_G\frac{(rr')^2}{|x-y|^3}.
\]

If the material supports have axial gap at least \(\ell>0\), then

\[
 \begin{aligned}
 E_{12}(t)
 &\le \frac{C_G}{\ell^3}
       \left(\int r^2\zeta_1\,d\nu\right)
       \left(\int r'^2\zeta_2\,d\nu\right)\\
 &=\frac{4C_G}{\ell^3}P_1(t)P_2(t).
 \end{aligned}
\]

Nonnegativity and conserved total impulse give

\[
 P_1P_2\le\frac{(P_1+P_2)^2}{4}=P_\varepsilon^2,
 \qquad
 \boxed{E_{12}(t)\le\frac{4C_GP_\varepsilon^2}{\ell^3}}.
\]

No fictitious conservation of individual impulse enters the bound.

### Source-derived qualitative modulus

For \(R\ge R_0\), define

\[
 \begin{aligned}
 d_R(\zeta,\Sigma_\varepsilon)
 =\inf_{a\in\mathbb R}\big(&
 \|\zeta-\zeta_\varepsilon(\cdot+a e_z)\|_{L^1(d\nu)}
 +\|\zeta-\zeta_\varepsilon(\cdot+a e_z)\|_{L^2(d\nu)}\\
 &+P(|\zeta-\zeta_\varepsilon(\cdot+a e_z)|)\big)
 \end{aligned}
\]

and

\[
 m_{\varepsilon,R}(\eta)=\inf\left\{
 S_\varepsilon-\mathcal F_\varepsilon(\zeta):
 \begin{array}{l}
 \zeta\in\mathcal R(\zeta_\varepsilon),\
 \operatorname{supp}\zeta\subset\{r\le R\},\\
 d_R(\zeta,\Sigma_\varepsilon)\ge\eta
 \end{array}\right\}.
\]

For every \(\eta>0\),

\[
 \boxed{m_{\varepsilon,R}(\eta)>0}.
\]

If the infimum were zero, Theorem 2.9 would give, after axial translations,
strong \(L^2(d\nu)\) convergence to an element of \(\Sigma_\varepsilon\).
Equimeasurability bounds the measure of the union of the two supports by
\(2\nu(\operatorname{supp}\zeta_\varepsilon)\), so

\[
 \|\zeta_n-\zeta_*\|_1
 \le \sqrt{2\nu(\operatorname{supp}\zeta_\varepsilon)}
       \|\zeta_n-\zeta_*\|_2.
\]

Radial support in \(r\le R\) then gives

\[
 P(|\zeta_n-\zeta_*|)
 \le\frac{R^2}{2}\|\zeta_n-\zeta_*\|_1.
\]

This contradicts \(d_R\ge\eta\). The modulus is positive but has no assigned
rate and no Hessian interpretation.

### Derived autonomous interval

The three-dimensional vorticity is

\[
 \boldsymbol\omega=r(\zeta_1+\zeta_2)e_\theta .
\]

Set \(R=2R_0\). While both supports lie in \(r\le R\), disjoint label
transport gives

\[
 \|\boldsymbol\omega\|_\infty\le R\Lambda_\varepsilon.
\]

Cauchy–Schwarz and total mass/impulse conservation give

\[
 \begin{aligned}
 \|\boldsymbol\omega\|_1
 &=2\pi\int r(\zeta_1+\zeta_2)\,d\nu\\
 &\le2\pi\sqrt{(2\kappa)(4P_\varepsilon)}
 =4\sqrt2\pi\sqrt{\kappa P_\varepsilon}.
 \end{aligned}
\]

Splitting the three-dimensional Biot–Savart integral at radius \(b\) gives

\[
 \|u\|_\infty
 \le b\|\boldsymbol\omega\|_\infty+
 \frac{\|\boldsymbol\omega\|_1}{4\pi b^2}.
\]

Optimizing at
\(b^3=\|\boldsymbol\omega\|_1/(2\pi\|\boldsymbol\omega\|_\infty)\) yields

\[
 \|u\|_\infty\le
 C_{BS}\|\boldsymbol\omega\|_1^{1/3}
 \|\boldsymbol\omega\|_\infty^{2/3},
 \qquad
 C_{BS}=\frac{3}{2^{4/3}\pi^{1/3}}.
\]

Define

\[
 U_R=C_{BS}
 \left(4\sqrt2\pi\sqrt{\kappa P_\varepsilon}\right)^{1/3}
 (R\Lambda_\varepsilon)^{2/3},
\qquad
 T_*=\min\left\{\frac{R_0}{U_R},\frac{g_0}{4U_R}\right\}.
\]

A continuity/bootstrap argument gives, for \(0\le t\le T_*\), radial support
in \(r\le2R_0\) and axial gap at least \(g_0/2\). Consequently

\[
 0\le\Delta_j(t)\le\Delta_1(t)+\Delta_2(t)
 \le\frac{32C_GP_\varepsilon^2}{g_0^3}=:B(g_0).
\]

For any desired orbital tolerance \(\eta>0\), the source-level condition

\[
 \boxed{B(g_0)<m_{\varepsilon,2R_0}(\eta)}
\]

implies

\[
 \boxed{
 d_{2R_0}(\zeta_j(t),\Sigma_\varepsilon)<\eta,\quad
 j=1,2,\quad 0\le t\le T_*}.
\]

Equivalently, it suffices that

\[
 g_0>
 \left(\frac{32C_GP_\varepsilon^2}
 {m_{\varepsilon,2R_0}(\eta)}\right)^{1/3}.
\]

**Route B1 verdict — established as stated.** Under the exact source
admissibility and displayed separation condition, each interacting Cao label
stays within the requested qualitative orbital tolerance through the derived
positive interval \(T_*\). The statement retains individual impulse exchange
and all shape deformation.

## Route B2 — closing the separation exit with cross energy alone

The exact deficit identity gives \(E_{12}(t)\ge E_{12}(0)>0\). Thus, whenever
the supports remain axially ordered with gap \(\ell(t)>0\),

\[
 \ell(t)\le
 \left(\frac{4C_GP_\varepsilon^2}{E_{12}(0)}\right)^{1/3}.
\]

This excludes unbounded axial scattering of the ordered pair. It does not
prevent \(\ell(t)\downarrow0\): the cross-energy upper bound becomes weaker
at small separation, while the nonnegative deficits may grow.

**Route B2 verdict — refuted with mechanism.** Cross-energy positivity and
the moment upper bound cannot close the lower-separation exit. An all-time
pair tube needs additional deformation-aware dynamics, a lower/renormalized
interaction bound, or a stable exact pair construction.

## Route C1 — published arbitrary-\(T\) concentrated-ring theorem

For the published theorem, each azimuthal-vorticity label must satisfy

\[
 |\log\rho|\int\omega_{i,\rho}\,dr\,dz=a_i,\qquad
 \|\omega_{i,\rho}\|_\infty\le
 \frac{M}{\rho^2|\log\rho|},
\]

and the fixed ring radii must obey
\(|r_i-r_j|\ge2D\). For every prescribed \(T<\infty\), it then gives sharp
radial support localization and axial concentration of the mass near

\[
 (z_i(t),r_i)=
 \left(z_i+\frac{a_i}{4\pi r_i}t,r_i\right).
\]

The intensity mismatch with a fixed Cao ring is exact:

\[
 \int\omega_\theta\,dr\,dz
 =\int\zeta_\varepsilon\,d\nu=\kappa,
\]

whereas the ZAMP theorem requires \(a_i/|\log\rho|\). Euler amplitude scaling
\(\widetilde\zeta=\alpha\zeta_\varepsilon\) with
\(\alpha=a/(\kappa|\log\rho|)\) fixes this intensity and preserves a
one-ring traveling solution with speed \(\alpha c_\varepsilon\), but it
changes the fixed-\(\kappa\) rearrangement class used in Route B.

More decisively, the two rings frozen in 0006 are axial translates of the same
Cao member. Hence \(r_1=r_2\), violating the published condition
\(|r_1-r_2|\ge2D\). The published Remark 1.1 expressly retains that condition:
only radial support is sharply localized; axial tails are not, so equal-radius
full supports cannot be kept separated by the theorem. It also states that
inter-ring interaction vanishes in this scaling limit, so its limiting
translations contain no finite interaction law.

**Route C1 verdict — refuted for the frozen same-Cao-family transfer.** The
published arbitrary-\(T\) theorem is constructive for a different-radius,
vanishing-circulation family, but its exact radius and intensity hypotheses
exclude the two identical fixed-circulation axial translates in Route B.

## Route C2 — use the reduced theorem to close the equal-radius exit

Theorem 2.2 would apply to either label if the field generated by the other
were globally replaceable on its support by a divergence-free external field
with size and Lipschitz constant \(O(|\log\rho|^{-1})\). Lemma 2.3 constructs
that replacement only while the two full supports are separated. At equal
radii, Theorem 2.2 controls the radial support but only axial mass
concentration, not the axial support. Therefore using it to prove continued
support separation assumes the very property needed to construct the
external field.

**Route C2 verdict — blocked by a circular missing estimate.** The reduced
theorem cannot continue Route B past its separation exit without a separate
axial support or tail-interaction estimate for the equal-radius pair.

## Strongest result, next executable construction, and spin boundary

The strongest 0006 result remains Route B1: a non-Hessian,
qualitative-modulus persistence theorem on an explicit autonomous interval
for two copies of the actual Cao carrier, with exact deficit bookkeeping,
moment-controlled cross energy, impulse exchange, and full-field deformation.

The new source sharpens the next construction. Audit and transfer the
short-time predecessor identified by the ZAMP authors
(Buttà–Marchioro, arXiv:1904.04785), whose stated distinction is strong
localization of the support in both radial and axial directions. The
executable test is: place two Euler-scaled Cao polynomial profiles at equal
radius and fixed axial gap; verify the predecessor's exact multi-ring
hypotheses; then combine its support-gap interval with the Route B deficit
identity and qualitative modulus. If that source also requires radial
staggering, pivot to the exact García–Hassainia–Hmidi periodic pair and prove
a pair-orbit compactness modulus rather than inventing a Hessian gap.

All carriers here are axisymmetric without swirl:

\[
 u=u_r e_r+u_z e_z,\qquad
 \boldsymbol\omega=\omega_\theta e_\theta,\qquad
 u\cdot\boldsymbol\omega=0.
\]

Thus the kinetic helicity is zero whenever integrable, and axial rotations
act trivially on the axisymmetric field. No result in 0006 supplies a
nontrivial rotation orbit, helicity, intrinsic spin, or spin-\(\tfrac12\)
license. The attempt also earns no full three-dimensional perturbation
stability, all-time equal-radius pair stability, restoring frequency,
electron or neutrino identification, parent completion, or global no-go.
