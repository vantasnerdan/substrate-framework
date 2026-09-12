# P253/0002 — exact R1 carrier/restoring/interaction transfer

## Frozen boundary and notation

This attempt addresses only P0/P2/P3 on R1 at base `v0.183.0` / `b6fc902`.
It neither changes the electron-first ordering nor reduces the conjunctive electron-and-neutrino
parent objective. It uses no target history, fitted particle datum, numerical stability claim, or
surrogate force law.

All four papers use constant-density incompressible Euler on \(\mathbb R^3\), restricted to
axisymmetric flow without swirl. Write cylindrical coordinates as \((r,\theta,z)\), azimuthal
vorticity as \(\omega_\theta\), and relative or potential vorticity as
\(\zeta=\omega_\theta/r\). Then

\[
  \partial_t\zeta+K[\zeta]\cdot\nabla\zeta=0,
  \qquad
  K[\zeta]=\frac{-\partial_zG\zeta}{r}e_r+
              \frac{\partial_rG\zeta}{r}e_z .
\]

In the Choi/Cao Green-operator normalization, the source functionals used below are

\[
 \Gamma[\zeta]=\int_{\mathbb R^3}\zeta\,dx,
 \quad
 P[\zeta]=\frac12\int_{\mathbb R^3}r^2\zeta\,dx,
 \quad
 E[\zeta]=\frac12\int_{\mathbb R^3}\zeta\,G\zeta\,dx .
\]

The source PDFs and their hashes are pinned in `source-provenance.yaml`.
Full downloaded papers and extracted text are retained outside the campaign tree in the
ephemeral primary-source cache named there. This artifact uses original derivations and
paraphrase rather than extended source excerpts.

## Source-theorem table

| Source result | Exact object and hypotheses | Exact conclusion transferable to P253 | Boundary of transfer |
|---|---|---|---|
| Choi, arXiv:2011.06808v2, Theorems 1.1–1.2 and Remark 1.5 (pp. 4–5, 14) | Hill relative vorticity \(\zeta_H^{\lambda,a}=\lambda 1_{B_a}\); nonnegative axisymmetric no-swirl data. The non-patch theorem assumes \(\zeta_0,r\zeta_0\in L^\infty\) and small \(L^1\cap L^2\) plus \(r^2L^1\) distance. | For every fixed \(\lambda,a>0\), the translating Hill vortex is orbitally stable for all \(t\ge0\), modulo axial translation, in the stated weighted vorticity metric. The cohesion mechanism is the energy-maximizer/impulse/conserved-distribution argument, not a Hooke force. | Axisymmetric, swirl-free, nonnegative relative vorticity only; the orbit removes axial translation only. The carrier is a spherical patch, not either thin-ring interaction family. No pair interaction or full 3-D perturbation theorem transfers. |
| Cao–Lai–Qin–Zhan–Zou, arXiv:2206.10165v2, Proposition 1.1, Theorems 1.2, 1.6, 1.8 (pp. 4–7) | Axisymmetric relative vorticity; the carrier/maximizer and Theorem 1.2 classes are nonnegative. For the explicit family: \(\kappa,W>0\), \(p\ge2\), sufficiently small \(\varepsilon\), compact thin core, and polynomial profile \(s_+^p\). Initial perturbations in Theorem 1.8 lie in \(L^1\cap L^\infty\cap L_w^1\), have \(r\omega_0\in L^\infty\), and are close in \(L^1\cap L^2\) plus impulse. | A family \(\zeta_\varepsilon^{(p,\kappa,W)}\) of classical steadily translating thin rings exists, is unique up to axial translation in the theorem's asymptotic class, and is nonlinearly orbitally stable for all positive time in the declared metric. The restoring mechanism is constrained maximization of \(E-W\log(1/\varepsilon)P\), uniqueness of the maximizing orbit, and conservation/compactness. | Still axisymmetric and swirl-free; no azimuthal or general 3-D modes. Stability is around one ring, not a two-ring relative periodic orbit. The theorem supplies no interaction dynamics, quantum identity, or selected \(\kappa,W,\varepsilon,p\). |
| Dávila–del Pino–Musso–Wei, arXiv:2207.03263v4, (1.14) and Theorem 1 (pp. 4–6) | A prescribed collision-free solution \((q_1,\ldots,q_k)\) of reduced system (1.14) on a fixed scaled interval \([0,T]\); sufficiently small \(\varepsilon\). The exact Euler vorticity is a smooth deformation of identical positive Kaufmann–Scully cores \(U(y)=8/(1+|y|^2)^2\). | For each such reduced path, an exact smooth axisymmetric no-swirl Euler solution exists whose ring centers and cores obey the stated asymptotics uniformly on \([0,T]\). Thus it rigorously transfers finite-window same-family interaction and its leading Hamiltonian law. | The Euler solution is globally defined, but the leapfrogging approximation is asserted only for the chosen finite \(T\); the authors explicitly do not assert long-time persistence or nonlinear stability. Because the reduced history is supplied to the theorem, this route cannot by itself serve P253's autonomous persistence mechanism. Its smooth decaying core is not Cao's compact polynomial-profile ring. |
| García–Hassainia–Hmidi, arXiv:2603.21644v1, Theorem 1.1 and Remarks 1.1–1.3 (pp. 6–7) | Fixed \(0<a<b\), \(\kappa>0\); sufficiently small \(\varepsilon\); \(\lambda\) in a Borel/Cantor set \(C_\varepsilon\subset(a,b)\) with \(|C_\varepsilon|\to b-a\); sufficiently large Sobolev index. Potential vorticity is the equal-strength patch pair \(\varepsilon^{-2}(1_{D_1(t)}+1_{D_2(t)})\). | Defines the named same-family interaction construction \(\mathcal G_{\varepsilon,\kappa,\lambda}\): two simply connected patch rings, related by a half-period time shift, giving a global exact Euler solution periodic in a translating frame. This is an autonomous all-time existence result for interacting carriers, not merely a finite-window reduced trajectory. | It is an existence theorem for selected periodic solutions, not Lyapunov/orbital stability of an initial-data neighborhood. The patch family is discontinuous and distinct from Cao's smooth \(p\ge2\) family and Dávila et al.'s Kaufmann–Scully family. Speed modulation and Nash–Moser invertibility construct the orbit; they do not establish dynamical restoration after arbitrary physical perturbations. |

## Exact calculations and useful transfers

### 1. Hill carrier invariants expose unsolved scale selection

For \(\zeta_H^{\lambda,a}=\lambda 1_{B_a}\), direct cylindrical integration gives

\[
 \Gamma=\lambda\,\mathrm{Vol}(B_a)=\frac{4\pi}{3}\lambda a^3,
 \qquad
 P=\frac{\lambda}{2}\int_{B_a}(x^2+y^2)\,dx
   =\frac{4\pi}{15}\lambda a^5.
\]

Choi's exact source formulas also give

\[
 W_H=\frac{2}{15}\lambda a^2,
 \qquad E_H=\frac{8\pi}{315}\lambda^2a^7.
\]

Hence stability holds leaf-by-leaf after choosing \((\lambda,a)\); it does not select either
parameter. In particular \(P/\Gamma=a^2/5\) and \(W_H=\Gamma/(10\pi a)\), so the source
contains a continuous scale/circulation family rather than a parameter-free particle mass scale.

### 2. The Cao family supplies the strongest carrier/restoring theorem

The family has exact circulation \(\kappa\), exact traveling speed
\(c_\varepsilon=W\log(1/\varepsilon)\), and concentration circle

\[
 r_*\longrightarrow\frac{\kappa}{4\pi W}.
\]

Eliminating \(W\) yields the source-level relation

\[
 c_\varepsilon=\frac{\kappa}{4\pi r_*}\log(1/\varepsilon)
 \quad\text{at the stated concentration limit.}
\]

This is a genuine relation among carrier variables, but \(\kappa\), \(\varepsilon\), and one
overall length/speed choice remain inputs. More importantly, orbital stability is a Lyapunov
statement: closeness to the axial-translation orbit remains small. It supplies no local spring
constant or internal oscillation unless a further quantitative second-variation theorem is proved.

### 3. The finite-window two-ring law has an exact reduced restoring geometry

For Dávila et al.'s symmetric two-ring reduction \(q_1=-q_2=q=(x,y)\), convention
\(q^\perp=(-y,x)\), equation (1.15) is

\[
 \dot x=\frac{2y}{x^2+y^2},\qquad
 \dot y=-\frac{2x}{x^2+y^2}-\frac{2x}{r_0^2}.
\]

It conserves

\[
 H(q)=-4\log(2|q|)-\frac{2x^2}{r_0^2},
 \qquad \dot H=0.
\]

Writing \(q=\rho(\cos\theta,\sin\theta)\) gives

\[
 \dot\rho=-\frac{\rho}{r_0^2}\sin(2\theta),
 \qquad
 \dot\theta=-\frac{2}{\rho^2}-\frac{2\cos^2\theta}{r_0^2}<0,
\]

and

\[
 \partial_\rho H=-\frac4\rho-\frac{4\rho\cos^2\theta}{r_0^2}<0.
\]

Thus every finite energy level is a unique smooth positive radial graph
\(\rho=\rho_H(\theta)\), and the nonvanishing vector field traverses it periodically. The
interaction logarithm and curvature term give exact bounded radial exchange in the reduced
law. This is a restoring geometry for the reduced centers, not a nonlinear stability theorem for
the Euler vorticity field. The source's error estimate transfers that law only on each fixed
\([0,T]\).

### 4. Exact same-family continuation from the stable Cao carrier

The direct published-theorem stitch fails below, so retain the stable family and build its
interaction without prescribing a future path. Fix one Cao member
\(\zeta_\varepsilon=\zeta_\varepsilon^{(p,\kappa,W)}\) and define the two-copy initial family

\[
 \mathcal C_{\varepsilon,d}(0)=
 \zeta_1(0)+\zeta_2(0),\qquad
 \zeta_j(0,r,z)=\zeta_\varepsilon(r,z-z_j),\quad |z_1-z_2|=d,
\]

with disjoint cores. This is nonnegative, compactly supported, axisymmetric no-swirl data in
the global source class. Transport the two initial labels by the single total Euler velocity.
Exactly, without a pair potential or target history,

\[
 \partial_t\zeta_j+
 \bigl(K[\zeta_1]+K[\zeta_2]\bigr)\cdot\nabla\zeta_j=0,
 \qquad j=1,2.
\]

Linearity and self-adjointness of the Green operator give the exact interaction decomposition

\[
 E[\zeta_1+\zeta_2]
 =E[\zeta_1]+E[\zeta_2]
  +\underbrace{\int\zeta_1G\zeta_2\,dx}_{E_{12}},
\]

and the exact cross-advection acting on ring \(j\) is
\(K[\zeta_{3-j}]\cdot\nabla\zeta_j\). Each material label conserves its circulation. For its
axial centroid \(Z_j=\Gamma_j^{-1}\int z\zeta_j\,dx\),

\[
 \dot Z_j=\Gamma_j^{-1}\int
   \bigl(K_z[\zeta_j]+K_z[\zeta_{3-j}]\bigr)\zeta_j\,dx.
\]

These identities are an exact full-Euler, same-smooth-family interaction construction. They
also expose the unclosed term: the cross velocity deforms each core, so the two translation
coordinates alone do not close. Cao's one-ring orbital theorem cannot simply be applied twice,
because an equal second ring is not a small perturbation of either one-ring state and the cross
energy changes the constrained variational problem.

Root's independent `attempts/0001/material-balances.md` now supplies the complementary exact
material identity behind this remainder: pressure moments and ambient pressure work drive shape,
spin, and internal energy. Its radial-swirl example proves that identical local tag moments can
have different accelerations because of the ambient pressure field. Therefore the two-ring route
must retain the full cross-generated pressure/shape state; the centroid equations above cannot be
promoted as an autonomous two-coordinate mechanics by omission. No 0001 result is duplicated or
modified here.

`verify_exact_transfer.py` checks the Hill integrals, the reduced Hamiltonian/polar identities,
and the bilinear cross-energy expansion symbolically; its first-run output is preserved in
`verify-exact.stdout`.

## Preserved failed route and constructive continuation

### R1-A — direct theorem stitch

**Route verdict: refuted with mechanism (composition only).** The claim that the published
Cao stability theorem can be attached directly to either published leapfrogging theorem is
false at the theorem boundary. The carrier classes do not coincide:

- Cao: compact smooth polynomial-profile relative vorticity, \(p\ge2\), steady one-ring
  maximizing orbit;
- Dávila et al.: smooth positive Kaufmann–Scully cores plus a correction, with a finite-window
  prescribed reduced path;
- García et al.: constant-potential-vorticity patches on two moving domains, forming a selected
  global periodic orbit.

Regularity, vorticity distribution, support type, and one-ring versus pair-orbit variational
problems are load-bearing hypotheses. Shared axisymmetry, circulation sign, thinness, or a
filament limit does not identify the families at finite \(\varepsilon\).

### R1-B — exact all-time same-family interaction

**Route verdict: established as stated at existence scope.**
\(\mathcal G_{\varepsilon,\kappa,\lambda}\) is a named exact same-patch-family construction:
the two members are related by a half-period shift and interact through the full Euler contour
dynamics for all time in a translating frame. This closes the earlier finite-window existence
gap for that selected orbit. It does not close P2 because no neighborhood stability/restoring
theorem is supplied.

### R1-C — stable-family interaction continuation

**Route verdict: established as an exact construction, blocked at persistence.**
\(\mathcal C_{\varepsilon,d}\) keeps the Cao carrier family and derives its exact cross-energy,
cross-advection, and centroid hierarchy from Euler. The missing construction is a quantitative
two-ring modulation/coercivity theorem controlling deformations around this manifold (or,
alternatively, an orbital-stability theorem for the García periodic patch pair).

## Attempt verdict and next executable construction

`route_verdict: blocked` for R1 as a P2+P3 mechanism.

`evidence_scope: SOURCE_THEOREM_SCOPED`.

Strongest exact result: R1 contains (i) a globally orbitally stable smooth classical one-ring
family, (ii) exact finite-window interacting smooth same-family solutions, (iii) exact global
periodic interacting patch-pair solutions, and (iv) the exact autonomous two-copy Cao
interaction identities above. No cited theorem supplies nonlinear persistence and restoring
control for the same interacting family.

The next executable analytic construction is the **Cao two-ring modulated energy route**:

1. Freeze \(\mathcal M_{\varepsilon,d}=\{\zeta_\varepsilon(\cdot-z_1)+
   \zeta_\varepsilon(\cdot-z_2)\}\) with common \((p,\kappa,W,\varepsilon)\) and separated cores.
2. Localize the exact functional into two single-ring energy-impulse pieces plus \(E_{12}\),
   and impose orthogonality to the two axial translation zero modes.
3. Upgrade Cao's compactness stability to a quantitative constrained coercivity estimate on
   each localized rearrangement leaf; compute the first and second variations of \(E_{12}\)
   rather than postulating a force.
4. Derive modulation equations for \((z_1,z_2)\) and a retained shape remainder from the exact
   labeled transport equations. Prove a closed error bound on the deformation and identify a
   periodic or scattering orbit in the resulting controlled system.
5. Test perturbations that exchange impulse between the two labels and perturb off the
   symmetric submanifold. Success earns P2/P3 only at the demonstrated axisymmetric
   no-swirl scope; general 3-D stability remains a separate continuation.

Any later soft eigenvalue, force, or small cross-energy verifier must first instantiate the
`small-ratio-numerics` prescriptions. At this boundary the remainder is analytic, so no
production numerical claim is made.

P4, the electron observables in P5, and the neutrino observables in P6 remain wholly unearned;
no parent completion or global no-go follows from this attempt.
