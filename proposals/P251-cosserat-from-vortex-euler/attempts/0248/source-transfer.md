# Source-transfer audit and the common-frame repair

This is the substantive 0248 construction record.  It was written by the
fresh non-author `herdr geometry-review pane w3:p4` after the coordinator
added the central 0248 entry and the schema receipt (`WORKFLOW VALID: 269
claims`, exit 0).  It does not re-review 0245 or the completed periodic join.

## 1. Exact licenses found in the archived bodies

| Input | Exact reusable statement | Boundary that matters here |
| --- | --- | --- |
| `0211/same-ring-beltrami-torus.md`, independently reviewed in `0215/review.md` | For each sufficiently large finite major radius (R), one smooth whole-space stationary ring has compact toroidal vorticity, actual far velocity (-U_R n), (U_R=\Gamma_\delta\log(R/a_{\rm ref})/(4\pi R)+O(R^{-1})>0), and a literal nonzero constant-curl inner solid torus with strict flux-action twist and an elliptic closed core. | The field is stationary in that particular translating frame.  Its velocity and pressure are not compactly supported, and the review does not license nonlinear superposition or a stationary array. |
| Enciso--Peralta-Salas, *Acta Math.* 214 (2015), Theorem 7.6 and equations (7.28)--(7.32), archived as `sources/1210.6271.pdf.txt` | An analytic measure-preserving return map with a Diophantine invariant boundary and nonzero normal torsion retains a nearby invariant circle after the local Moser identification of the two positive section measures. | This is a local persistence theorem.  It preserves a tube already carried by a nearby divergence-free field; it neither constructs the nearby stationary Euler field nor changes a noncontractible periodic core into a Euclidean ring. |
| The same primary source, Theorem 8.3 | A local constant-λ Beltrami field on a compact set with connected complement can be approximated in any fixed (C^k) norm by a decaying global constant-λ field. | This is a same-curl geometry approximation.  It does not preserve the generalized-Beltrami exterior of 0211, construct a positive-density stationary array, or transfer Euler/Lin histories by itself. |
| `0236/periodic-image-forcing.md` | For a compact-vorticity ring in a cubic cell, the periodic image field has leading uniform velocity (-2I/(3\rho P^3)).  After it is removed, the local velocity and gradient costs are (O(|\Gamma|R^4/P^5)) and (O(|\Gamma|R^3/P^5)).  On the circular filament the displayed quartic contribution is also a uniform axial drift. | Uniform image drift is a frame row, not a deformation.  The remaining source still requires a zero-frequency stationary inverse; the positive-frequency optical inverse is not that inverse. |
| `0204/finite-ring-transfer.md` and `0204/density-preserving-packing.md` | Full-pressure finite-window Euler/Lin transfer and exact positive-density disjoint packing are proved for the 0198 Constantin--La--Vicol compact-velocity annular shell. | The marked response lies in a quiet cavity.  This is not a transfer theorem for the nonzero 0211 inner core, and disjoint support is unavailable for 0211's long-range velocity. |
| `0226/density-normalization.md` and `0226/compact-pressure-recursion.md` | On the straight 0211 literal-curl supplier the phase and energy forms are positive, the nonlinear normalization (H=\nu\beta) can be solved exactly, and the scaling (a=\lambda^{-1/2}), (R=\bar R/\lambda), (P=dR) keeps (M_{\rm phase}/P^3\) finite and positive.  Finite full-pressure recursion rows can be constructed. | The document expressly leaves the stationary array and the curvature/error diagonal open.  Its available (O(\log\bar R/\bar R)) curved error is not (o(K^2)) under the proposed density scaling. |
| `0169/localization-and-interaction.md` and the Constantin--La--Vicol source receipt | If an already stationary field has advected pressure, (U=\chi(p)u), (P'(p)=\chi^2) is an exact compact localization.  The CLV construction instead solves the extra speed-profile equation and gives compact annular shells; disjoint copies then sum exactly. | The source does not permit prescribing the 0211 Lundquist core and the 0200 annular mode independently.  The extra equation changes the core jets, and the proved compact shell excludes its nonsmooth central circle. |
| `0145/global-field.md`, `0147/toroidal-transfer.md`, and `0153/same-field.md` | Same-curl Herglotz/rational-mode approximation can preserve fixed local topology and finite-time packet margins with full pressure; 0147 attaches an actual finite-action packet to one unknotted torus; 0153 can make a prescribed finite insertion have small normalized bulk energy. | These are not bare-existence gaps.  They do not supply the completed 0246 coupled history/current jet on the compact ring.  In particular 0153 explicitly withholds a uniform (C^2) Bloch derivative as its period grows and a fixed positive density in that limit. |
| `0057`, `0059`, and the Gaussian primary source `source-2006.15033` | The Gaussian Beltrami law is smooth and stationary, has full local support, and has positive spatial density of robust knotted-torus patches; the campaign already constructed stationary action density and a double-divergence repair for selected compact force columns. | Full local support does not make a nonlocal Euler response local.  These records do not transfer the completed C016 coupled histories/current to the 0211 ring, so citing the Gaussian theorem again would not close the present same-field join. |

The decisive source mismatch is therefore exact: 0204 owns a complete
finite-window transfer on the compact annular/quiet-cavity carrier, while
0211 owns the nonzero constant-curl Euclidean core.  No reviewed file gives
the implication from one carrier to the other.

## 2. The balanced-orbit candidate loses every ring's self frame

Write the rigidly oriented 0211 seed in its actual steady frame as

\[
 u_m=-U_R n_m+w_m,\qquad \omega_m=\operatorname{curl}w_m,
 \qquad [u_m,\omega_m]=0 .
\]

Here ([v,\omega]=(v\mathbin\cdot\nabla)\omega-
(\omega\mathbin\cdot\nabla)v).  Since (n_m) is constant, the isolated
identity is

\[
 [w_m,\omega_m]=U_R\,\partial_{n_m}\omega_m.                 \tag{1}
\]

Now take the preregistered balanced orbit, with disjoint vorticity cores and
(\sum_m n_m=0).  Summing the *full* seed velocities cancels their uniform
parts and leaves (W=\sum_m w_m).  In the (m)-th core its exact transport
residual is

\[
 [W,\omega_m]
 =U_R\partial_{n_m}\omega_m
   +\Big[\sum_{\ell\ne m}w_\ell,\omega_m\Big].              \tag{2}
\]

Thus total-impulse or orientation cancellation is not local frame
cancellation.  Projecting (2) onto the translation neutral row exposes the
scale mismatch without estimating a shape inverse.  At separation (P=dR),

\[
 U_R=\frac{\Gamma_\delta}{4\pi R}\log(R/a_{\rm ref})+O(R^{-1}),
 \qquad
 |v_{\rm distant}|=O\!\left(\frac{|\Gamma_\delta|}{d^3R}\right). \tag{3}
\]

For fixed admissible (d), the first term in (2) is larger in its translation
row by a factor of order (d^3\log(R/a_{\rm ref})).  The distant interaction
cannot be treated as the small correction which restores the missing
(-U_Rn_m); it is parametrically too small.  Signed copies or a cubic group
do not alter (1).

**Route verdict for preregistered Candidate A:** refuted as stated, by the
local self-frame deficit (2).  This is not a no-go for stationary ring
arrays, and none is inferred.

## 3. Failure-derived Candidate A1: one orientation and one common frame

The minimum representation repair is to stop cancelling frames that the
individual rings require.  Put one 0211 orientation class in a periodic
Bravais cell, or equivalently begin with one ring and its parallel periodic
images, and retain one global mean velocity (-V_Pn).  Let the periodic
image field near the reference core be

\[
 u_{\rm im}=c_Pn+\widetilde u_{\rm im}.
\]

Choose

\[
 V_P=U_R+c_P.                                               \tag{4}
\]

Then the trial field (-V_Pn+w_R+u_{\rm im}) has exactly the isolated
translation row (-U_Rn+w_R).  Its residual is only

\[
 [\widetilde u_{\rm im},\omega_R],                         \tag{5}
\]

not the logarithmically larger term in (2).  The 0236 Green expansion gives

\[
 |\widetilde u_{\rm im}|
   \le C|\Gamma|R^4/P^5,\qquad
 |\nabla\widetilde u_{\rm im}|
   \le C|\Gamma|R^3/P^5.                                  \tag{6}
\]

For the circular leading ring the explicitly computed quartic term is also
absorbed into (c_P), so the first displayed anisotropic forcing is later;
(6), which retains finite-core corrections, is the safe bound.  Relative to
the self-frame speed, its velocity scale is

\[
 \frac{|\widetilde u_{\rm im}|}{U_R}
   =O\!\left(\frac{1}{d^5\log(R/a_{\rm ref})}\right).       \tag{7}
\]

Equations (4)--(7) are a concrete analytic advance: they identify a genuinely
small stationary source after the correct neutral parameter is used.  A Haar
rotation of the *whole* oriented field and lattice restores isotropy of the
law; an orientation-balanced molecule inside each realization is neither
needed nor licensed.  The uniform background contributes
(\rho V_P^2/2) to the physical kinetic-energy density and must remain in
the mean/material action.  Whole-field time reversal changes the common
velocity and circulation together.

The exact next equation is nevertheless not yet solved.  If
(\mathcal F(\omega,V)= [u(\omega,V),\omega]), the required theorem is a
bordered zero-frequency inverse for

\[
 D_\omega\mathcal F(\omega_R,U_R)\,\delta\omega
      +D_V\mathcal F(\omega_R,U_R)\,\delta V
   =-[\widetilde u_{\rm im},\omega_R],                     \tag{8}
\]

on three-dimensional periodic perturbations, after quotienting translation,
rotation and circulation/profile kernels.  The axisymmetric bordered inverse
used to create the isolated 0211 ring does not cover (8), because the cubic
image field is not axisymmetric.  The positive-frequency contour inverse of
0220 also does not cover frequency zero.  Flat vorticity edges require the
weighted/tame mapping property named in 0236.

Pressure is another explicit row of the same repair, not a postprocessing
formality.  Once (8) gives an exact divergence-free (u_P,omega_P) with
([u_P,omega_P]=0),

\[
 \operatorname{curl}(u_P\times\omega_P)=0.
\]

One must additionally kill the three harmonic periods of
(u_P\times\omega_P) in the border.  Only then does a periodic Bernoulli
function (B_P) exist with

\[
 \nabla B_P=u_P\times\omega_P,
 \qquad p_P=B_P-\tfrac12|u_P|^2,                           \tag{9}
\]

which is the full stationary pressure including all image and correction
cross terms.  A solution of the vorticity equation without these period rows
would not yet be the claimed periodic pressure pair.

## 4. Density and normalization costs of A1

For one core circle of actual radius `r_c(d)` and one positive-radius
axisymmetric tube with meridional section `D_T(d)` in a cubic cell of side
`P=dR`, the exact geometric densities are

\[
 \ell_{\rm core}=\frac{2\pi r_c(d)}{P^3},\qquad
 \phi_T=\frac{V_T(d)}{P^3},\qquad
 V_T(d)=2\pi\!\int_{D_T(d)}(R+x)\,dx\,dz.                \tag{10}
\]

Thus `ell_core = 2*pi/(d^3 R^2)[1+o(1)]` and
`V_T=2*pi*R*A_T[1+O(a/R)]` for a core of normal scale `a`; the exact
weighted volume, rather than unweighted normal area, is what enters the
density.

If (j_\ell(d)>0) is the continued measured coefficient per physical core
length and (g_{\rm law}>0) is the already declared whole-law detector Gram,
then

\[
 j_{\rm vol}(d)
  =g_{\rm law}\,f_{\rm tag}(d)\,
    \frac{J_{\rm ring}(d)}{P^3}>0,
 \qquad J_{\rm ring}=2\pi r_c j_\ell                    \tag{11}
\]

for any one fixed finite construction on which persistence and the response
margin prove (A_T,j_\ell,f_{\rm tag}>0).  No (d\to\infty) limit is used.
For the 0226 density scaling,

\[
 a=\lambda^{-1/2},\quad R=\bar R/\lambda,\quad P=dR,
 \quad M_{\rm phase}/P^3\asymp\rho\lambda^2a^4\asymp\rho, \tag{12}
\]

while the positive tag fraction is (O(\lambda)), the tagged mass density
vanishes like the inertia density divided by (R^2), and the background
energy density is

\[
 e_{\rm bg}\asymp
 \rho\Omega_0^2/(\lambda^2\bar R^2)+\rho V_P^2/2.          \tag{13}
\]

Here the first term denotes the zero-mean periodic ring part and the second
the retained common mean frame; their kinetic cross is zero in the
zero-mean periodic Biot--Savart gauge.

These are distinct physical costs.  Equation (12) is not permission to call
the unconstructed array stationary, and a dilute tagged mass is not the
ambient translational mass `rho`.

Put `S_im=[tilde u_im,omega_R]` in the actual weighted source space `Y`.
The core bounds and (6) give
`||S_im||_Y <= C_core,k(R) |Gamma| R^4/P^5`, with the gradient term in
(6) retained inside `C_core,k`.  If the inverse in (8) has norm
`C_st(R,d)`, then the correction and topology/observation losses have to be
budgeted as

\[
 \|\delta u\|_{C^k(\mathcal T)}
 \le C_{\rm st}(R,d)\|S_{\rm im}\|_Y,
 \qquad
 |j_\ell(d)-j_\ell(\infty)|
 \le C_{\rm obs}\|\delta u\|_{C^k(\mathcal T)}.           \tag{14}
\]

The fixed tube persists only when the first quantity is below its reviewed
0211/0215 (C^k) margin; measured positivity requires it also be below
(j_\ell(\infty)/(2C_{\rm obs})).  Pressure costs follow from the same
solution and (9), not from a compact-support assumption.

## 5. Candidate B and the dynamical transfer do not disappear

Candidate B remains blocked, rather than refuted.  The exact localizable
axisymmetric system contains both Grad--Shafranov balance and

\[
 |\nabla\psi|^2+F(\psi)^2=2r^2A_{\rm speed}(\psi).         \tag{15}
\]

The 0169/CLV source proves a particular smooth annular shell after excluding
the nonsmooth central circle.  It does not give a free profile theorem that
matches the 0211 Lundquist ratio to the 0200 constraint
(W/\Omega=s/\sqrt2), with a flat compact exterior and all derivatives.
The minimum repair for B is still an actual solution of (15) and the companion
elliptic equation with those inner, transition, annular and outer data.  A
formal buffer or a citation to compact Euler existence is insufficient.

Even a successful A1 stationary inverse would not by itself close the parent
dynamical row.  The audit found two noninterchangeable transfer licenses:

1. 0204 controls the full pressure and histories on the compact annular
   quiet-cavity ring, not on the 0211 nonzero core;
2. 0226 supplies the straight 0211 action/current normalization, but its
   curved (O(\log\bar R/\bar R)) error is not (o(K^2)) on the density
   diagonal.

The minimum same-ring dynamical repair after (8) is therefore to finish the
0226 compact pressure recursion through the finite parent spatial order on
the curved 0211 carrier, including its exterior pressure moments, material
angle/spin/current rows and the independent quadratic action normalization,
and then add the stationary correction loss (14) to that same error diagonal.
The completed C016 periodic join cannot simply be substituted: its local
normal dynamics and its noncontractible core are different source objects.

## 6. Strongest supported statement and supplier license

The audit establishes the exact failure mechanism of the balanced molecule,
the correct common-frame reduction (4)--(7), the full pressure-period rows
(9), and the physical density/normalization ledger (10)--(14).  It identifies
a small forcing for a genuinely new one-orientation periodic-ring candidate;
it does **not** prove the zero-frequency inverse (8) or the curved same-ring
dynamical transfer.

`route_verdict`:

- balanced Candidate A: **refuted as stated** by (2)--(3);
- failure-derived Candidate A1: **blocked on the named bordered
  zero-frequency inverse (8)**, with its correct small source established;
- Candidate B: **blocked on the common localizable profile theorem (15)**.

`evidence_scope`: exact source-transfer audit, local-frame calculation,
periodic-image/frame reduction, pressure-period specification, and
finite-density cost accounting.  No production numerics and no exhaustion
claim.

The precise parent supplier license is currently negative-edged but useful:

> Issue200 may import from 0248 that an impulse-balanced orbit of 0211 rings
> is not a valid small-interaction starting point, and that the viable
> stationary continuation starts with one parallel orientation class and a
> common frame (V_P=U_R+c_P), leaving only the image shape source (6).
> It may also import the density and pressure rows (9)--(14) as the contract
> for that continuation.  It may not yet import an exact stationary compact
> Euclidean array, a transferred coupled history/action/current, or the
> compact Euclidean geometry row of the parent theorem.

Next construction: prove (or replace) the weighted bordered inverse (8) for
the one-orientation periodic ring, then execute the finite 0226 curved
pressure/current recursion on that corrected field.  Failure of that inverse
would generate another stationary representation; it would not establish
obligation-level exhaustion.

## 7. Supplemental 0251 source-license boundary

The later root-owned `0251/source-transfer.md` adds two compatible primary
licenses.  They refine the continuation map but do not change any 0248 route
verdict.

Peralta-Salas--Slobodeanu, arXiv:2606.13462v1, Theorem 1.1 concerns an
analytic steady velocity in a bounded domain with a `C1` extension, velocity
tangent to the smooth boundary, advected pressure `u dot grad p=0`, constant
Bernoulli on every boundary component, and nonzero boundary Bernoulli
gradient.  Under those full hypotheses its conclusion is axisymmetry and a
toroidal domain with convex poloidal boundary curves.  Candidate B was
preregistered as axisymmetric, so this theorem is compatible structural
guidance for its allowable boundary geometry.  It does not refute B, produce
the inner-Lundquist/outer-annular matching profile, establish compact support,
or supply either steady zero-frequency inverse needed by the actual join.

Baldi, DOI 10.1007/s00021-023-00836-1, Theorem 1.1 gives an analytic
action-angle chart for a toroidal shell of the Gavrilov flow.  Its particle
frequencies satisfy the exact relations recorded in 0251,

\[
 \Omega_1=\chi(K(I))K'(I),\qquad
 \frac{\Omega_2}{\Omega_1}=\sqrt I\,R(I),                 \tag{16}
\]

with strictly increasing frequency ratio, and its action measures enclosed
fluid volume up to a fixed constant.  This is a valuable fixed-geometry
transport-frequency license for 0250.  It is not a field-changing Euler/Lin
mode: it supplies neither the perturbation pressure feedback nor the physical
action/current observation gain.  Consequently it cannot be substituted for
0248's stationary matching/inverse equations (8), (9), or (15), and it cannot
by itself transfer the completed coupled histories to the compact flow.

The resulting typed separation is:

- the 2026 symmetry theorem constrains the geometry and boundary hypotheses
  of Candidate B;
- Baldi supplies an exact passive-particle action-angle/frequency chart for
  the fixed compact carrier used by 0250;
- 0248 still requires an exact stationary matching or bordered steady inverse;
- the parent still requires the subsequent full Euler/Lin pressure,
  action/current and measured-output transfer on that same field.

These are complementary licenses, not competing verdicts and not a completed
source-to-parent implication.
