# P253/0033: actual GHH low/center Hamiltonian and finite-core period map

## Preregistration and activation boundary

This README is the sole preregistration artifact for attempt 0033. Its
positive direction is robust control of the same García--Hassainia--Hmidi
(GHH) finite-core periodic pair, with the actual low and center physical
degrees of freedom settled before any infinite-dimensional reducibility
scheme is resumed.

Before substantive work, root must:

1. register 0033 centrally with the ownership and scope below;
2. validate the central proposal and this frozen contract; and
3. place `activation-schema.exit` containing exactly `0` in this directory.

Until then, no primary-source body, derivation, verifier, sampled coefficient,
period integration, multiplier calculation, or result artifact may be opened
or produced for 0033. After activation, this worker owns only append-only
files in `attempts/0033`. Root retains ownership of the central proposal,
source/API/test changes, memory, commits, review, promotion, and downstream
replay. Attempts 0030--0032 are separate carrier routes and are not inputs.

Accepted base: release `v0.183.0`, source baseline
`substrate-framework@b39f392f8feee21c93a6052a5f971bc2df862779`.
Target kind: a fixed finite-dimensional Hamiltonian/Floquet theorem embedded
in the actual full two-contour IVP. There is no empirical comparator and no
production numerical verifier at freeze.

## Frozen actual object

Fix one GHH pair at the exact quantifiers of source Theorem 1.1:
\(0<a<b\), \(\kappa>0\), sufficiently small \(\epsilon\), and
\(\lambda\in\mathcal C_\epsilon\). Put \(L=|\log\epsilon|\), use physical
time, and write the translating-frame period as

\[
 \mathcal T_\epsilon
 =\frac{T(\epsilon,\lambda,\kappa)}{L}.                     \tag{1}
\]

The actual labeled vorticity and translating Hamiltonian are

\[
 q_*(t)=\epsilon^{-2}
   \bigl(\mathbf 1_{D_{1,*}(t)}+\mathbf 1_{D_{2,*}(t)}\bigr),
 \qquad
 \mathscr H_c=E-c_\epsilon P.                               \tag{2}
\]

Both labels evolve by the same self-consistent Euler velocity. No prescribed
center history, target-history surrogate, annular wall, reversible
perturbation restriction, or half-period label constraint is permitted.

Use the actual P253/0021 full-contour space

\[
 X^s=H^s_0(\mathbb T)\oplus H^s_0(\mathbb T)
 \simeq\mathbb R^4_{\rm ctr}\oplus S_1^s\oplus S_2^s,       \tag{3}
\]

where the two individual area constraints are exact and only tangential
reparametrization is gauge. The centroid/shape chart must be constructed
without double-counting: physical first Fourier normal variations and center
motions form one coupled sector, and a source first-mode constraint may not
delete either label's physical center tangent.

Retain the exact invariant flag

\[
 G=\operatorname{span}\{v_{\rm ph},v_z\}
 \subset K=\ker d\mathscr H_c\cap\ker dP,
 \qquad \mathcal N^s_{\rm red}=K/G,                         \tag{4}
\]

together with normalized companions \(w_E,w_P\). The phase vector, common
axial translation, energy covector, impulse covector, and companion shears
must be carried through the finite block before the restriction to \(K\) and
quotient by \(G\). A unit multiplier in the unreduced flag is not
automatically a reduced shape multiplier.

The physical convention is fixed:

\[
 E_{\rm kin}^{\rm phys}=2\pi\rho_m E,\quad
 I_z^{\rm phys}=2\pi\rho_m P,\quad
 \Omega^{\rm phys}=2\pi\rho_m\Omega,\quad
 \mathscr A^{\rm phys}=2\pi\rho_m\mathscr A,                \tag{5}
\]

with material density \(\rho_m\) distinct from \(\rho=r^2/2\). The factor
\(2\pi\rho_m\) must appear in both the KKS matrix and Hamiltonian Hessian and
therefore cancels from the Hamiltonian vector field.

## Exact finite sector to be constructed

Let \(\mathcal F_{\rm ctr}\) denote all real two-label center/first-mode
directions plus \(v_{\rm ph},v_z,w_E,w_P\). For a cutoff \(J_0\ge2\), let

\[
 \mathcal W_j
 =\operatorname{span}_{\mathbb R}
 \{\cos(j\theta)e_a,\sin(j\theta)e_a:a=1,2\},
 \qquad
 \mathcal L_{J_0}
 =\mathcal F_{\rm ctr}\oplus\bigoplus_{j=2}^{J_0}\mathcal W_j. \tag{6}
\]

The value of \(J_0\) is not fitted. It is the smallest cutoff above every
center/shape resonance or indefinite Krein block exposed by the exact
finite-core calculation. Modes cannot be removed by imposing the annular
source's \(m\)-fold symmetry.

Choose a periodic real frame
\(E_L(t):\mathbb R^{N_L}\to X^s\) spanning (6) after its exact
centroid/area slice is derived. Its physical KKS and Hessian matrices are

\[
 W_L^{\rm phys}(t)
 =E_L(t)^*\Omega_{q_*(t)}^{\rm phys}E_L(t),\qquad
 L_L^{\rm phys}(t)
 =E_L(t)^*D^2\mathscr H_c^{\rm phys}(q_*(t))E_L(t).         \tag{7}
\]

Neither matrix may be inferred from frequency reality or uniqueness.
Starting from the boundary KKS definition, the calculation must reproduce
the already exposed shape block

\[
 W_j^{\rm phys}
 =-\,2\pi\rho_m\frac{\epsilon^4}{j}
   \operatorname{diag}(J_2,J_2),\qquad
 J_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix},                \tag{8}
\]

in the normalized real sine/cosine basis, then derive every center--center,
center--shape, and label-cross entry rather than assuming block diagonality.
The moving-frame generator is defined directly from the full Jacobi operator:

\[
 A_L(t)
 =E_L(t)^{-1}\bigl(A_*(t)E_L(t)-\dot E_L(t)\bigr),          \tag{9}
\]

only after the full-contour invariant-subspace bridge below makes the inverse
meaningful. In a nonconstant symplectic frame it must satisfy

\[
 \dot W_L+A_L^TW_L+W_LA_L=0.                               \tag{10}
\]

Equation (10), not a guessed \(A=-W^{-1}L\) in a moving chart, is the sign and
connection check.

## Filament block and first finite-core correction

The GHH source contains a four-dimensional canonical filament Hamiltonian.
Common axial translation and the conserved radial-center sum reduce it to a
one-degree-of-freedom Hamiltonian with a nonconstant periodic orbit. Attempt
0033 first derives its real unsymmetrized two-label variational system and
period map, including the directions suppressed by the source's label-delay
chart.

In canonical action--phase coordinates for that reduced periodic orbit, the
filament center variational return has the form

\[
 \begin{pmatrix}\delta\varphi(\mathcal T_0)\\
                 \delta I(\mathcal T_0)\end{pmatrix}
 =
 \begin{pmatrix}1&\mathcal T_0\Omega_0'(I_0)\\0&1\end{pmatrix}
 \begin{pmatrix}\delta\varphi(0)\\\delta I(0)\end{pmatrix}. \tag{11}
\]

The source's strict period twist can therefore produce a genuine filament
phase/energy-companion Jordan shear. This is not yet an instability of the
reduced finite-core pair: 0033 must show whether (11) disappears after
\(\delta I=0\) and phase quotient in (4), and whether finite-core
center--shape coupling creates a different reduced block.

Next derive, rather than prescribe, the first nonzero finite-core scale
\(\delta_\epsilon\) and expansion

\[
 A_L^\epsilon(t)
 =A_L^{(0)}(t)+\delta_\epsilon A_L^{(1)}(t)
   +R_L^\epsilon(t),\qquad
 \|R_L^\epsilon\|=o(\delta_\epsilon).                       \tag{12}
\]

The audit must decide whether \(\delta_\epsilon\) is an
\(\epsilon\), \(\epsilon^2\), \(L^{-1}\), mixed logarithmic, or other source
scale. It must separate:

1. the finite-core Green-kernel correction to the filament Hamiltonian;
2. the source periodic-orbit and period correction;
3. the shape-slaving and vertical-modulation correction;
4. center/first-mode chart connection terms; and
5. self and cross terms between the two physical labels.

If \(\Phi_0(t)\) is the exact filament fundamental matrix, the first return
correction must include both the Duhamel term and the period shift:

\[
 M_L^{(1)}
 =\Phi_0(\mathcal T_0)
   \int_0^{\mathcal T_0}
    \Phi_0(t)^{-1}A_L^{(1)}(t)\Phi_0(t)\,dt
   +\mathcal T_1 A_L^{(0)}(\mathcal T_0)
      \Phi_0(\mathcal T_0).                                 \tag{13}
\]

Every entry in (13) must be traced to an actual GHH coefficient or an exact
Euler variation. No coefficient, discriminant, or hyperbolic example from
the planar annulus may be substituted.

## Common full-contour bridge

A projected matrix \(P_LMP_L\) is not the period map of an invariant Euler
sector. Both competing verdict routes must earn one of the following actual
bridges.

### Bridge C1: exact periodic invariant graph

With \(P_L+P_H=I\), decompose the full generator as

\[
 A_*=
 \begin{pmatrix}A_{LL}&A_{LH}\\A_{HL}&A_{HH}\end{pmatrix}.
                                                               \tag{14}
\]

Construct a periodic bounded graph \(Y(t):\mathcal L_{J_0}\to\mathcal H\)
solving

\[
 \dot Y=A_{HL}+A_{HH}Y-YA_{LL}-YA_{LH}Y,                    \tag{15}
\]

with the area, phase, translation, energy, and impulse conditions respected.
Then

\[
 A_{\rm eff}=A_{LL}+A_{LH}Y                                \tag{16}
\]

is the actual induced finite generator. The graph embedding and its inverse
must be uniformly bounded in the declared Sobolev topology.

### Bridge C2: separated hyperbolic dichotomy or actual Riesz block

If the calculated finite-core block has an off-circle multiplier separated
from the rest of the actual monodromy, construct its exponential dichotomy
or Riesz projection and prove persistence under the measured low--high
coupling. The separation and coupling norm must be actual bounds. The
conditional statement
\(\sigma_{\rm ess}(M_{\rm red})=\mathbb S^1\) from P253/0024 `CT` cannot be
used as an unconditional gap theorem.

An elliptic low block has no analogous automatic bridge in the presence of
dense high phases. A positive low metric is therefore a completed low/center
dependency, after which the high-mode reducibility or metric construction
remains active; it is not full-pair stability by itself.

## Route A: exact reduced positive/Krein metric

After applying (4) and an earned bridge, let \(A_{\rm low,red}(t)\) be the
actual finite reduced generator. Construct a periodic symmetric matrix
\(G_{\rm low}(t)\) such that

\[
 \dot G_{\rm low}
 +A_{\rm low,red}^TG_{\rm low}
 +G_{\rm low}A_{\rm low,red}=0,\qquad
 0<mI\le G_{\rm low}(t)\le MI.                              \tag{17}
\]

Equivalently, its exact return obeys

\[
 M_{\rm low,red}^TG_{\rm low}(0)M_{\rm low,red}
 =G_{\rm low}(0).                                          \tag{18}
\]

Route A must calculate the actual KKS signature and Hamiltonian quadratic
form on every retained block, prove semisimplicity of all unit multipliers,
and give a uniform conditioning bound. The strict-leading sign from
P253/0029 is only the zeroth rung. The final GHH conjugation is not presumed
KKS-symplectic on the unsymmetrized full-label space, so its diagonal metric
cannot be relabeled as (17) without an explicit lift.

**Route-A success:** (17)--(18) for the actual reduced low/center invariant
sector, including finite-core corrections and the complete companion
reduction. This supplies the missing finite-sector input to later
full-pair control; it does not by itself settle infinite high modes.

## Route B: actual hyperbolic or Jordan alternative

For each actual retained real symplectic block, compute its characteristic
polynomial and Krein form. In a two-dimensional reduced block,

\[
 \det M_B=1,\qquad
 \Delta_B=(\operatorname{tr}M_B)^2-4.                       \tag{19}
\]

The alternatives are kept distinct:

| Exact actual result | Route-specific conclusion |
|---|---|
| \(\Delta_B<0\) | elliptic block; construct its positive metric and return to Route A |
| \(\Delta_B>0\) | reciprocal real off-circle multipliers; transfer by C1 or C2 |
| \(\Delta_B=0\), \(M_B\ne\pm I\) | nontrivial parabolic/Jordan block; construct its generalized full-contour Jacobi vector and exact power law |
| \(M_B=\pm I\) | semisimple scalar return; no Jordan claim |

For blocks of dimension greater than two, use the full symplectic
characteristic polynomial and collision-specific Krein form rather than
applying (19) componentwise. A leading filament Jordan shear, an annular
hyperbolic example, real diagonal frequencies, or a sampled multiplier is
not an actual finite-core verdict.

**Route-B success:** a named actual GHH finite-core hyperbolic or Jordan block,
with a C1/C2 full-contour bridge and correct reduction status. If the exact
block is elliptic, Route B continues constructively into Route A rather than
reporting only the absence of growth.

## Frozen source and dependency inventory

No body listed here is opened for new 0033 work before activation.

| Source/input | Frozen provenance | Intended post-activation use | Non-inheritance boundary |
|---|---|---|---|
| García--Hassainia--Hmidi, *Time-periodic leapfrogging vortex rings in the 3D Euler equations* | `arXiv:2603.21644v1`; cached PDF SHA-256 `7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6`; text SHA-256 `f6461cfc68e2aa3b08cadc3409a36f4f0267c5a7c73fdb03e46f9636c121c871` | Theorem 1.1 quantifiers; pp. 21--24 Hamiltonian/contour convention; pp. 24--34 filament Hamiltonian, orbit, period and twist; pp. 60--72 exact self/cross linearization; pp. 87--100 first/normal split; source expansions defining the first finite-core correction | No source theorem supplies an unsymmetrized full-label low monodromy, Krein sign, positive metric, or stability verdict |
| P253/0015 construction | SHA-256 `90eec03cd419ca92e6796e6e4d3839b6c6d2152ffa2b04ef49048eacd3a5bff9` | Actual labeled KKS action, physical normalization, periodic Jacobi equation, phase and translation modes | Action criticality and weak KKS conservation do not imply a positive Hessian |
| P253/0021 construction and source receipt | construction SHA-256 `d1c27273e34903e76126233076fdcb87cbe41f47ffb37bb8c8700b7087c301aa`; receipt SHA-256 `d09dd90b2f3c70ca953464bcca583dc17d306d73405c9fc169b825d9da54e143` | Actual full independent-label IVP, centroid/shape chart, one-period propagator, KKS transport, and exact invariant flag | One-period bounded invertibility gives no multiplier or power control |
| P253/0024 construction and source receipt | construction SHA-256 `8f078bf1d895167bed15813a0a492b92242c287e1f25a66b985e2563e91799ad`; receipt SHA-256 `21397f05f237e2376a77f6af875088c2204092fdd3e000a6c80d6f82dd29548d` | Exact companion power identities and actual-versus-conditional essential-spectrum boundary | `CT` and \(\sigma_{\rm ess}=\mathbb S^1\) remain conditional |
| P253/0029 construction, source receipt, and verdicts | construction SHA-256 `ff0b18f27d76a2778ac70fb26043b72cb32a1159ca42793db3da5b75469f622a`; receipt SHA-256 `4f7eabc9ab86577a6f0fa815366d2e1a53b12b2f2c1946058904ebea07319a34`; verdicts SHA-256 `9ed79edb24559719cb63391068330a00acac65bd79c6147f7ffc4d2a65378874` | Exact label-block representation, shape KKS matrix, strict-leading Hessian sign, principal metric bounds, and high-mode remaining dependency | Its finite-cutoff and principal results do not classify the actual low blocks |
| Hassainia--Hmidi--Roulley, *Invariant KAM tori around annular vortex patches for 2D Euler equations* | `arXiv:2302.01311v1`; PDF SHA-256 `394af598aa77ed3e8646e6a6d28e85c6fd3de9b02725781ccf5045f17665775d` | Method comparator only: pp. 34--40 show how an explicit two-interface discriminant and uniformly conditioned block transform are proved | Nested planar annulus, different interface transports, external parameters, and \(m\)-fold symmetry; none of its coefficients or low-mode verdicts transfer to GHH |
| Accepted `C-FLO-001` | `governance/claims.yaml` | Exact finite-matrix Jordan-sensitive power criterion after the actual block and bridge exist | No Euler realization, domain, KKS form, or full-contour bridge |

The Buttà--Cavallaro--Marchioro concentrated-ring theorem is not a default
input: its published quantifiers concern a singular distinct-radius family
and finite-time concentration, not the GHH unsymmetrized variational period
map. It may be added after activation only if a named bridge estimate is
shown to match every hypothesis.

## Frozen analytic ladder and oracle

After activation, work in this order:

1. audit the listed GHH pages and issue a source-theorem/variable table that
   separates filament, finite-core, source-restricted, and physical full-label
   objects;
2. reconstruct the complete real center/first-mode canonical coordinates,
   individual-area slice, physical KKS matrix, and energy/impulse companions;
3. derive the filament variational return and reduce its twist shear through
   (4);
4. identify \(\delta_\epsilon\), \(A_L^{(1)}\), the period correction, and
   the exact first-return formula (13);
5. build C1 or C2 before treating a finite matrix as an Euler block;
6. execute both outcome routes: positive metric/Krein classification and,
   if exposed, the actual hyperbolic/Jordan alternative; and
7. give exactly one verdict for each tested route, with the strongest positive
   construction and remaining dependency separated.

The strongest oracle is exact Hamiltonian reduction and symbolic/asymptotic
evaluation of the actual source coefficients, followed by an exact
finite-matrix characteristic/Krein calculation and a proved full-contour
bridge. Short exact symbolic checks may verify signs, determinants,
symplecticity, Duhamel terms, quotient dimensions, and Jordan powers.

No production numerical work is licensed by this README. If the analytic
ladder leaves a genuinely irreducible small finite-core splitting,
near-unit multiplier, Krein collision, soft Hessian eigenvalue, or
discriminant sign, 0033 must first load `small-ratio-numerics/SKILL.md` and
append a numerical preregistration specifying:

- the exact continuum matrix entry or scalar remainder and its analytic
  bridge to (17) or (19);
- the dominant filament scale and expected finite-core splitting scale;
- zero-mode/phase/translation gauges and energy/impulse constraint handling;
- precision, quadrature/ODE method, interval or residual error budget, and
  eigenpair residuals;
- observed-order refinement, period and parameter sensitivity, and jitter
  sign tests; and
- a verdict ceiling of the actual bridged finite block only.

Any long run is coordinated with root and executed in scripts pane `w3:p2`.
A sampled finite projection without C1/C2 remains exploratory and cannot be
retrospectively promoted.

## Verdict and continuation contract

Route A and Route B each receive exactly one of `established`, `refuted` with
the actual mechanism, or `blocked` with the missing construction. Bridge C
receives its own verdict and cannot be hidden inside either spectral outcome.

The maximum 0033 positive verdict is a correctly reduced, physically
normalized, actual invariant GHH low/center sector with a uniformly positive
periodic metric and bounded powers. The maximum alternative verdict is a
named actual finite-core hyperbolic or nonsemisimple Jordan block with a
proved full-contour realization. Neither verdict supplies nonlinear orbital
stability, infinite high-mode control, swirl/helicity, intrinsic or quantum
spin, an electron or neutrino identification, parent completion, or a global
no-go.

If the low/center sector is positively controlled, the next achievement is to
join it to the P253/0029 high-mode block construction. If an actual bridged
growing block is found, preserve that route result and continue the parent
robust-carrier search through the other registered same-field carriers; do not
generalize one pair's failure to Euler or to the parent particle objective.
