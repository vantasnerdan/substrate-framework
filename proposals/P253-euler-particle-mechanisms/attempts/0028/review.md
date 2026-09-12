# P253/0028 independent exact-existence review

## Frozen transaction and independence

This review adjudicates root-owned attempt `0027` as one coherent theorem:
existence of a small axisymmetric solitary Euler excitation on a fixed smooth
unbounded pure-swirl column, including its threshold radial mode, same-fluid
exterior, logarithmic axial scale, nonlinear continuation, smooth axis and
finite literal kinetic excess. Stability and particle interpretation are not
part of the unit.

Central registration is present and `0028/activation-schema.exit` is zero. The
frozen `0028` README hash is
`a90ffcb1d0e05ccdbcf6b97db4941af430eab5c457ba04c96ae1094dd848887a`;
the activation exit and stdout hashes are respectively
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`
and `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`.
The accepted base is `v0.183.0`, and the observed review head is
`bd61f489484701c849bc309f4c1df361b85b7f13`. The reviewer authored or
implemented neither the proof nor its API/tests.

The substantive review boundary is:

| Artifact | SHA-256 |
| --- | --- |
| `0027/README.md` | `6cd66e3acb69c39d4e1347fde56cb8f49c406b421091a6fd11b3534a1fe48d41` |
| `0027/source-access.md` | `b6fbb03d64a6017fe123387fb96089ac2eaac825ff5b863eadf610ba1e04ae47` |
| `0027/solitary-wave-construction.md` | `c18d394a757f9fd967975ba23a486e62087c285c8fda4f08d8aeaa0159c6f850` |
| `0027/exterior-construction.md` | `0c2cda195dd5bea9163a7a6c582bfcff746a608a2e48d0a9cd7a007c1b260fc3` |
| `0027/result.yaml` | `50c022db1a6a9f542b10e732f3aa043b1443b2832a18e5ae5ff241ad0df43c29` |
| `0027/validation.md` | `fbfd4c387d1e1a5b39d82b2edf0eb64b2fc6fd97b4559cabb3453992f3ee18cc` |
| `0027/verify_column_structure.py` | `f5c7a0b21a942622d269975232f82cf2e570f6df3ac0a45ac036c76aaa122b3a` |
| `0027/structure-check.stdout` | `147ccf21de0fda94613e965e8a89fabbdff51f703751544c17bb83400de430aa` |
| `src/substrate_framework/euler_column_wave.py` | `1f43a3fdd0539c5bac9305a189796d2e82cf49cfbb303dac6c30183e8557e106` |
| `tests/test_euler_column_wave.py` | `27bdd5bd78d11b1bec1dfac81958e62e5eaeac4e4669f9ad87f8a5924a592816` |
| `0027/column-api.stdout` | `ea67cf02f51ac98ad7fc948e600af9a25a19d653a2e02fbbc631649858aa61da` |

Both captured exits are zero. The structure script reports six exact identities
and the focused API receipt reports four passing tests in 3.64 seconds. Those
receipts are corroboration at this unchanged boundary, not substitutes for the
operator convergence and implicit-function proof audited below.

## Strongest supported positive statement

Let `w(r^2)` be smooth and nonnegative, constant and strictly positive near
`r=0`, zero for `r>=R_0`, and therefore nontrivially transitioning between
those regions. Put

    L(r)=integral_0^r s w(s^2) ds,
    U=(L(r)/r)e_theta.

Then `U` is a smooth stationary constant-density Euler column on all of
`R^3`. It has compact transverse axial vorticity and the physical circulation
tail `L(infinity)/r`; hence its total kinetic energy is infinite. There is a
positive critical speed `c_0` determined by the first weighted radial mode.
For every `c>c_0` sufficiently close to `c_0`, there is a nonzero smooth
axisymmetric traveling solution

    u(t,r,z)=U(r)+v_c(r,z-c t)

of the full Euler equations. The perturbation belongs to `H^m(R^3)` for every
finite `m`, tends to zero relative to the column, has an exact irrotational
same-fluid exterior rather than a wall, and has absolutely convergent literal
kinetic excess

    (rho_m/2) integral_R3 (|U+v_c|^2-|U|^2) dx.

The branch is one-sided and small: if `mu` is the first interior radial
eigenvalue, its streamfunction amplitude is `O(mu)` and its axial length is
`L_mu`, where

    mu L_mu^2=f_0(R)^2 log L_mu

on the unique large-root branch. This is an exact existence theorem for one
traveling profile at each small parameter, not stability of nearby data.

## Source boundary

The inaccessible Sun 1995 body is not used. Its publisher extract supports
only provenance for a rigid-cylinder solitary-wave theorem and cannot supply
the present hypotheses, operator domain, or unbounded exterior. The cached Sun
1994 source, hash
`a06730b4f84611bd53390b5450d12684ad39a7715a23e93dd9f4e4dcf3ec50cc`,
is correctly classified as weakly nonlinear/interface motivation: its free
surface or two-density/cylinder assumptions do not prove this same-fluid
whole-space theorem. The load-bearing proof in `0027` is direct.

## Full Euler reduction and background

In the frame moving at speed `c`, use

    u_r=-psi_z/r,  u_z=psi_r/r,  u_theta=F(psi)/r,
    p=B(psi)-(|grad psi|^2+F(psi)^2)/(2r^2).             (R1)

Direct substitution into cylindrical Euler gives

    Delta_* psi+F F'(psi)-r^2 B'(psi)=0,
    Delta_*=partial_rr-r^-1 partial_r+partial_zz.        (R2)

The focused Cartesian/cylindrical residual test independently confirms that
the Euler residual is `-(R2) grad psi/r^2`, so no momentum row or pressure term
is omitted. Density is handled consistently by treating `p` as specific
pressure; the physical kinetic statements restore `rho_m`.

For an axial column `(W,V)` with `q=W-c`, `psi_0'=rq` and `L=rV`, differentiating
the two material labels gives

    F(psi_0)=L,
    B'(psi_0)=W'/r+L L'/(r^3 q),
    Q_c=2LL'/(r^3q^2)-(W''-W'/r)/q.                    (R3)

For the selected pure-swirl background, `W=0`, `q=-c`, and
`Q_c=Phi/c^2`, `Phi=(L^2)'/r^3>=0`. Smoothness of `w(r^2)`, its constant core,
and flat compact support make the labels smooth: near the axis `F` is linear
and `B'` constant, while beyond the vorticity support both labels are constant.
Thus the nonlinear remainder vanishes near the axis and in the exterior and is
supported in a fixed transition annulus. Adding `c e_z` after reconstruction
returns to the laboratory frame and cancels the translating background
velocity with the correct sign.

## Critical radial mode and exact exterior

On `H=L^2((0,R),dr/r)` let `V` be the completion under
`integral f_r^2 dr/r`, with `f(0)=0`. The weighted Hardy estimate
`|f(r)|<=r||f_r||_H/sqrt(2)` controls the axis, makes the trace at `R`
continuous, and gives the compact embedding needed by the variational problem.
The critical value

    c_0^2=sup_f
      [integral Phi f^2 dr/r]/[integral f_r^2 dr/r]       (R4)

is attained. Sturm theory gives a positive simple maximizer `f_0`; its natural
boundary condition is `f_0'(R)=0`. The identity
`(f_0'/r)'=-Phi f_0/(c_0^2 r)` makes `f_0` nondecreasing and nonconstant in the
core, hence `h_0=f_0(R)>0`. For `c>c_0`, the first eigenvalue
`mu=lambda(c)` is positive and

    lambda'(c_0)=2 integral Phi f_0^2 dr/(c_0^3 r)>0.      (R5)

The remaining finite-interval eigenvalues retain a positive gap. This is not
misstated as a whole-space `L^2` eigenmode: its constant zero-frequency
exterior is a threshold resonance.

For axial Fourier number `k!=0`, the unique decaying same-fluid exterior is

    fhat(k,r)=g_hat(k) r K_1(|k|r)/(R K_1(|k|R)),          (R6)

with outward derivative `-|k|K_0/K_1` at `R`. Its positive minimized energy
trace is

    T_R(k)=|k|K_0(|k|R)/(R K_1(|k|R)),
    T_R(k)=k^2[log(2/(|k|R))-EulerGamma]+o(k^2).           (R7)

The sign and `1/R` factor follow both from differentiating (R6) and from the
boundary flux of `integral(f_r^2+k^2f^2)dr/r`. Matching `f` and `f_r` matches
both poloidal velocity components; constant exterior `F,B` and (R1) then match
swirl and pressure. The Rankine comparison correctly distinguishes this
Neumann threshold (`J_0=0`) from a rigid-wall Dirichlet mode (`J_1=0`).

Consequently the full interior-plus-exterior form operator is

    L_c(k)=A_c+k^2 I+T_R(k)b tensor b,                    (R8)

where `b(f)=f(R)`. It is positive with lower bound `mu+k^2`. The stated
rank-one Schur inverse is the exact Sherman--Morrison/Feshbach inverse on
`V -> V*`; the trace is never treated as an `L^2` delta. Positivity of the
radial complement gives `d_k>=0`, so its denominator cannot vanish.

## Focusing coefficient and logarithmic operator limit

Twice differentiating the actual labels at fixed `r` gives

    J_c=2Q_c'/(r q).

At `c=c_0`, direct integration by parts using the radial equation yields

    integral Q_c' f_0^3 dr/r^2
      =-(1/2)integral [f_0^3/r^5]
          h(3h^2-8h+8)dr <0,
    h=r f_0'/f_0.                                         (R9)

All endpoint terms vanish: `f_0=O(r^2)` at the axis, and `Q_c=0`,
`f_0'=0` at the matching radius. Since
`3h^2-8h+8=3(h-4/3)^2+8/3>0`, `h>=0` and `q=-c_0<0`, the projected coefficient

    beta=(1/2)integral J_c f_0^3 dr/r

is strictly positive. The structure verifier derives this identity rather
than assuming a sign.

The functional-analytic bridge does not follow from those exact checks; it is
supplied by the proof. With `mu L_mu^2=h_0^2 log L_mu`, set

    K_mu(kappa)=mu L_c(kappa/L_mu)^-1.                    (R10)

The Schur denominator divided by `mu` converges on bounded `kappa` to
`1+kappa^2`; the complement vanishes and the dressed critical vector tends to
`f_0`. This gives `K_mu -> P_0/(1+kappa^2)` locally. The proof also controls the
entire Fourier line:

- the coercive resolvent estimate makes both `K_mu` and its derivative vanish
  uniformly for `|kappa|>=L_mu^delta`;
- for `M<=|kappa|<=L_mu^delta`, small-`k` Bessel bounds keep the complement
  uniform and bound the critical part by `C/(1+kappa^2)` and its derivative by
  `C/(|kappa|(1+kappa^2))`;
- bounded `kappa` uses the differentiated Bessel expansion, continuously at
  zero.

The exact exterior minimization identity

    T_R'(k)=2k integral_R^infinity |f_k|^2 dr/r,
    0<=kT_R'(k)<=2T_R(k)                                  (R11)

supplies the derivative control rather than an assumed differentiability of a
threshold eigenvalue. Taking first `M` large and then `mu` small proves uniform
`H -> V` convergence of both the multiplier and its first `kappa` derivative.
Thus the logarithm and exterior tail are present in the actual operator limit.

## Weighted nonlinear map and implicit-function theorem

On the real even spaces

    Z^s={F, XF in H^s(R_X;V)},
    Y^s={G, XG in H^s(R_X;H)},  s>1/2,                    (R12)

the first multiplier derivative is exactly what controls the spatial weight:
Fourier transformation of `XK_mu G` produces `iK_mu' Ghat+K_mu Fourier(XG)`.
Hence (R10)--(R11) give operator-norm convergence `Y^s -> Z^s`, not merely
pointwise symbol convergence.

Because the nonlinear remainder is smooth, quadratic, and supported in the
fixed radial transition annulus, the radial `V` embedding there and the
one-dimensional `H^s` algebra estimate give

    n_mu(F)=mu^-2 N_c(r,mu F)
      ->(J_c0/2)F^2                                      (R13)

in `C^1` on bounded subsets `Z^s -> Y^s`, including the single `X` weight.
Thus the rescaled full equation is exactly

    F-K_mu n_mu(F)=0.                                    (R14)

At `mu=0`, its nonzero solution is

    F_*=f_0 A_*,
    A_*(X)=3 sech^2(X/2)/(2 beta),
    A_*-A_*''=beta A_*^2.                                (R15)

The derivative is the identity on the radial complement and is triangular
over the critical component. The scalar Schrodinger representative is
`-partial_XX+1-3sech^2(X/2)`. Its only zero mode is the odd translation
`A_*'`; on the frozen even space zero is absent, while the negative even mode
does not obstruct invertibility. The essential spectrum starts at one and the
remaining discrete spectrum stays away from zero. Its inverse preserves one
spatial weight because the restricted zero-energy Green kernel decays
exponentially. Therefore the full derivative is an isomorphism on (R12).

The `C^1` operator convergence, not the algebraic verifier, now invokes the
Banach-space perturbative implicit-function theorem and produces

    F_mu -> f_0 A_* in Z^s,
    f_mu(r,z)=mu F_mu(r,z/L_mu)                            (R16)

for every sufficiently small `mu>0`. Since `lambda'(c_0)>0`, this is precisely
a one-sided interval of speeds `c>c_0`. Nontriviality follows from the nonzero
limit. Evenness fixes the translation phase; smooth dependence at `mu=0` is
not required.

## Axis, whole-space reconstruction, and energy

Writing `f=r^2 phi` conjugates the singular radial part exactly to
`phi_rr+3phi_r/r`, the radial four-dimensional Laplacian inside a
five-dimensional elliptic equation including `z`. The `V/H` weights become
the local radial `H^1/L^2` weights for `phi`. Because coefficients are smooth
in `r^2` and the nonlinear forcing vanishes near the axis, elliptic regularity
gives a smooth even Cartesian extension. In particular `f/r^2` and `f_r/r`
are `O(mu)`, so every reconstructed velocity component and the pressure are
smooth and the no-stagnation/label-range condition persists.

At `R`, the exact exterior DtN relation gives matching first derivatives.
Both sides solve the same homogeneous equation in a radial margin where the
labels are constant, so elliptic transmission bootstraps the match to
smoothness; `R` is only an integration surface. The Fourier exterior energy is
exactly (R7), including its low-frequency tail.

The weighted space has `F in L^1_X(V)` by Cauchy--Schwarz with
`(1+X^2)^-1`. Thus the radially supported swirl perturbation is absolutely
integrable in `z`, making the only background cross term finite. The poloidal
perturbation is orthogonal to the swirl background, its interior kinetic norm
is the `V` form, and its exterior kinetic norm is (R7). With the cylindrical
measure `2 pi r dr dz`, both `integral |v|^2` and
`integral |2U dot v|` converge absolutely. This proves the literal excess
claimed in the result without subtracting two undefined infinite energies.
The background itself remains an explicitly infinite-total-energy column.

Elliptic bootstrapping at arbitrarily high finite regularity gives
`v in H^m(R^3)` for every finite `m`, while the one spatial weight and Sobolev
embedding give decay to the background. Replacing `z` by `z-ct` and adding the
Galilean `c e_z` reconstructs the advertised global-in-time traveling Euler
solution.

## API and oracle audit

The importable API exposes only the generic column coefficients and exact
exterior Bessel solution. Its documentation explicitly says it is not a
numerical wave constructor or certification of the IFT proof. The tests derive
the full Euler residual through the existing cylindrical field API, generic
`Q,J`, exterior equation/DtN/energy flux, the finite-dimensional algebra of the
rank-one inverse, and the homoclinic/odd translation mode. They are sensitive
to the principal signs, missing exterior factor, rank-one denominator, and
homoclinic coefficient. They do not claim operator convergence.

The analytic proof above is the strongest oracle. The exact script and four
tests are correctly classified as exposing regression/corroboration. Their
sources and inputs are unchanged, so the captured first executions were
inspected rather than rerun. No sampled eigenvalue, fitted profile, soft
stability edge, or production numerical observable is used; small-ratio
numerical prescriptions do not bind.

## Findings and verdict

| Finding | Direct evidence | Disposition |
| --- | --- | --- |
| The whole-space critical object is a threshold resonance, not an isolated radial `L^2` eigenfunction. | Constant `k=0` exterior and logarithmic divergence of the auxiliary full radial `H` norm. | Correctly scoped and explicitly handled by eliminating the exterior before using the finite-interval gap. |
| Pointwise exact algebra cannot prove the multiplier limit or IFT. | The structure/API checks contain no infinite-dimensional norm or nonlinear-map argument. | Correctly treated as corroboration; the Schur estimates, three-region Fourier proof, weighted `C^1` convergence and even-sector isomorphism supply the missing analytic bridge. |
| Sun 1995 is inaccessible and Sun 1994 has different physical boundaries. | `source-access.md` and the pinned Sun 1994 audit. | Neither is imported; the direct proof is self-contained. |

No concrete correction is required. The route verdict is **established as
stated** at `evidence_scope:
DIRECT_OPERATOR_AND_NONLINEAR_IMPLICIT_FUNCTION_EXISTENCE_PROOF`.

- Verification: `exact_analytic_existence_verified`; executable algebra is
  corroborating evidence.
- Review: `established_as_stated`.
- Compatibility: `compatible_proposal_evidence`; the background's infinite
  total energy and the perturbation's finite literal excess are kept distinct.
- Epistemic: exact for the declared smooth axisymmetric one-sided branch and
  its same-fluid whole-space exterior.
- Relationship: this establishes an exact carrier/excitation supplier only.

The next scientific dependency lies outside this review: a restoring or
interaction theorem for nearby states with the full ambient pressure/action
must be proved separately. This verdict does not review or imply spectral,
orbital, asymptotic, or nonlinear stability; an all-time open neighborhood;
particle identity; physical or quantum spin/statistics; electron/neutrino
completion; parent completion; or a global Euler no-go.

## Correction check

No correction was requested. One bounded correction check remains available
only if the frozen theorem statement or a directly supporting edge changes.
