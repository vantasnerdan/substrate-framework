# P253/0035 independent review of the 0030 column action and propagation claim

## Frozen transaction

This is a non-author/non-implementer review of root-owned `0030` at commit
`1d378f9aebd10a8ff0dd91f26bf15d371f3e243e`. Central activation is valid:
`0035/activation-schema.exit` is zero. The frozen README hash is
`f1d6d5036585219cdfb058f365a8ee52553b0ba66bb59ad7a6ccd0441653269b`.
The activation command, exit, stdout, and empty stderr hashes are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
and `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.

The reviewed boundary is:

| Artifact | SHA-256 |
|---|---|
| `0030/README.md` | `1728735f7d6aa91179f3b264ee7b1b52c1430c186b9df03efdc80fabbb35180e` |
| `0030/source-inventory.md` | `d47f5faf690c99f4e0703e184f4e015e135479a094a734d04e880ba2aeceff15` |
| `0030/action-construction.md` before correction | `56ee231804c1ccfbd60b44546a7943f6518731f52bdd7a8ab95f777ce1aa1fe6` |
| `0030/action-construction.md` after correction | `ad59399682c9de070a4b0102b0b893d4a2408b9dccb996eaf6f9228fa197c24e` |
| `0030/review-correction-0035.md` | `537b54751a7c98257267478e034bbe1cba360e69a28165c077be6fbe7a25add8` |
| `0030/column-propagation.md` | `d50facb41126c1b0f31609dee04d733dfff720d34df4e83b59c5ca707f69de2e` |
| `0030/column-profile-transfer.md` | `c4c9acc3ec29e709a0124a402e85ad688055cc945365b1788ceea23a95bd5fba` |
| `0030/result.yaml` | `43b741c6550807a0c1ae97307daa75c6abb7e849d1b21587e9162cc128debac8` |
| `0030/validation.md` | `c13fb4c529e456fd48210b55d47d4a2380cc990b6ebb99a5ed1cf2682976e305` |
| `src/substrate_framework/euler_column_wave.py` at `1d378f9` | `b66328cfafe82fb1f56eb2b52911786ddb5686718befc10ad34dfb0dc8c11dc0` |
| `tests/test_euler_column_wave.py` at `1d378f9` | `d106dc0af6b1db2cb9575df7ea4ba5b856583280da66fa0f6358c35b76ca7974` |
| `0030/column-action-api.stdout` | `ff374897f35fdad44fb3cfd82d5aa0733a061d861049fafcb35d74c750d66ff8` |
| `0030/checkpoint-validation.head` | `26b4b082e5ecd03ba894394d51bbd68614b4befa40f21e3eb3376c1c2767c816` |
| `0030/checkpoint-validation.stdout` | `0d2408c1c2fb1896f3a8ebf66a46c197ba73a83d272d4e11612c783b35e45e37` |

Both receipt exit files contain zero and have SHA-256
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`;
both stderr files are empty. The focused receipt records seven tests passing
in 4.03 seconds. The checkpoint head is exactly `1d378f9` and records 2712
tests passing in 455.12 seconds. These were reused and not rerun. The tests
expose the cylindrical Euler residual, both linearized rows, pressure/Hodge
reconstruction through `K`, weighted energy flux, the Bessel exterior, and a
centrifugal-factor mutation; they do not prove the functional domain,
unitarity closure, speed differentiation, or profile-class transfer.

The primary bodies are Gallay--Smets,
[*Spectral stability of inviscid columnar vortices*](https://arxiv.org/abs/1805.05064v3),
`arXiv:1805.05064v3`, PDF SHA-256
`081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`,
and [*On the linear stability of vortex columns in the energy
space*](https://arxiv.org/abs/1811.07584v2), `arXiv:1811.07584v2`, PDF
SHA-256
`a895efa1acf6066c0b50881fe08002be5d0a73f92f413ec8c02bebcf31b370d8`.
Accepted `0027/0028` supplies only the unchanged exact background/solitary
existence and exterior estimates. Attempts `0031`, `0034`, `0032`, and `0033`,
including all later source edits, were excluded.

## Strongest supported result

After the bounded action-domain wording correction checked below, the coherent
`0030` claim is established at its stated restricted scope.

For the explicit smooth pure-swirl Euler column constructed in `0030`, whose
axial vorticity is constant near the axis, strictly decreasing on one annulus,
and flat compactly supported, the following hold:

1. The translating solitary field has the displayed physical
   energy--impulse--Casimir first and second variations on every regular
   invertible-label patch and on the corresponding compact-generator
   dynamically accessible tangents. The full-leaf traveling Hessian has
   accessible high-axial-frequency directions of both signs, so it cannot be
   a uniformly positive kinetic metric. This is a route obstruction, not a
   linear or nonlinear instability theorem.
2. The stationary column has a distinct positive energy--Casimir quadratic
   form. Its full-pressure axisymmetric linearization generates an exact
   unitary group on the weighted accessible column-energy completion. The
   result retains the divergent flat-edge swirl weight and is not an
   unrestricted unweighted `L2` or nonaxisymmetric theorem.
3. Every nonzero axisymmetric radial wave branch has phase speed strictly
   below `c0`; its group speed has magnitude no larger than its phase speed
   and hence no larger than `c0`. The accepted solitary branch has `c>c0`, so
   the exact Doppler separation is `(c-c0)|k|`, with no uniform gap at `k=0`.
4. The explicit compact-vorticity profile lies in the global `C_b^1` closure
   of the strict Gallay--Smets profile class. Therefore Remark 1.9 of
   `1805.05064v3` gives imaginary spectrum for each fixed integer `m` and each
   fixed `k!=0` in the solenoidal enstrophy space `X_{m,k}`. This does not
   transfer the all-mode velocity-space group bounds of `1811.07584v2`.

The nonlinear background inequality is also correct at its stated conditional
scope: for classical axisymmetric evolution with `0<=xi<=L_infinity`, finite
relative energy/Casimir integrals, and vanishing fluxes, the conserved convex
relative functional bounds the physical kinetic perturbation norm. It is not
global regularity with swirl and is not stability of the nonzero solitary
wave.

## 1. Physical variation and accessibility

With `dnu=r dr dz`, `xi=r u_theta`, `zeta=omega_theta/r`, and
`psi_lab=K zeta`, direct variation gives

    delta E/delta zeta=psi_lab,  delta E/delta xi=xi/r^2,
    delta I/delta zeta=r^2/2.

Thus for `psi_frame=psi_lab-c r^2/2=Psi(xi)` the functional

    A=E-cI-integral zeta Psi(xi)dnu
      -integral B(Psi(xi))dnu

has exactly the two first-variation rows printed in `0030`. Substituting
`zeta=F F'/r^2-B'` and `Psi'=1/F'` cancels the `xi` row on every regular label
patch. A second variation gives

    Q=integral[eta K eta-2Psi' eta chi
      +(r^-2-zeta Psi''-(B composed Psi)'')chi^2]dnu,

with one common physical multiplier `2*pi*rho_m`; no density or cylindrical
measure is missing.

The axisymmetric coadjoint tangent

    chi={xi,g},  eta={zeta,g}+{xi,a}

is the correct dynamically accessible class. On a compact subpatch where
`xi_r!=0` and `Psi'!=0`, high-`z` generators independently realize the leading
amplitudes of `chi` and `eta`. Taking axial frequency `N`, `eta` of size
`N chi`, and using that `K` has order minus two leaves the positive kinetic
terms of order one while the mixed term is of either sign and order `N`.
This validates the claimed noncoercivity on actual accessible directions.

For the background `L=r^2 Omega`, `C'(L)=-Omega` and

    r^-2+C''(L)=2L/(r^3 L')>0.

The condition that the axial vorticity is nonincreasing gives
`2L/(rL')>=1`, so the quadratic form controls physical kinetic energy. At the
flat edge `L'->0`, the coefficient diverges; retaining that weighted domain is
essential and is done explicitly. The exterior sign in the nonlinear convex
inequality is also correct because `r>=R0` and `chi=xi-L_infinity<=0` make
`L_infinity(1/r^2-1/R0^2)chi>=0`.

Finally, on the pure-swirl background,
`Psi(L(r))=-c r^2/2` and `(B composed Psi)'=Omega`. Hence

    A=H_C-cP_C,
    P_C=I-(1/2)integral zeta r(xi)^2 dnu,

where the second term is a Casimir. Expanding at the column gives the stated
leading term `-integral (r/L') eta chi dnu`. The sign and physical role as the
same translation generator on the leaf are correct.

### Bounded correction and check

Equation (4) does not vanish pointwise as an unrestricted field derivative in
the constant-label exterior: there `xi=L_infinity` and `zeta=0`, while
`psi_frame` varies, so no single-valued inverse value `Psi(L_infinity)` can
make `A_zeta=0` throughout that region. The result remains true on its actual
accessible domain. For a compact smooth generator, flatness gives
`chi={xi,g}=0` and `eta={zeta,g}+{xi,a}=0` in the exterior and removes an
interface distribution, so `delta A` pairs to zero even though the coefficient
`A_zeta` need not vanish there. Endpoint Casimirs are therefore orbitwise or
extension-dependent, and `Q` plus the high-frequency construction belongs
strictly inside the regular label set.

The one authorized correction changes the source hash from
`56ee231804c1ccfbd60b44546a7943f6518731f52bdd7a8ab95f777ce1aa1fe6`
to `ad59399682c9de070a4b0102b0b893d4a2408b9dccb996eaf6f9228fa197c24e`;
its append-only receipt hash is
`537b54751a7c98257267478e034bbe1cba360e69a28165c077be6fbe7a25add8`.
The diff adds exactly the exterior accessible-tangent identity, absence of an
interface distribution by flatness, orbitwise/declared endpoint treatment,
and regular-interior support for the Hessian argument. It removes the
unrestricted full-field wording and adds no broader claim. This fully closes
the finding. No displayed pre-existing equation, API, test, source-transfer,
or other oracle changed, so the unchanged receipts were correctly not rerun.

## 2. Positive weighted energy and exact unitary propagation

The two full-pressure axisymmetric linearized equations are

    eta_t=(2L/r^4) chi_z,
    chi_t=(L'/r)(K eta)_z.

They follow directly from cylindrical Euler after eliminating pressure by the
poloidal Hodge reconstruction `psi=K eta`; the focused test independently
reconstructs both rows from the full radial and axial momentum equations.
With `a=2L/(r^3L')`, self-adjointness of `K` and integration in `z` make

    (1/2)integral(eta K eta+a chi^2)dnu

exactly conserved.

For each `k!=0`, let `p=K_k^(1/2)eta`, `q=a^(1/2)chi`, and

    D_k=K_k^(1/2)(2L/r^4)a^(-1/2).

The identity `a^(1/2)L'/r=(2L/r^4)a^(-1/2)` gives the block generator
`ik[[0,D_k],[D_k^*,0]]`. Its squared singular values are the generalized
radial eigenvalues

    c_n(k)^2=B[f]/A_k[f],

and `||D_k||=c_1(k)<=c0`. Thus each block is bounded self-adjoint after
removing the factor `i`, and the measurable direct integral of
`k[[0,D_k],[D_k^*,0]]` is self-adjoint on its natural graph domain. Stone's
theorem gives the exact two-sided unitary group. The measure-zero `k=0` fiber
is stationary; no isolated radial `L2` mode is invented there.

This construction includes the nonlocal exterior pressure in `K_k`. It also
explains the precise topology: smooth accessible data are completed in the
positive weighted form, exterior poloidal energy remains, and independent
exterior swirl perturbations are absent because `L'=0`. The theorem is not
unitarity on all axisymmetric physical-energy perturbations and does not rely
on the weaker Gallay--Smets energy-space estimate.

## 3. Exact phase and group-speed ceiling

The exterior minimizer at fixed trace is

    f_k(r)=r K_1(k r)/(R K_1(kR)),

and its energy trace is

    T_R(k)=k K_0(kR)/(R K_1(kR)),  k>0.

The full form is

    A_k[f]=integral_0^R(f'^2+k^2 f^2)dr/r
           +T_R(k)|f(R)|^2.

Since `A_k[f]>A_0[f]` for every nonzero mode and the positive compact
generalized problem attains its eigenvalues,

    0<c_n(k)<=c_1(k)<c0,
    c0^2=sup B[f]/A_0[f].

The group bound is independently earned. The exterior envelope identity gives

    T_R'(k)=2k integral_R^infinity f_k(r)^2 dr/r,
    0<=kT_R'(k)<=2T_R(k).

Therefore `0<=k A_k'[f]<=2A_k[f]`. Differentiating a simple generalized
eigenpair at its stationary Rayleigh quotient yields

    d(k c_n(k))/dk
      =c_n(k)[1-k A_k'[f_n]/(2A_k[f_n])],

which lies in `[0,c_n(k)]` on the positive temporal branch; the negative
branch has the opposite group velocity. Reflection handles negative `k`.
Thus both speed conclusions, including the sign-separated temporal branches,
are exact. The Doppler estimate loses one factor `|k|` at the threshold, so it
does not supply a uniform inverse or solitary stability.

## 4. Compact profile and Gallay--Smets closure

On `a<r<b`, set

    J=exp(1/(r-a)-1/(b-r)),
    alpha=4/(1+sqrt(1+4J/r^2)),
    Omega=Omega0 exp(-integral_0^r alpha(s)ds/s),
    W=Omega(2-alpha).

Extending `alpha` by zero inside and two outside is smooth because the endpoint
errors are flat. The resulting `W` is constant near the axis, strictly
decreasing on the transition annulus, flat zero for `r>=b`, and smooth as a
Cartesian radial vorticity. Direct algebra gives

    Phi/(Omega')^2=2r^2(2-alpha)/alpha^2=J.

For `epsilon>0`, the displayed

    J_epsilon=J/(1+epsilon r^2J)+epsilon/r^2

has `J_epsilon'<0` and `(J_epsilon/r^2)'<0`; the corresponding
`alpha_epsilon` is strictly increasing and
`W_epsilon'=Omega_epsilon[-alpha_epsilon(2-alpha_epsilon)/r
-alpha_epsilon']<0`. Near zero it has
`alpha_epsilon=O(sqrt(epsilon)r^2)` and `W_epsilon'(0)=0`; beyond `b`,
`W_epsilon=O(epsilon r^-6)`. Hence circulation is finite,
`r^3W_epsilon'->0`, and `rJ_epsilon'->0`.

The flat-coordinate argument `q=J^-1/2` at `a`, ordinary flatness at `b`,
the explicit inner bound, and the integrable exterior tail prove

    ||W_epsilon-W||_infinity
      +||W_epsilon'-W'||_infinity ->0.

After `Omega0=1`, every approximant has the source normalization `W(0)=2`.
This is exactly global `C_b^1` convergence, not local ODE convergence.

Remark 1.9 of `1805.05064v3` states that
`W -> L_{m,k}^W` is continuous from `C_b^1` into bounded operators on the
fixed enstrophy space and explicitly extends Theorem 1.3 to the closure,
including profiles that may be compactly supported. It therefore supplies the
claimed fixed-`(m,k)`, `k!=0`, imaginary-spectrum conclusion directly; no new
Evans-function theorem is being assumed by `0030`.

By contrast, Theorem 1.1 of `1811.07584v2` assumes a strictly positive,
strictly decreasing `C2` profile and proves only
`||e^{tL}||<=C_epsilon e^{epsilon|t|}` in the whole-space velocity energy
space. It supplies neither a closure statement for this compact profile nor
exact unitarity, uniform boundedness, polynomial growth, or nonlinear
stability. `0030` keeps this distinction correctly.

## Verdict, evidence roles, and next dependency

The direct calculus/operator/source audit establishes the joined fixed claim
after the one bounded accessible-domain wording repair. Verification is exact
analytic; review is established at the restricted domains above;
compatibility is proposal evidence on the accepted `0027/0028` background;
epistemic status is exact for the column and fixed-mode statements and
explicitly open for solitary propagation/stability.

The next scientific construction is the actual nonzero solitary-wave Leray
generator and its zero-frequency translation/speed modulation, with a
continuum-complement estimate retaining the Bessel pressure kernel. Nothing
here establishes nonlinear or orbital stability of the solitary excitation,
an all-mode three-dimensional unitary group, global regularity for arbitrary
axisymmetric swirl, a particle, physical/quantum spin, electron/neutrino
completion, or parent-campaign completion.
