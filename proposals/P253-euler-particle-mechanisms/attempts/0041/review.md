# P253/0041 independent review of the exact 0037 threshold subclaims

## Frozen transaction and independence

This transaction independently reviews five exact or conditional subclaim
groups from root-owned attempt `0037`. The reviewer authored and implemented
none of the target work and opened no target scientific body before central
activation. The frozen review README has SHA-256
`5e9e77e2e79dd355e3b5bef3546dfb7bf06c5b8c8c766764e16e17b7725a40b4`.
Central registration names `particle-balance-review`; activation contains
exactly `0`, with SHA-256
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`,
activation stdout
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
and empty stderr
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
The accepted authority base is release `v0.183.0`; the observed
preregistration head was
`ca929ce584b462c033fde7866e6231076445ecd5`.

The substantive-pass boundary was the content-blind inventory and hashes in
the frozen README. Four bodies received the one authorized bounded correction:

| Artifact | Substantive-pass SHA-256 | Corrected SHA-256 |
|---|---|---|
| `0037/threshold-reduction.md` | `0cfe43b2138ce910940189c0ecc56648aca6625ae841baf6150ffe6f93edf900` | `5a2b619df9ec69f3f9865a6efcad11c9090d6b3ce0bab543051fa1c07df13a23` |
| `0037/ray-exclusion.md` | `f7fa76057492c242e2e211aff6fd7036654d0af7bbb8a9d544e303a559f8fc56` | `fe3a0aa0b2e1c2070ad6848b9d5cbe05de49e188e4f71919271e04254819fc52` |
| `0037/scattering-block-reduction.md` | `994e160d8c785f338b0502e1a8f63dc32c597b2259fc8e6e44402499ec0e73ba` | `35fbb57dc23435fcf7158249c4041a86ca5173d7227c365f796b8120cbb9f5ac` |
| `0037/block-scattering-lemma.md` | `e68ff5c3213c412e82fc78afc5861481c709c2771a76a5ee262e1890f7c00195` | `ccc65503594b67ad88c0aa9641f67995f3288b48bd5b9d5ac4c4b050e7100d90` |

The append-only correction receipt is
`0037/review-correction-0041.md`, SHA-256
`78e5846f4726bde8cfa0d6a496e4813e04a2a7f6ab157e65c571abad2202d429`.
The remaining reviewed bodies stayed at their preregistered hashes:

| Artifact | SHA-256 |
|---|---|
| `0037/source-transfer.md` | `705bfefc03b7d0905660e1538aef175d79cd3a94c1056117367fe987794389ff` |
| `0037/mixed-casimir-flux.md` | `28ab32fa54abcb58afb96162a4cd2d2e502b66d967952e9d8e4902db049b36c6` |
| `0037/supervisor-projection-analysis.md` | `2158c90dfad0a0fc6a2a395945a9a5ed1adde787625816cce0164ccffd917728` |
| `0037/verify_threshold_reduction.py` | `7d69f44f6d64b7f39d1b56a20812e2acfa572aec52790aa07f7a23a4c89b8acd` |
| `0037/verify_mixed_casimir_flux.py` | `7d3bcaf5a666f54a97cd6f3b337bf4070417499713bea41f54e4048a12caf930` |

Files added to `0037` after preregistration, including an adjoint-transfer
body and verifier, are outside this transaction and were not opened. Active
`0040` was likewise not opened. No target source, oracle, or canonical API was
changed by the correction, so the unchanged exact scripts were correctly not
rerun.

## Strongest supported positive result

After the bounded correction, all five frozen groups have meaningful exact
content:

1. The actual same-fluid Bessel exterior and critical radial slope determine
   the right-moving threshold dispersion and logarithmic axial scale. The
   exact branch then fixes the unique compatible local conservative KdV
   normal form. Its profile, Hessian spectrum, invariant formulas, and
   constraint slopes are exact scalar results. A time-dependent full-Euler
   adjoint projection and long-time remainder theorem are not supplied.
2. Every free column branch has the stated exact local-in-space, integrated-
   time Plancherel estimate. Once the critical first radial branch and the
   stationary kernel are separated, the remaining radial branches have a
   uniform noncritical group-speed gap.
3. Every material trajectory of each sufficiently small `0027` solitary
   member moves strictly in one axial direction and meets the localized
   excitation only once. The actual weighted exterior makes the added strain
   integrable, and the principal Euler bicharacteristic-amplitude cocycle has
   zero exponential exponent on each ray. This is not a full-operator
   spectral or semigroup theorem.
4. Every smooth mixed label `D` has the displayed exact local axial
   conservation law on the declared axis/exterior domain. The `D=1`
   coordinate is unbounded even on compact-generator accessible graph-core
   data; a compact regular-label `D` supplies a bounded low-frequency
   critical coordinate with positive overlap. Conversion to the actual
   adjoint oscillator pair and the scaled cancellation remains open.
5. The two-channel Kato-supersmooth argument is a correct conditional abstract
   lemma once its generator realization, two-sided inhomogeneous estimates,
   and small Birman--Schwinger norm are assumed. Neither the P-channel
   supersmooth estimate nor the full Euler factorization is proved by this
   lemma.

These are strong threshold, transport, and functional-analytic building
blocks. They do not close the all-time axisymmetric Euler propagation target
in `0037/README.md` and do not support nonlinear P2 or a particle/quantum
claim.

## Source and dependency applicability

Accepted `0027/0028` supplies the exact smooth solitary family, positive
critical finite-interval radial mode, same-fluid modified-Bessel exterior,
logarithmic scale, one-spatial-weight reconstruction, and smallness needed for
strict axial drift. The review rechecked every new use rather than inheriting
stability or propagation. In particular, the exterior symbol is

    T_R(k)=|k| K_0(|k|R)/(R K_1(|k|R))
          =k^2[log(2/(|k|R))-EulerGamma]+o(k^2),

and the accepted weighted reconstruction is enough for the axial `L1` strain
bound used by the ray calculation.

Reviewed `0030/0035` supplies the positive weighted column oscillator and the
exact branch bounds

    0<c_n(k)<=c_1(k)<c_0,
    0<=d(k c_n(k))/dk<=c_n(k).

It does not supply a solitary-wave group. Reviewed `0034/0036` supplies the
full axisymmetric linearized Euler expression and the exact mixed-Casimir/end
rows; it does not supply `0037`'s projection or resolvent estimates.

The two historical source statements are correctly limited. Crossref's
primary metadata/abstract for S. Leibovich, *Weakly non-linear waves in
rotating fluids*, JFM 42 (1970), DOI
`10.1017/S0022112070001611`, states that the tube KdV method fails as the wall
moves to infinity and is replaced by a singularly derived
integro-differential equation. Cambridge's primary abstract for F. J. Higuera
and J. Jimenez, *Solitary waves on a vorticity layer*, JFM 264 (1994), DOI
`10.1017/S0022112094000674`, concerns contour-dynamics layers/tubes and reports
KdV with small logarithmic corrections for its axisymmetric family. Those
papers corroborate the singular-limit phenomenon but supply neither this
smooth pure-swirl carrier nor its dynamic complement estimate. The direct
source-transfer verdict is therefore correctly blocked.

## Unit A: threshold dispersion, compatible KdV, and scalar restoring calculus

Let `mu=lambda(c)>0`, `lambda(c_0)=0`, and
`sigma=1/lambda'(c_0)>0`. With `f_0` normalized in `L2(dr/r)`, the exact
rank-one exterior eigenvalue row gives

    lambda(c_ph(k))+k^2+T_R(k)f_0(R)^2
       +o(k^2 log(1/|k|))=0.

Taylor expansion at `c_0` therefore yields

    c_ph(k)-c_0
      =-sigma f_0(R)^2 k^2 log(1/|k|)+O(k^2).

The signs agree with the exact `0030` result that every nonzero column phase
speed is below `c_0`. Choosing the large root

    mu L_mu^2=f_0(R)^2 log L_mu

and setting `X=(z-c_0t)/L_mu`, `T=mu t/L_mu`, `k=kappa/L_mu` gives

    (omega-c_0k)L_mu/mu -> -sigma kappa^3.

The `kappa^3 log|kappa|/log L_mu` remainder is uniformly small on compact
`kappa` sets, including zero by continuity. Hence the limiting linear scalar
equation is `A_T+sigma A_XXX=0`.

The exact column perturbation equations and critical vector are correctly
written, and the accepted steady Euler reduction gives

    beta=(1/2) integral J_c0 f_0^3 dr/r>0.

The original body called this an executed time-dependent projection, but it
did not display the adjoint pairing or complement normalization. The bounded
correction now states the exact result: assuming a local conservative KdV
form, the already established steady profile fixes the compatible equation

    A_T+sigma partial_X(A_XX+beta A^2)=0.

A traveling profile `A(X-sigma T)` reduces this equation to
`A-A_XX=beta A^2`, so the coefficient and sign are fixed without fitting.
This is an exact compatibility/solvability identity, not a completed dynamic
Euler projection.

The scalar calculus is exact. The positive profile

    A_*=3 sech^2(X/2)/(2 beta)

is critical for `E_KdV+P_KdV`. Its Hessian

    L_*=-partial_XX+1-3 sech^2(X/2)

has the simple eigenvalues `-5/4`, `0`, and `3/4`, with the zero generated by
translation and essential spectrum `[1,infinity)`. The threshold value `1`
listed by the script is not mislabeled as a fourth discrete eigenvalue in the
prose. For the speed-`v` family,

    integral A_v dX=6 sqrt(v/sigma)/beta,
    P_KdV(A_v)=3(v/sigma)^(3/2)/beta^2,

and both derivatives in `v` are positive. One negative direction, the
translation zero, and the positive momentum slope yield the declared scalar
constrained coercivity. The two proposed Euler constraint pullbacks and the
time-dependent complement estimate remain open.

The threshold verifier correctly confirms the profile, traveling identity,
two integrals, slopes, and Pöschl--Teller values. Its terminal string says six
checks while the source contains several assertion nodes; that tally is only
provenance. More importantly, it contains no Bessel dispersion, radial
projection, or Euler remainder predicate and is credited only as scalar exact
corroboration.

**Unit A verdict: established after bounded scope/sign correction at
`EXACT_THRESHOLD_DISPERSION_AND_SCALAR_KDV_COMPATIBILITY`.** The dynamic
full-Euler projection and modulation approximation are blocked with their
adjoint-pairing and long-time complement construction named.

## Unit B: branchwise local energy and noncritical gap

For a positive-frequency column branch in the solitary frame,
`nu_n(k)=k c_n(k)-ck` and the exact group-speed ceiling gives

    nu_n'(k)=v_g,n(k)-c<=-(c-c_0)<0.

The monotone change of variable `k -> nu_n(k)` in time Plancherel is valid.
With the Fourier convention printed in the body it gives exactly

    integral_R |u(z,t)|^2 dt
       =2pi integral |a(k)|^2/|nu_n'(k)| dk,

and hence the local-window estimate with loss `(c-c_0)^-1`. Negative temporal
branches have the opposite group velocity and a stronger one-way separation.

For the noncritical complement, min--max monotonicity of the radial form gives
`c_n(k)<=c_n(0)` for `n>=2`; simplicity and the positive finite-interval gap
give `c_2(0)<c_0`. Thus one may take a fixed positive
`delta_2<=c_0-c_2(0)` and obtain `c-v_g,n(k)>=delta_2` for all noncritical
oscillating branches and all sufficiently small solitary members. The
stationary kernel is transported by `-c partial_z` and obeys the same type of
local estimate directly. The measurable spectral resolution requirement is
stated; no unjustified formal radial-mode sum is used.

The `O(mu L_mu)` coefficient accumulation on a fast crossing is consistent
with the accepted `0027` profile scale, while the critical branch retains the
explicit `O(mu^-1)` residence-time loss. The body correctly leaves the
actual solitary `Pi/Q` graph and off-diagonal estimate open.

**Unit B verdict: established as stated at
`EXACT_FREE_COLUMN_BRANCHWISE_LOCAL_ENERGY_AND_NONCRITICAL_GAP`.** It licenses
the fast-channel input to a future block theorem, not an all-time solitary
group.

## Unit C: one-pass geometry and the ray exponent

For the exact solitary member in its translating frame,

    (u_c)_z=-c+f_r/r<=-c_0/2

for sufficiently small `mu`. Therefore `z(t)` is strictly decreasing and
covers the localized axial region at most once; no material/covector orbit can
return to the excitation. The accepted weighted exterior reconstruction gives

    integral_R sup_r |grad delta u_mu(r,z)| dz
       <=C mu L_mu,

including the Bessel tail. Dividing by the strict axial speed produces the
same bound along every ray, and `mu L_mu ->0` follows directly from the scale
equation.

The column flow has fixed radius and angular rate `Omega(r)`. In a co-rotating
frame its deformation, cotangent, and principal full-pressure amplitude
systems are polynomial shear systems and have zero exponential exponent. The
solitary coefficients approach column coefficients at both axial ends and
their added strain is integrable along a trajectory. This preserves zero
exponential rate for every fixed sufficiently small member:

    limsup_|t|->infinity |t|^-1 log ||BAS_mu(t)||=0.

The original prose promoted this input to a literal finite scattering matrix.
That requires polynomially conjugated or weighted integrability and was not
proved. The corrected body now withholds it. It also correctly separates a
zero principal-ray exponent from high-frequency point spectrum,
Birman--Schwinger poles, finite/intermediate modes, and the full nonlocal Euler
semigroup.

**Unit C verdict: established after bounded scope correction at
`EXACT_ONE_PASS_MATERIAL_GEOMETRY_AND_ZERO_PRINCIPAL_RAY_EXPONENT`.** A
high-frequency full-operator exclusion remains open.

## Unit D: exact mixed-Casimir flux, counterexample, and repair

With `G_D'(xi)=xi D(xi)` and `q_D=zeta D(xi)`, the exact Euler equations give

    partial_t q_D+{Psi,q_D}=(2/r^4)partial_z G_D(xi).

The cylindrical divergence identity

    r{Psi,q}=partial_z(Psi_r q)-partial_r(Psi_z q)

has the printed signs. Axis regularity and exterior irrotationality remove the
radial trace, up to the stated `z`-independent background subtraction, and
yield

    partial_t M_D+partial_z F_D=0,
    M_D=integral_0^infinity r zeta D(xi)dr,
    F_D=integral_0^infinity
          [Psi_r zeta D(xi)-2G_D(xi)/r^3]dr.

Linearization retains all three variations in the flux, including the Hodge
term `(K eta)_r`, the `D' chi` contribution, and
`-2 xi_s D(xi_s)chi/r^3`. The executable symbolic oracle checks the transport,
divergence, chain-rule, and regular-label integration-by-parts identities; it
does not stand in for topology or continuity.

The `D=1` supervisor counterexample is exact and decisive. For
`psi_epsilon=epsilon^2 F(r/epsilon)g(z)`, `F(s)=s^2 chi(s)`, `g=G'`, direct
change of variables gives

    ||v_e||_rho^2=2pi rho_m[epsilon^2 C_F||g||_2^2
                              +epsilon^4 D_F||g'||_2^2],

while radial integration gives

    M_1(v_e)=2g-epsilon^2 B_F g''.

Thus the input tends to zero in both energy and the column graph norm while
the coordinate does not. The explicitly constructed compact azimuthal
generator is smooth at the axis and satisfies
`curl(a_e cross omega_0)=delta omega_theta e_theta`; the counterexample lies
in the true compact-generator coadjoint core and has zero integrated mixed-
Casimir variations. The `D=1` bounded-coordinate route is therefore refuted,
not merely unsupported.

For nonnegative `D` supported strictly inside a regular label annulus,
integration by parts gives

    M_D=integral (psi_r/r) partial_r[D(L)]dr
          -partial_zz integral psi D(L)/r dr.

Both radial weights are compactly supported away from the axis and flat edge.
On a fixed low-frequency window this is bounded by the physical energy, and
on the one-dimensional zero-frequency critical space

    m_D,0=integral Phi f_0 D(L)/(c_0^2 r)dr>0.

This is a genuine same-space repair. It supplies a bounded critical coordinate
and an exact output derivative, but not the physical left/right adjoint
oscillator pair, a uniformly bounded spectral projection, or cancellation of
the independent `mu^-1` loss.

**Unit D verdict: the `D=1` route is refuted by the named accessible
axis-concentration mechanism; the compact regular-label repair and exact local
flux are established; conversion to the true adjoint pair and scaled
cancellation is blocked with that construction explicitly named.** No part of
this ledger is a global no-go for the carrier.

## Unit E: abstract Kato-supersmooth block lemma

The algebraic sandwiched-resolvent identity has the correct order and sign.
For `A_0=A_P direct-sum A_Q` and `V=B^*CB`,

    B(lambda-A_0-V)^-1B^*
      =[1-B(lambda-A_0)^-1B^*C]^-1
        B(lambda-A_0)^-1B^*.

If the diagonal groups are uniformly bounded, the translation kernel is
modulated, the perturbation has a closed-generator or valid Kato-resolvent
realization, both forward and backward inhomogeneous supersmooth estimates
hold, and

    sup_(Re lambda !=0)||B(lambda-A_0)^-1B^*||<=K^2,
    K^2||C||<1,

then the Neumann inverse is uniform. Standard Duhamel/resolvent machinery
under exactly those hypotheses gives the conditional bounded wave operators,
uniform group bound on the declared subspace, and absence of spectrum in the
covered open half-planes.

The bounded correction adds the previously omitted realization and
inhomogeneous hypotheses. It also repairs the scalar generator convention:
in the right-moving profile coordinate `y=X-sigma T`, linearization of the
printed KdV equation gives

    A_P,0=+sigma partial_y L_*.

The spectrum and constrained coercivity of `L_*` imply scalar energy
boundedness but do not themselves prove Kato supersmoothness for the
nonselfadjoint generator `partial_y L_*`. The corrected document now treats
that scalar resolvent estimate as an additional open hypothesis. It likewise
withdraws the invalid ray-to-pole inference.

**Unit E verdict: established as a conditional abstract
Kato-supersmooth perturbation lemma after bounded domain/sign correction.**
Its P-channel supersmooth estimate, Euler generator realization, full-pressure
factorization, high-frequency operator exclusion, and `O(mu L_mu)` coupling
bound are not established and cannot be inherited from Units B--D.

## Evidence map and oracle audit

| Evidence | Proposition supported | Role | Limit |
|---|---|---|---|
| Accepted `0027/0028` operator/exterior proof | Exact Bessel symbol, log scale, weighted solitary reconstruction | Accepted supplier | No dynamic projection, ray-to-spectrum theorem, or stability |
| Accepted `0030/0035` column proof | Positive column group, phase/group-speed bounds, radial gap inputs | Accepted supplier | Column, not solitary operator |
| `threshold-reduction.md` | Threshold asymptotic and exact scalar KdV compatibility/calculus | Direct exact/asymptotic proof after correction | No time-dependent Euler projection or remainder |
| `verify_threshold_reduction.py` and receipts | Profile, invariants, slopes, Pöschl--Teller values | Exact scalar corroboration | Does not derive dispersion, beta, projection, or Euler transfer |
| `scattering-block-reduction.md` | Free branchwise Plancherel and gap reduction | Direct exact proof | `Pi/Q` graph and coupling remain open |
| `ray-exclusion.md` | Strict no-return and zero per-ray BAS exponent | Direct analytic proof after correction | No finite scattering limit or full-operator spectral exclusion |
| `mixed-casimir-flux.md` | Local conservation and regular-label density formula | Direct exact proof | No spectral-dual conversion |
| `supervisor-projection-analysis.md` | Accessible `D=1` discontinuity and compact-label repair | Exact counterexample and constructive repair | Repair is low-frequency/critical-coordinate scoped |
| `verify_mixed_casimir_flux.py` and repaired receipt | Four local algebra/integration identities | Exact symbolic corroboration | Does not prove accessibility, continuity, or adjoint status |
| `block-scattering-lemma.md` | Small supersmooth block implication under explicit hypotheses | Conditional abstract proof after correction | Hypotheses are not proved for the Euler blocks |
| Leibovich and Higuera--Jimenez primary abstracts | Different wall/contour singular-limit results | Source applicability | No theorem transfer to this field/domain |

The recorded threshold run first failed because the symbolic hyperbolic
integral path did not simplify; the repaired exact run exited zero and printed
the scalar tally. The mixed-Casimir script's repaired run exited zero and
prints four checks. These histories are provenance. No numerical observable or
small ratio is used, so `small-ratio-numerics` does not bind.

## Findings and bounded correction

| Finding | Direct evidence | Minimum correction | Upgrade path |
|---|---|---|---|
| The time-dependent quadratic projection was asserted but not executed. | The original A body jumped from the exact steady coefficient and compatible profile to the dynamic KdV equation without an adjoint pairing/complement calculation. | State exact branch compatibility only. | Construct the physical adjoint mode, calculate the projected quadratic source, and prove the complement normalization/remainder. |
| The original KdV block generator had the wrong sign in the natural right-moving coordinate. | Direct linearization in `y=X-sigma T` gives `a_T=+sigma partial_y L_*a`. | Use the plus sign in both block documents, or declare a reflected coordinate. | The corrected documents now use the plus sign consistently. |
| Integrable one-pass strain did not alone prove a literal scattering matrix between polynomial shear backgrounds. | Such a limit needs polynomially conjugated or weighted integrability; the body supplied only the unweighted `L1` bound. | Retain zero exponent and withhold the matrix limit. | Prove the missing weighted conjugated integral. |
| The abstract block paragraph omitted realization/inhomogeneous hypotheses and overread the ray result as a pole exclusion. | A sandwiched identity for a formal sum does not create a generator; zero BAS exponent does not exclude full-operator point spectrum. | Declare the generator/Kato realization and two-sided estimates; leave scalar supersmoothness, Euler factorization, and high-frequency spectral exclusion open. | Prove the scalar and Euler resolvent estimates on their exact domains. |

All four findings were repaired in the one bounded correction package. No
equation-level oracle changed and no unchanged execution needed replay.

## Correction check

The one bounded correction check passes. The supplied post hashes match the
filesystem and the append-only receipt pins all pre/post hashes. Direct
recalculation confirms that the right-moving coordinate produces the new plus
sign. The corrected A wording no longer claims a missing projection; C keeps
strict no-return, integrable strain, and zero ray exponent while explicitly
withholding a scattering matrix; and E now assumes or constructs a closed
generator, requires both inhomogeneous time directions, labels P-channel
supersmoothness open, and disclaims ray-based pole exclusion. The exact D body
was unchanged. The receipt records scoped `git diff --check` exit zero.

No new scientific claim was introduced by the correction. Postfreeze
adjoint-transfer artifacts remain excluded and supply no evidence here.

## Four-axis decision and frontier

The final decision is:

- Verification: exact analytic/asymptotic derivation, exact accessible
  counterexample, constructive compact-label repair, and conditional abstract
  operator proof; scripts are limited corroboration.
- Review: established after one bounded scope/domain/sign correction package
  at the A--E component scopes above.
- Compatibility: compatible active-proposal evidence on the accepted
  `0027/0028` solitary family and reviewed column/canonical suppliers.
- Epistemic: exact or controlled asymptotic for A--D at the declared scalar,
  free-column, ray, and local-flux domains; explicitly conditional for E.
- Relationship: no component inherits another component's missing projection,
  adjoint, resolvent, or full-pressure bridge.

The strongest next construction is the actual bounded change from the compact
regular-label density/flux pair to the physical two-branch adjoint oscillator
pair, together with P-channel supersmoothness and the full-pressure Bessel
`Pi/Q` factorization showing an `o(1)`—proposed `O(mu L_mu)`—sandwiched coupling.
That construction would determine whether the abstract lemma applies to the
linear axisymmetric Euler operator. The separate time-dependent Euler
projection/remainder is also required before the scalar KdV dynamics becomes
an Euler modulation theorem.

Full all-time Euler stability, nonlinear P2, nonlinear global regularity,
active `0040`, nonaxisymmetric control, particle mechanics, spin, statistics,
quantization, electron/neutrino identification, and parent-campaign completion
remain outside this review and are not implied.
