# P253/0104 independent final review of P253/0095

Reviewer: `particle-balance-review`

Target owner: `root`

Activated review README SHA-256:
`acf0af844fb373ce85669f788a08e1af665e2849289f4325fc06c1eb3a3905b0`

Frozen pre-review target manifest SHA-256:
`b305d2ea6ded594d24a20384789f03c86d282d175a4cf3c749a52fbc53280211`

Final post-correction target manifest SHA-256:
`22f0aeffa4c0b3a74b64a266eb3185beb066fc1178cdebd171c2883bcdc9f5c6`

## Boundary and provenance

This was a content-blind, independent, non-author review of final P253/0095.
The activated README matched its centrally registered hash, and
`activation-schema.exit` contained exactly `0`. P253/0088, P253/0094 and
P253/0095's parent, particle, P2 and P4 conclusions were outside the review.

The original target selectors all matched the activated inventory. The review
used the cached Hattori--Fukumoto primary PDF only after activation; its
SHA-256 is
`c6a35c44d55d26b8bb3e8fc9e55b29b638d0983b646b17eb522a4192b6444fb9`,
matching the target's source receipt. The Recovery-7 exact verifier has fourteen
derived predicates plus its final tally line, empty stderr and exit zero. The
focused Recovery-7 receipt has three passing API tests, empty stderr and exit
zero. They support algebra and regression respectively; neither proves the
compact-interior transfer, graph-domain Egorov estimate or nonlinear leaf
completion. No unchanged oracle or production numerical run was repeated.

The substantive audit found one bounded defect in Unit G and the combined
nonlinear verdict. The correction receipt
`0104-unit-g-bounded-correction-receipt.md`, SHA-256
`1e6afd6f20c9d41e850f1a3d8c42bab06226d1de74f9c71ba7c9817b37c82a5b`,
applied the minimum scope repair. The correction-only check verified its six
post-correction core hashes. Units A--F, the API semantics, tests and exact
verifier did not change.

Checkpoint staging subsequently exposed surplus blank lines at EOF in
`packaging-test-failures.md`, `recovery-audit-6-receipt.md` and the importable
API. The byte-only cleanup is pinned by
`checkpoint-whitespace-receipt.md`, SHA-256
`c67e2653b3867837adf04fc8bea72ae6605d75d44dad93b029ef4142a63bc338`.
It changes no source token, API behavior, equation, predicate or verdict. The
pre-cleanup manifest was
`dc15211dd0a2b4799a72c89446592202cf3256898f896118ad31915067bfaba4`;
all 63 entries of the final refreshed manifest verify.

## Source and authority audit

Hattori--Fukumoto equations (7)--(9) supply the trace-free amplitude system and
Hill elimination. Their equation (10) identifies the returned-covector
condition. Equations (22)--(24) supply the general leading-flow potential,
resonance condition and first-curvature growth formula. Their discussion after
equation (20) is load bearing: a merely perturbatively periodic covector drifts
out of resonance and produces only saturation. It does not license repeated
circuit growth by itself.

The reviewed Cao and charged Euler--Maxwell suppliers provide the fixed-member
carrier, compact-core regularity and the coupled equations only at their pinned
scopes. They do not supply the first-integral phase, finite-row submersion,
relative-total-momentum domain, nonlinear leaf path or global P2 conclusion.
Those bridges therefore had to be checked in P253/0095 itself.

## A — HF elimination, smooth-center coefficient and saturation boundary

For the source trace-free system

    p'=A p+B q,
    q'=C p-A q,

with the material coefficient `C` constant on the orbit, elimination gives

    q''+[A'-A^2-BC]q=0.

The signs agree with the source Hill equation. On the smooth-center branch,
`zeta/Omega=2+O(s^2)` and first resonance gives
`cos^2(chi)=1/16+O(s^2)`. Substitution in the source general-profile
coefficient gives forcing `15 Omega_c/128` and the paired physical exponent

    gamma_hat=delta |U_theta^(0)(s)| 15/256
              +O(delta s^3+delta^2).

The target does not turn this local coefficient into repeated growth before
repairing covector return. It accurately preserves the source warning that a
nonreturned perturbative covector exits the tongue and saturates.

**A verdict: established at the source-derived smooth-center first-curvature
scope.**

## B — exact first-integral phase and returned Floquet pair

In a regular action tube the phase

    phi=F(P)+ell theta,        ell in Z minus {0},

is globally circle-valued after exponentiation by integer `N`. Because the
no-swirl relative flow preserves `P` and has no toroidal component,

    (partial_t+W dot grad)phi=0,
    k=dphi,
    kdot=-(grad W)^T k.

Thus `k` and the physical action frame return exactly after a closed circuit.
This directly removes the post-equation-(20) HF saturation mechanism rather
than assuming it away. The continuous parameter `F'(P)` tunes the resonant
angle; the displayed nonzero detuning derivative permits the diagonal
first-order term to be canceled while periodic gauge changes leave the
off-diagonal Melnikov pairing invariant. At `g=0` the exact returned monodromy
has a reciprocal real hyperbolic pair with logarithms
`+/-delta |b_HF|T_*+O(delta^2)`. The correction correctly limits exact
reciprocity to this uncharged pair. At nonzero small `g`, equation (31) gives
only `det M=1+O(g^4)`; continuity still preserves two real hyperbolic
multipliers and an expanding one, but not exact reciprocity.

**B verdict: established, with exact reciprocal return at `g=0` and only
continuity-persistent hyperbolicity at nonzero `g`.**

## C — fixed-member compact-interior Cao transfer

The transfer is performed on a tube compactly contained in the positive
regular core. The normalized Cao profile converges there in the differentiable
topology needed by the Hill coefficients; affine logarithmic and recentering
rows are removed before comparison. The source Appendix-C cancellation then
gives the fixed-tube asymptotic

    gamma_hat,epsilon/(delta s)
      =15 A^p/512+O(s^2)+o_epsilon(1).

This establishes a positive first-curvature coefficient for every sufficiently
thin member on the selected compact regular orbit. It neither uses nor claims
endpoint, free-boundary, collar or uniform global `C^1` control.

**C verdict: established at fixed-member compact-interior scope.**

## D — DA/tag/Gauss packet and weighted input topology

At a returned resonant covector, the displacement polarization

    d=(k cross a)/(k dot bar_omega)

gives the packet

    xi_N=N^(-1) chi_tube Re[d exp(iN phi)],
    q_N=curl(xi_N cross bar_omega).

The principal Hodge velocity is `a`; an equivariant lower-order Bogovskii
correction preserves the nonzero toroidal character. The tag and electric
rows

    eta_N=-xi_N dot grad chi_g,
    e_N=-(g/epsilon_EM)grad(-Delta)^(-1)eta_N

satisfy the linear material and Gauss constraints. The tag has zero mean, so
its electric field starts at dipole order. In three dimensions
`<x>^(2 sigma)` is an `A_2` weight precisely for `|sigma|<3/2`; combining this
with monopole exclusion gives the nonempty interval
`1/2<sigma<3/2`. Compact curl supplies the fluid zero Fourier row. After
complete-state normalization, vorticity is the order-`N^s` row, tag and Gauss
are lower WKB order, and the sequence is weakly null in the declared graph
domain.

This establishes the linear DA/tag/Gauss input packet. It does not assert that
the nonlinear curve fixes every conserved finite row.

**D verdict: established at the weighted global linear-input and raw
coadjoint-tangent scope.**

## E — constrained transverse-Maxwell block and observed Egorov lift

After longitudinal Gauss reconstruction, the two transverse Maxwell
polarizations have roots

    -i c_g k_z +/- i c_EM |k|.

On the slow first-integral phase, strict subluminality gives the exact gap

    dist(0,spec D_EM)>=(c_EM-|c_g|)|k|.

The displayed Sylvester equations therefore have order-minus-one solutions
which remove both slow/Maxwell order-zero off-diagonal blocks. The target keeps
the longitudinal electric component as a slaved order-minus-one row and uses
the same compact cutoff, Leray/Hodge continuation and action-flow graph for
the remainder. Smooth separated kernels plus nonstationary phase control the
collar/exterior contribution. What is proved and used is the localized
fixed-time estimate

    sup_(0<=t<=T)||J_loc[S_g(t)Q_N-Q_N^BAS(t)]||=O_T(N^-1),

not propagation of a positive global weight by free Maxwell evolution and not
a growing-time uniform theorem.

**E verdict: established at the strict-noncharacteristic, fixed-time,
fixed-observation graph scope.**

## F — fixed observed propagator

`J_loc` is fixed once, includes the declared instantaneous local fluid, tag
and field rows, and is independent of circuit count. For each integer `j`, the
packet tube and frequency threshold may depend on `j`, while the central
returned symbol ratio does not. The weak-null sequence and fixed-time Egorov
estimate give

    ||J_loc S_g(jT_*)||_ess,
      X_in^s -> X_DA,loc^s/sym >= c_obs lambda_+^j,

with one `c_obs>0` and one `J_loc` independent of `j`. This is exactly the
frozen global-weighted-input to local-same-`H^s` observed operator. It is not a
claim about `||S_g||_ess`, and it does not use a single packet uniformly as
time grows.

**F verdict: established for the fixed observed propagator.**

## G — exact leaf path and finite-row submersion

The raw curve

    omega_tau=(Phi_tau)_*omega_g,
    chi_tau=chi_g composed Phi_tau^(-1),
    B_tau=B_g,
    E_tau=E_g-(g/epsilon_EM)grad(-Delta)^(-1)(chi_tau-chi_g)

is exact for coadjoint/Kelvin-circulation data, the tag distribution, charge,
divergence constraints, the difference Gauss law and zero monopole. This is a
useful exact same-field construction.

It is not an exact curve in the README's complete fixed-row leaf. The original
proof assumed continuity and independence of the impulse, center and
total-momentum rows and then invoked density; it exhibited no
carrier-specific compact generator witnesses or invertible row matrix. Density
cannot create the missing independence.

There is also a genuine domain obstruction. The inherited physical row

    P_tot=rho_m integral u dx
          +epsilon_EM integral E cross B dx

exists only under sufficient convergence. A compact mean-zero vorticity
perturbation has `qhat(k)=O(|k|)` and generically
`uhat=i k cross qhat/|k|^2=O(1)`, corresponding to an `O(|x|^-3)` Hodge tail.
That tail need not be `L^1`, so `integral u` is not a continuous functional on
the declared `sigma<3/2` completed topology. Relative first-moment/impulse
cancellation producing `qhat=O(|k|^2)`, plus the field integrability ledger,
would be one repair. Alternatively one may prove total momentum automatic or
redundant for this curve. Either route still needs explicit compact smooth
witnesses and a nonzero retained-row determinant after the actual symmetry
quotient.

The bounded correction now states this boundary and leaves its formal witness
matrix as the precise next construction. It does not delete `P_tot` and thereby
silently enlarge the frozen leaf.

**G verdict: blocked at the relative-total-momentum domain and
carrier-specific finite-row submersion.**

## H — fixed-thinness small-charge conclusion

For each fixed sufficiently thin `delta`, finite-dimensional continuity gives

    ||M_(delta,g)-M_(delta,0)||
      <=C_mon(delta)g^2+o_delta(g^2).

Hence a real expanding multiplier and the Unit-F observed essential-norm bound
persist for signed `0<|g|<g_0(delta)` whenever
`C_mon(delta)g^2<c_HF delta`. No bound uniform in `delta` is proved, so there is
no `g^2<=c delta` wedge and no result in the complementary regime.

The linear observed theorem is established. Turning a hypothetical nonlinear
uniform-Lipschitz estimate into a bound on `J_loc S_g` requires a smooth exact
path in the same frozen leaf and same-target directional differentiation.
The differentiation estimate is a valid conditional bridge, but Unit G does
not supply its antecedent. The nonlinear fixed-member refutation therefore
does not follow.

**H verdict: established for fixed-`delta` small-charge linear observed
growth.** Unit G leaves the separate downstream nonlinear implication
unlicensed; that does not narrow Unit H's verdict.

## I — complementary charged-column route

The column reduction gives the exact slow block with square
`-x^2 Z_tau R_tau I` and the dimensionless classifier

    K(s,tau)=4 Z_tau R_tau/Omega_tau^2.

Thus `Z_tau R_tau<0` is direct column hyperbolicity, while `K>=1` permits the
first integer curvature resonance at `|x|=1/sqrt(K)`. The ordered Melnikov
functional `b(tau)` has `b(0)!=0`; on any connected analytic charged-column
cell this nonvanishing persists locally. A restoring candidate would instead
need common positive-Krein signs and a uniform
`0<kappa_0<=K<=1-eta` throughout the declared core. Pointwise `0<K<1` is not a
global gap, and the `x->0` frequency degeneration still needs graph/resolvent
control.

No delta-uniform charged Cao-to-column free-boundary branch, curvature jet or
whole-core coverage is constructed. The equations are a concrete
continuation target, not a result in the complementary charge regime.

**I verdict: blocked at the charged-column/free-boundary transfer and
whole-core spectral coverage, with its exact classifier and Melnikov target
established.**

## Combined verdict

The strongest supported theorem is the A--F linear statement: for every fixed
sufficiently thin Cao member and sufficiently small signed nonzero charge
satisfying the displayed fixed-member gap inequality, an exact returned DA
phase carries a real expanding hyperbolic multiplier into a fixed-`J_loc`
observed-propagator essential-norm lower bound from the weighted global
complete-state graph domain to the local same-`H^s` quotient.

The full source-to-exact-leaf chain is not complete. Consequently the frozen
fixed-member nonlinear uniform-Lipschitz refutation is blocked, not refuted or
established. This does not weaken A--F and does not decide the complementary
charge regime, Hill's full constrained Hessian, nonlinear instability, P2,
P4, particles or the parent campaign.

The next load-bearing construction is Unit G: type a relative-total-momentum
domain or prove that row redundant, then exhibit carrier-specific compact
generator witnesses with an invertible impulse/center/remaining-row matrix,
including the Gauss-reconstructed electromagnetic contribution wherever that
row is retained. Independently, the complementary route needs the charged
column/free-boundary branch, uniform curvature jet, whole-core `K` coverage,
and `b` evaluation over every resonant point.
