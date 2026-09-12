# P253/0045 independent review of the 0042 retained Euler theorem

## Frozen transaction and correction provenance

This transaction independently reviews root-owned attempt `0042`. The
reviewer authored and implemented none of its construction, API, tests, or
correction. The frozen review README has SHA-256
`6743a12978fcad28d15ca3c900975f0135b897312b10414b1775c1a8f706314f`.
Central registration names `particle-balance-review`; the activation exit is
exactly `0`, with SHA-256
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`,
stdout
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
and empty stderr
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
The accepted authority base is release `v0.183.0`; the preregistration head
was `ca929ce584b462c033fde7866e6231076445ecd5`.

The substantive pass used these frozen target hashes:

| Artifact | Substantive-pass SHA-256 |
|---|---|
| `0042/README.md` | `9c37f046edde581edd3e03e0c69fadaffb5d5ef91c7b1001c87bc1797c013a87` |
| `0042/construction.md` | `6e9993547fea6617150e4a7f8f72cefa75dc443026b4215f1e173e2767d1d277` |
| `0042/result.yaml` | `d23ebb23b0dc755ecb83926580ab49a932d689bfafa6e78f8037e80e6e8d6424` |
| `0042/validation.md` | `e8efa13afae3422c4e4526394c9dfb35da08fed4709fe66233c97a612758de28` |
| `0042/readme-correction.md` | `3d48d6f4d23cee710c0a61de0646232e787a0f19c225d7f90a8dc60606b8f4f1` |
| `euler_retained_memory.py` | `75ab59d079de41481f1bb85d6f9a97a7d230ab00d3458236e43456f4f9362f0d` |
| `test_euler_retained_memory.py` | `46159071e1c7b90952643c2f24efb626ec3e2e26e63b435abea1f0dddf422e92` |

One bounded analytic-domain correction was requested. Its append-only receipt
is `0042/review-correction-0045.md`, SHA-256
`4dd1802f582a1033f242e0e37715b8bf8d175f8eec62f9d9e494fd6d0c39b0cf`.
The single bounded check confirms its complete before/after map:

| Artifact | Before SHA-256 | Corrected SHA-256 |
|---|---|---|
| `0042/README.md` | `9c37f046edde581edd3e03e0c69fadaffb5d5ef91c7b1001c87bc1797c013a87` | `416423a720c879f9fd98666292305e2d9492b685fad424a4f53e437f7310fa84` |
| `0042/construction.md` | `6e9993547fea6617150e4a7f8f72cefa75dc443026b4215f1e173e2767d1d277` | `b573d04f0d3a46819e569568c73c69c4d79a70cabb845f35e15dc56f93fa7061` |
| `0042/result.yaml` | `d23ebb23b0dc755ecb83926580ab49a932d689bfafa6e78f8037e80e6e8d6424` | `319b248ee41e84a9753198b57de016dcf98003c83744bb80d575ef8416a7ebd7` |
| `0042/validation.md` | `e8efa13afae3422c4e4526394c9dfb35da08fed4709fe66233c97a612758de28` | `226bfd5cdf9330a3e2adc55780d5b52adbb2b8e091757153f13a437db7ee7934` |

The correction changes no equation, API, test, or exact oracle. The unchanged
12-test receipt was therefore correctly not rerun. Active `0043`, P2/P4, and
all particle or quantum extensions remained outside the review.

## Strongest supported result

Let `s>5/2`, let `u` be a local divergence-free `H^s(R^3)` Euler solution,
and let the material tag define a bounded transported region or an
integrable transported weight with the finite moments and boundary
regularity needed by `0001`. Let `Pi` be a finite-rank `L2`-orthogonal
solenoidal projection with smooth range, bounded on `H^s` and smoothing from
`H^(s-1)` to `H^s`, and set `Q=I-Pi`, `v=Pi u`, and `w=Q u`.

Then `(v,w,tau)` obeys an exact Markov system equivalent to literal Euler.
For any prescribed `C^1` resolved history and `w_0 in QH^s`, the complementary
equation has a unique local response `W[v,w_0]` by the forced-Euler
transport/energy method. Requiring the resolved history to satisfy its
projected equation gives the exact causal retained-history representation.
The same full field determines the physical pressure, transports the tag,
and supplies every finite material balance reviewed in `0001/0004`.

The unresolved initial datum is load-bearing. It cannot generally be deleted
or inferred from the listed local tag moments. This is an exact
infinite-dimensional retained representation, not a finite-dimensional
autonomous particle law.

On a separate observable space, the full-observable Dyson identity is exact
for a bounded linear idempotent `P_obs`, `Q_obs=I-P_obs`, on a common domain
where `L`, `Q_obs L`, and all displayed compositions generate the required
evolutions. It is exact for the bounded/Galerkin and reviewed shear examples.
No generic infinite-dimensional Euler theorem generating
`exp(t Q_obs L)` is supplied, and the result correctly leaves that route
blocked rather than importing a formal semigroup.

## Pressure and retained-state calculus

With physical pressure and density `rho_m`, divergence of Euler gives

    -Delta p=rho_m partial_i partial_j(u_i u_j)
            =rho_m partial_i u_j partial_j u_i.

Thus, for nonzero Fourier frequency,

    p_hat(k)=-rho_m k_i k_j (u_i u_j)_hat(k)/|k|^2.

The sign, density factor, derivative order, and crossed contraction are all
correct. The API's `trace(Du*Du)` is exactly
`partial_j u_i partial_i u_j`; it is not the Frobenius norm
`trace(Du^T Du)`. Pressure is computed from `u=v+w`, never from `v` alone.

Writing `B(a,b)=-P_L[(a dot grad)b]`, direct projection gives

    v_t=Pi B(v+w,v+w),
    w_t=Q B(v+w,v+w).

Adding the rows reconstructs Euler because `Pi+Q=I`. Orthogonality, not mere
idempotence, gives the kinetic split. Integration by parts yields

    d/dt [rho_m ||v||_2^2/2]
      =rho_m integral w dot ((v dot grad)v)
       +rho_m integral w_i w_j partial_j v_i
      =-d/dt [rho_m ||w||_2^2/2].

Internal exchange therefore cancels in the full energy. The full ambient
field is retained for all stated `H^s` data. Its integrated momentum is a
finite scalar only under the corrected additional `L1` or equivalent decay
hypothesis.

**Retained-state verdict: established after the bounded domain correction.**

## Forced-Euler response and material-history substitution

For prescribed `v in C^1([0,T];Ran Pi)` and `w_0 in QH^s`, set
`u=v+w`. The complement equation is equivalent to

    u_t=-P_L[(u dot grad)u]+v_t+Pi P_L[(u dot grad)u],
    u(0)=v(0)+w_0.

The quadratic transport term has the usual one-derivative Euler structure.
The corrected smoothing bound controls the finite-rank term in `H^s`, so the
standard Euler energy/compactness argument gives
`u in C H^s intersect C^1 H^(s-1)` on a time controlled by the initial norm
and the prescribed `C^1` history norm. This is not an invalid Picard theorem
for a locally Lipschitz vector field on `H^s`. Applying `Pi` to the equation
gives `partial_t Pi(u-v)=0`; the compatible initial condition therefore
preserves `w in QH^s`.

Substitution produces

    v_t=Pi B(v+W[v,w_0],v+W[v,w_0]).

An arbitrary prescribed `v` does not become an Euler solution: literal Euler
is recovered only when this resolved equation also holds. On that common
local interval, the flow of `v+W` transports the same finite-moment tag and
the pressure is reconstructed from the same `v+W`. Consequently the moving
centroid, centering terms, pressure force/moment/torque, kinetic covariance,
and pressure-energy flux are retained rather than replaced by a fitted force.

The reviewed pressure-quadrupole example proves that identical local tag
moments can have different initial pressure acceleration. Independently, the
reviewed periodic shear pair has the same resolved initial datum and energy
but opposite resolved initial derivatives through its nonzero unresolved
initial term. These show missing-state dependence, not nonuniqueness of the
full local Euler solution and not an impossibility theorem for special
invariant families.

**Causal-history and load-bearing-state verdict: established after the
bounded domain correction.**

## Full-observable Dyson identity and domain boundary

Use the convention `(e^(tL)A)(U_0)=A(S_tU_0)`, so
`d e^(tL)A/dt=e^(tL)LA`. With the corrected observable projection,
the right-Dyson formula is

    e^(tL)=e^(t Q_obs L)
      +integral_0^t e^((t-s)L) P_obs L e^(s Q_obs L) ds.

Splitting `LA=P_obs LA+Q_obs LA` gives

    d/dt e^(tL)A
      =e^(tL)P_obs LA+e^(t Q_obs L)Q_obs LA
       +integral_0^t e^((t-s)L)P_obs L
          e^(s Q_obs L)Q_obs LA ds.

The operator order and signs are correct. At `t=0`, the right side reduces to
`P_obs LA+Q_obs LA=LA`. The noise term remains unprojected; placing `P_obs`
on the whole equation would erase it. Componentwise application to the full
material vector is supported only when those nonlinear pressure/tag
observables lie in the declared common domain.

**Dyson verdict: established under the explicit generator/domain hypotheses;
the generic infinite-dimensional `Q_obs L` Euler semigroup remains blocked
with its missing construction named.**

## API and oracle scope

`euler_pressure_poisson_source` implements the correct crossed pressure
source and states incompressibility as a caller hypothesis. The mutation test
distinguishes it from the Frobenius contraction. `linear_retained_memory`
implements, for a finite matrix generator and idempotent projection,

    P L P,
    P L Q exp(t QLQ) Q,
    P L Q exp(t QLQ) Q L P.

The symbolic two-block test derives the Markov, unresolved-initial, and
memory factors `a`, `b exp(dt)`, and `bc exp(dt)` with their signs and order.
The invalid-projection and density checks exercise bounded input contracts.
Together with the affected supplier tests, the pinned receipt is
`12 passed in 5.06s`, exit zero. These tests corroborate finite algebra only;
they do not prove forced-Euler well-posedness, tag integrability, or an
infinite-dimensional orthogonal-dynamics semigroup.

## Final verdict and next dependency

The joined `0042` claim is **established after one bounded analytic-domain
correction** at scope
`EXACT_LOCAL_RETAINED_EULER_STATE_AND_CAUSAL_HISTORY_WITH_CONDITIONAL_DYSON`.
All four requested domain repairs are supported and no further correction is
needed. The strongest result preserves an exact same-field local Euler state,
pressure, tag-observable map, and causal unresolved-history representation.

The next scientific dependency is not another closure disclaimer. It is to
instantiate this retained state and its unresolved initial datum on the same
independently supported persistent carrier used by subsequent interaction or
state/action work. A generic `Q_obs L` semigroup would require a separate
observable-space generator/domain construction; it is unnecessary for the
positive physical response map proved here.

This review does not license P2/P4, active `0043`, a finite-dimensional
autonomous material or particle law, all-time three-dimensional Euler
regularity, carrier persistence or stability, a reciprocal pair potential,
electromagnetic or weak current, quantum mechanics, electron/neutrino
identity, or parent-campaign completion.
