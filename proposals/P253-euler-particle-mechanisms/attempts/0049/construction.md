# Exact first-event probability and the missing physical detector

## 1. Analyzer energy is not yet an event probability

Let `psi in C^2`, `||psi||=1`, and let a physical analyzer, if one is
constructed, implement `U in U(2)`. The two positive quadratic output
fractions are

    I_i=|(U psi)_i|^2,   I_0+I_1=1.                    (1)

Equation (1) is exact linear-wave energy algebra. On the accepted compact
Euler background it can only be consumed after an invariant positive mode
pair and an actual analyzer evolution have been constructed; the current
finite-window prepared response is not such an autonomous device.

For a deterministic detector with unresolved initial state `xi` distributed
by one fixed law `mu_D`, let `B_i(psi,U)` be its exclusive outcome basins.
The Born rule is precisely the dynamical statement

    mu_D(B_i(psi,U))=I_i                              (2)

for every state and analyzer setting. Neither Hamiltonian determinism nor
energy conservation implies (2). For example, a hidden threshold
`lambda in [0,1]` with outcome zero when `lambda<=I_0` gives probability
`F(I_0)`, where `F(x)=Prob(lambda<=x)` is the distribution function; it equals
`I_0` only after choosing the uniform law. This exposes where a probability
postulate can hide in an apparently deterministic threshold rule.

The exact retained-state theorem `0042/0045` makes the missing construction
concrete: detector basins must be derived from the coupled Euler state
`(v,w,tau)` and the same unresolved `w_0`, not assigned directly from (1).
No such exclusive capture/reset flow is present in the accepted suppliers.

**Route A verdict:** blocked by the actual analyzer, exclusive detector basins,
and invariant basin measure. The equivalence (2) is exact and prevents a
classical energy split from being relabeled as a one-event probability.

## 2. Route B: exact first-event competition

Declare independent waiting times `T_i` with exponential rates

    lambda_i=kappa I_i,   kappa>0.                     (3)

For total rate `Lambda=lambda_0+lambda_1`, independence gives

    Prob(T_0>t,T_1>t)=exp(-Lambda t).                 (4)

The density that channel `i` fires first at time `t` is
`lambda_i exp(-Lambda t)`. Therefore, for `Lambda>0`,

    Prob(i fires first)
      =integral_0^infinity lambda_i exp(-Lambda t)dt
      =lambda_i/Lambda
      =I_i/(I_0+I_1).                                 (5)

A zero-intensity channel has zero rate and never wins. Independent continuous
waiting times tie with probability zero. On a finite observation window,
both the probability of outcome `i` and the total click probability acquire
the common factor `1-exp(-Lambda T)`, so the outcome conditional on a click
still obeys (5). For normalized (1), it is exactly `I_i`.

If every rate vanishes, `Lambda=0`: no finite event occurs and the
conditional winner distribution is undefined.  The normalized analyzer in
(1) with `kappa>0` does not enter this boundary, but the declared clock theorem
and its API must retain it for unnormalized or empty input.

If the event dynamics additionally resets the carrier state to output basis
ray `e_i`, repetition in the same basis is certain. A later analyzer `V`
then has transition probabilities

    Prob(j|i)=|V_ji|^2.                               (6)

Equations (5)--(6) reproduce the projective finite-dimensional probability
calculus exactly. The reset is load-bearing: the clock race alone destroys no
coherence and defines no post-event state.

An exact many-to-one projective reset cannot be the autonomous time map of
reversible, measure-preserving Euler on a finite invariant phase volume.  It
is therefore an open-system/coarse-grained rule or external repreparation
hypothesis here.  An autonomous first-exit or scattering event can retain the
escaped degrees of freedom and remain reversible, but by itself it does not
reprepare the captured carrier for the sequential law (6).

This is a positive conditional mechanism, but its physical premises are not
Euler consequences. The accepted material currents are signed momentum and
angular currents, not nonnegative event intensities. The rule
`lambda_i=kappa I_i`, independent temporal clocks, exclusive capture, and
open-system or external reset are new detector laws. The scale `kappa`
cancels from (5) and therefore
cannot select `hbar`.

**Route B verdict:** established as an exact probability theorem conditional
on four explicit new hypotheses: physical analyzer energy channels,
independent exponential clocks with rate (3), exclusive first capture, and
projective reset. It is not an Euler derivation of measurement.

## 3. Route C: deterministic rare-event realization remains a PDE problem

A same-substrate derivation could replace the declared clocks by hitting
times of a deterministic mixing detector flow. It would have to construct
two physical target sets and prove, uniformly over `psi,U`, that their joint
rescaled hitting process converges to independent exponentials with rates
proportional to (1). It must then couple the first hit to an exclusive
finite-energy capture/reset solution.

The accepted spatial Poisson marking in `C-CST-008` does not do this: it marks
good cells in a stationary random field. It neither evolves detector cells in
time nor proves a joint hitting law. The compact assembly is stationary, and
the periodic-pair supplier is recurrent rather than mixing. Thus no current
source meets the required hypotheses.

**Route C verdict:** blocked by the named joint temporal hitting-time theorem
and its physical capture coupling. This is the positive same-substrate repair
of Route B, not a generic appeal to chaos.

## 4. Action scale and finite propagation remain separate

The winner probability in (5) is invariant under a common rate rescaling, so
it contains no action quantum. The accepted Euler similarity continuously
rescales circulation and KKS action while preserving dimensionless carrier
topology. A detector that transfers one fixed packet must therefore derive a
nonzero minimum packet from a background-selected scale or add it as a new
law; merely calling a click one quantum is circular.

Likewise, (1)--(6) contain no propagation theorem. An analyzer may use an
effective finite group-velocity band, but the full same-field pressure is
elliptic and its detector coupling must be retained. The first-event law does
not turn a Galilean Euler background into a Lorentz cone.

## 5. Joined verdict and continuation

The exact exponential-race construction shows a minimal way in which a
classical substrate could produce the Born fractions and sequential
projective transition table. It also localizes every new assumption. Current
Euler evidence supplies channel-energy and retained-state ingredients, but
not the autonomous analyzer, temporal clock law, exclusive reset, selected
action packet, exchange character, or finite-speed band.

The same-substrate continuation is Route C on an actual persistent carrier:
derive rare-event detector statistics from deterministic unresolved Euler
dynamics and couple the event to an exclusive topological capture. In
parallel, the compact `S2` carrier orbit sought by `0048` is the natural state
space for analyzer rotations. These two lines must meet on the same carrier;
otherwise the construction remains a conditional detector model rather than
LP4.
