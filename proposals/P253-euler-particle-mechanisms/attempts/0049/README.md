# P253/0049: physical analyzer and first-event measurement bridge

Root preregisters a failure-derived P4 construction after `0043/0046` showed
that a KKS phase plane, geometric quantization, exchange topology, and complete
classical current tomography still do not make a physical two-state quantum
system. This attempt asks for an explicit state/analyzer/detector dynamics. It
does not allow the words measurement, random phase, or Poisson to hide a new
law.

The state input is a normalized two-component complex amplitude only if an
actual same-carrier positive mode pair and its physical energy norm have been
constructed. The analyzer must be a same-substrate evolution implementing a
declared `SU(2)` basis change or its controlled approximation. Its two output
channel energies must be derived from the physical Euler energy/current, not
assigned as probabilities.

## Route A: deterministic exclusive threshold detector

Couple the two analyzer outputs to two initially identical localized Euler
detector cells. Derive the full Hamiltonian/retained-state equations, including
the unresolved initial state from `0042`, and test whether the dynamics has
two exclusive capture basins. For an input ray `psi` and analyzer basis
`e_0,e_1`, compute the basin measure under one fixed substrate invariant law.
This route succeeds only if

    Prob(outcome i)=|<e_i,psi>|^2

follows for every ray and analyzer orientation, sequential analyzer updates
are physical, one event excludes the other, and the same fixed law is used for
all settings. A threshold model that is tuned separately for each amplitude
or merely splits energy between both outputs fails the target.

## Route B: first-event Poisson competition

Test the exact conditional mechanism in which each output drives an
independent first-event clock of rate

    lambda_i=kappa I_i,
    I_i=|<e_i,psi>|^2.

The exponential-race identity would give

    Prob(i is first)=I_i/(I_0+I_1).

The calculation must retain zero-intensity channels, ties, detector reset,
and sequential measurements. Then audit whether the required temporal clocks,
their independence, and rate-current law follow from the accepted stationary
Poisson marking/mixing construction. Spatial Poisson marks or a random choice
of initial tube do not automatically provide temporal detector clocks. Any
unproved clock law is recorded as a new stochastic detector hypothesis.

## Route C: rare-event deterministic mixing limit

Replace the clocks by an actual deterministic Euler detector flow with small
target sets. Prove a joint hitting-time limit, including independence across
channels and intensity-proportional rates, uniformly over analyzer settings.
The limit must be coupled to a physical exclusive capture/reset event. A
generic appeal to chaos, ergodicity, or mixing is not a construction. Compare
this route with Route B because it could derive the same probability law from
substrate dynamics rather than postulate it.

## Action and propagation obligations

For every route, separately test whether a detector event transfers one fixed
physical action/energy packet and whether that packet scale is selected by the
same substrate. Continuous Euler similarity is an exposing mutation: a rule
that changes continuously under `u_AB=A u(Bx,ABt)` does not select a universal
action without a background or new invariant. Also retain the finite-speed
obligation: an effective analyzer band may have finite group velocity, while
the same-field Euler pressure remains elliptic and must be included in the
detector interaction.

## Success and verdict boundary

The strongest positive result may be an exact probability theorem conditional
on a plainly declared detector clock law. That advances the design without
earning LP4 from Euler alone. LP4 requires the actual same-substrate carrier,
analyzer, detector/reset dynamics, selected action scale, exchange character,
and finite-speed band on one state space.

Each route receives one verdict and failure generates the next concrete
detector concept. Active `0048` owns carrier geometry and may be used only at
its independently supported scope. This attempt makes no electron, neutrino,
nonlinear carrier-stability, or parent-completion claim. No production
numerics are registered; exact probability and energy-current algebra comes
first.
