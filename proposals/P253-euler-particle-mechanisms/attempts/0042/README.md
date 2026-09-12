# P253/0042: exact nonlinear material observables with retained Euler memory

Root preregisters a fixed synthesis theorem for obligation P1.  The inputs are
the independently reviewed `0001/0004` material balances and `0009/0011`
projected-memory identities, together with the exact whole-space Euler
pressure reconstruction.  Earlier bodies have already been inspected during
this campaign; no comparator blinding is claimed.  This attempt asks whether
those results close P1 at an honest retained-state or retained-memory scope.
It does not ask for a particle, carrier persistence, or a finite-dimensional
autonomous constitutive law.

## Frozen physical state and observables

Work on the common smooth-existence interval of constant-density Euler on
`R^3`, with `u in H^s_sigma`, `s>5/2`, and the decay needed for kinetic energy
and the pressure Poisson solution. The material tag is the indicator/smooth
weight of a bounded material region, or is compactly supported/integrable with
the finite weighted moments and boundary regularity required by `0001`; it is
transported by the same velocity. Every material observable below is scoped
to that finite-moment class.
The resolved observable vector contains exactly:

- tag mass and centroid;
- centered inertia/shape and its first moment;
- tag momentum and centered angular momentum;
- kinetic covariance and tag kinetic energy;
- the pressure force, pressure moment, pressure torque, and pressure-energy
  flux that occur in their exact rates.

Pressure is the same-field nonlocal quantity

    p=rho (-Delta)^-1 partial_i partial_j(u_i u_j),

with the decaying normalization.  It is not an independent fitted stress.
The resolved map is denoted `R(U)`, where `U` includes the Euler field and
the material tag.

## Route A: exact resolved/complement Markov system

Choose a finite-rank `L2`-orthogonal solenoidal projection `Pi` with smooth
range, bounded on `H^s_sigma` and smoothing `H^(s-1)_sigma` to `H^s_sigma`,
and write `u=v+w`. Retain the material tag and `w` as state. Derive
the exact coupled equations, pressure reconstruction, total energy exchange,
and every observable rate above.  This route succeeds only if reconstructing
`u=v+w` makes the system literally Euler and if the pressure/current and
ambient momentum rows are explicit. The full ambient field is always retained;
its integrated momentum is a scalar observable only under the additional
`L1`/decay hypothesis making that integral finite. Calling the full complement an
effective force without its evolution does not close this route.

## Route B: exact nonlinear observable memory

On a function class where the local Euler-tag flow and its Liouville
operators are defined, freeze a bounded linear idempotent `P_obs` on
observables and `Q_obs=1-P_obs`. These are distinct from the velocity-space
operators `Pi,Q`. Require a common invariant domain on which `L` and
`Q_obs L` generate the displayed evolutions and all products are defined.
Derive the Dyson/Mori--Zwanzig identity with all three terms:

    d/dt e^(tL) A
      =e^(tL) P_obs L A
       +e^(t Q_obs L) Q_obs L A
       +integral_0^t e^((t-s)L) P_obs L
          e^(s Q_obs L) Q_obs L A ds,

or the exactly equivalent convention obtained by a proved Dyson identity.
Apply it to each component of `A=R` only when the complete material vector
lies in that common domain, and show that the noise term is the unresolved initial
Euler state and the memory kernel uses the orthogonal dynamics, including
pressure.  The identity must be nonlinear at the observable-flow level;
the periodic shear convolution in `0009` is an exposing exact specialization,
not the general proof.

The activated README initially put `P_obs` on the full left side and on the noise
term. Since `e^(t Q_obs L)Q_obs L A` remains in the `Q_obs` range, that projected noise would
vanish and would not be the required trajectory identity.  The displayed
full-observable Dyson formula is the corrected frozen target.  Pre-correction
README SHA-256: `2002c78268aa436824d5cf1beb7ae805a662da3c09ca86d42893df8301c47569`.

## Success statement and exclusions

P1 is earned if Route A supplies the exact retained-state closure and Route B
supplies its equivalent exact observable-memory representation on a declared
local-flow domain, while `0001` supplies every listed material balance.  The
result must state that two fields with identical resolved local tag data can
have different pressure accelerations, so a memory/noise or complement state
is load-bearing.

This theorem would license an honest state/observation map for later carrier
and interaction work.  It does not license a finite-dimensional particle
equation, all-time Euler regularity, persistent localization, reciprocal pair
force, electric charge/current, quantum mechanics, electron/neutrino identity,
or parent completion.  Those remain P2--P7.

The strongest oracle is an algebraic Dyson rederivation plus the existing
exact nonlinear shear specialization and pressure-quadrupole counterexample.
Any reusable general definitions belong in an importable module with exposing
tests.  No production numerics or small-ratio quantity enters this attempt.
