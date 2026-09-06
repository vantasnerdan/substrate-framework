# P253/0003 exact transfer, topology, scale, and bridge derivation

## 1. Objects and maximum verdict

This attempt concerns P0/P4 route R2 only. Its object is the map

```
F : phi |-> B=phi^*omega=i_u nu_g |-> u=(star_g B)^sharp          (1)
```

on Slobodeanu's spatial Riemannian domain. The maximum result available here
is an exact stationary correspondence plus an exact audit of what survives
under time dependence, topology, dimensional normalization, and physical
state-space quotienting. It cannot identify an electron or neutrino without
LP1--LP4 and the same-carrier P2/P3 inputs.

## 2. Actual time-dependent transfer test

Let `phi_t:M->N` be a smooth one-parameter family, let
`B_t=phi_t^*omega`, and define `u_t` by (1). Since `d omega=0`,

```
d B_t = 0,
div u_t = 0.                                                      (2)
```

If `W_t=partial_t phi_t` is the variation field along `phi_t`, the pullback
variation formula gives the exact kinematic tangent

```
partial_t B_t = d a_t,
a_t(X)=omega(W_t,dphi_t X),
partial_t u_t^flat = star_g d a_t.                               (3)
```

The time-dependent incompressible Euler equation in one-form notation is

```
partial_t u^flat + i_u d u^flat = -d Pi,
Pi=p+(1/2)|u|^2.                                                  (4)
```

### Test A: a family of instantaneous Slobodeanu critical points

At every time, Slobodeanu criticality supplies only

```
i_u d u^flat = -d(P_t o phi_t).                                  (5)
```

Equations (4)--(5) can both hold only if `partial_t u^flat` is exact.
It is also coclosed by (2). Hence it is an exact harmonic one-form. On a
closed connected `M` it vanishes; on Euclidean `R^3` it also vanishes under
the isolated finite-energy/decay conditions that exclude nonzero harmonic
gradients. Therefore

```
instantaneous sigma2 criticality + actual isolated Euler evolution
                         ==> partial_t u = 0.                     (6)
```

This rules out direct promotion of a time-indexed family of static critical
maps to nontrivial isolated Euler dynamics. It does not rule out a different
dynamical field map.

### Test B: treating `phi` as the material label used by the same velocity

The natural material-label condition is

```
partial_t phi + dphi(u)=0.                                       (7)
```

It implies frozen-in transport of the pullback area form,

```
(partial_t + L_u)B=0.                                            (8)
```

But the proposed self-identification is `B=i_u nu_g`. For incompressible
`u`, `L_u nu_g=0`, and the commutator identity gives

```
L_u(i_u nu_g)=i_[u,u] nu_g+i_u(L_u nu_g)=0.                      (9)
```

Substitution of (9) into (8) yields `i_(partial_t u)nu_g=0`, hence

```
material transport + B=i_u nu_g ==> partial_t u=0.               (10)
```

Thus the same two target labels cannot simultaneously be nontrivial material
labels and generate their transporting velocity through Slobodeanu's map.
This is the precise dynamical obstruction of the direct R2 extension.

### Test C: action and admissible variations

The actual constant-density Euler kinetic action is the geodesic action on
volume-preserving configurations `eta(t)`:

```
S_E[eta]=(rho/2) integral dt integral_M |partial_t eta(a,t)|^2 da,
eta(t) in SDiff(M).                                               (11)
```

After Euler--Poincare reduction, `ell(u)=(rho/2) integral |u|^2 nu_g`
is varied with the constrained tangent

```
delta u = partial_t xi + [u,xi],   div xi=0,                     (12)
```

which produces (4). Direct substitution of (1) instead gives

```
S_pull[phi]=(rho/2) integral dt integral_M |phi^*omega|^2 nu_g.   (13)
```

This density contains spatial derivatives of `phi` but no `partial_t phi`.
Independent map variations therefore impose a static `sigma2` equation on
each time slice, not (12) or the Euler evolution. Adding
`rho integral P(phi)` produces the paper's static potential equation on each
slice; it does not turn the Bernoulli function into the pressure multiplier
enforcing `eta in SDiff(M)`.

A Lorentzian Faddeev--Skyrme action instead contracts the spacetime two-form
components `F_mu_nu=omega(dphi_mu,dphi_nu)` with a Lorentz metric. Its
`F_0i` terms contain `partial_t phi` and define a hyperbolic field theory.
Those terms, its independent field variations, and its spacetime measure
`dt d^3x` have no equality to (11)--(12) in Proposition 2. Consequently a
Lorentzian completion is an added dynamical law until an action-preserving
reduction theorem is supplied.

## 3. The positive topology bridge that actually survives

For isolated boundary data `phi(infinity)=phi_infinity`, the configuration
space is a component

```
C_Q = Map_*^Q(S^3,S^2),   Q in pi_3(S^2)=Z.                      (14)
```

Krusch and Speight, *Fermionic Quantization of Hopf Solitons*,
arXiv:hep-th/0503067, prove for every component

```
pi_1(C_Q)=Z_2.                                                    (15)
```

The elementary vacuum-component form is the loop-space adjunction
`pi_1 Map_*(S^3,S^2)=pi_4(S^2)=Z_2`; their Hopf-fibration argument transfers
the result to all components. The same work proves that a `2 pi` spatial
rotation loop is noncontractible exactly for odd `Q`. Therefore the topology
permits a two-sheeted universal cover and a nontrivial sign character. With
the **additional Finkelstein--Rubinstein choice**

```
Psi(deck q) = -Psi(q),                                           (16)
```

odd-`Q` states transform spinorially, and localized identical solitons can
receive the associated fermionic exchange sign.

This is a useful, exact bridge: R2 has the correct configuration-space
topology to *permit* spinorial/fermionic quantization, which is stronger than
an integer Hopf label alone. It is not yet a physical Euler bridge, for three
precise reasons:

1. Slobodeanu's physical map is many-to-one. Target area-preserving
   relabelings of `phi` leave `B` and therefore `u` unchanged. One must prove
   that the nontrivial loop in (15) descends nontrivially to the quotient or
   image that represents actual Euler states.
2. An actual localized two-carrier Euler family and its exchange path must be
   embedded in that physical configuration space. The topology of the full
   mapping space does not prove that the exchange path remains within the
   dynamically admissible/stable carrier sector.
3. Equation (16) is a choice of quantum line-bundle holonomy. The classical
   static energy and Euler action do not select the minus character, provide
   a Hilbert-space measure/Hamiltonian, or set the phase normalization
   `exp(i S/hbar)`.

The accepted `C-TOP-001` API is consistent with this boundary: its
`(-1)^w` is only a mathematical character and explicitly does not identify
fermionic statistics. `C-REP-002` supplies finite-dimensional SU(2) and
projector algebra but no Lorentz chirality or physical state. Neither can be
used to fill items 1--3.

## 4. Classical scale and action normalization

Let a dimensionless Euler solution use `y=x/L`, `tau=Ut/L`, and

```
u_phys(x,t)=U u_hat(y,tau),
p_phys/rho=U^2 p_hat(y,tau).                                     (17)
```

For every positive `L,U` this is again an Euler solution. Its natural energy,
time, circulation, and action scales are

```
E_0 = rho U^2 L^3,
T_0 = L/U,
Gamma_0 = U L,
S_0 = E_0 T_0 = rho U L^4 = rho Gamma_0 L^3.                     (18)
```

Accordingly, the dimensionless equality in (S5) becomes a physical kinetic
energy only after multiplication by `E_0`; a history phase needs the further
independent ratio

```
S_0/hbar = rho U L^4/hbar.                                       (19)
```

Neither Proposition 2 nor Euler scaling selects `L`, `U`/circulation, or
`S_0`. The converse construction also chooses a target metric compatible with
the descended area form, leaving an additional normalization choice. Hence a
Hopf integer or flux-helicity unit is not an absolute physical action quantum.

This agrees with the exact scope of accepted `C-ACT-001`, `C-DIM-001`,
`C-DIM-002`, and `C-DIM-003`: those APIs manipulate declared action, speed,
length, and mass coordinates but select none of their physical values.

## 5. Relativistic transfer test

The microscopic equations used in R2 are Galilean and have the continuous
two-scale covariance (17). Their pressure is fixed at each time by an elliptic
constraint. Proposition 2 adds no time coordinate, invariant signal speed,
Lorentz metric, boost action, mass shell, or causal current. Therefore no
coefficient in the stationary map can presently be interpreted as a derived
`c_0`, and no exact implication to

```
E(P)^2 = c_0^2 |P|^2 + m^2 c_0^4                           (20)
```

exists. An emergent relativistic dispersion on a special carrier sector is
not globally ruled out; it requires a derived collective action or principal
symbol with a finite invariant speed and controlled corrections.

Accepted `C-WLN-001`--`C-WLN-003` prove conditional worldline identities only
after a metric, mass, and signal speed are supplied. `C-LOR-001`--`C-LOR-002`
prove conditional orbit/little-algebra statements. `C-SKY-002` proves
pointwise identities for a declared dynamical O(4) Skyrme action. None maps
those inputs to incompressible Euler or licenses (20).

## 6. Input/output ledger

### Derived exactly in this route

- `B=phi^*omega=i_u nu_g`, incompressibility, and the steady Euler equation
  under the source hypotheses.
- Equality of the quartic static density with unit-density kinetic energy.
- Direct time-dependent transfer forces stationarity under either
  instantaneous criticality plus isolated Euler evolution or self-advection
  of the same label area form.
- For global based maps, integer Hopf charge, source-normalized flux helicity,
  `pi_1(C_Q)=Z_2`, and nontrivial `2 pi` rotation for odd `Q`.
- Continuous Euler energy/action scaling (18).

### Still declared/imported rather than derived

- Target surface/area normalization and the global simple quotient in the
  converse direction.
- A specific stable isolated Euler carrier and its two-body exchange path.
- Descent of the `Z_2` loop through the field-to-fluid relabeling quotient.
- The nontrivial FR character, Hilbert space, Born rule, and `hbar`.
- A finite invariant speed, Lorentzian collective action, mass shell, and
  physical electric/weak currents.

No empirical electron or neutrino value is used as an input or comparator in
this attempt.

## 7. Failure classification and materially different executable route

The direct time-dependent extension of R2 has a **representation/action
failure**, not a global impossibility result. It uses the velocity-flux form
`i_u nu_g` where Euler's naturally advected two-form is the vorticity
`d u^flat`; enforcing material transport of the former freezes the flow.

The failure generates route **R2b: Euler vorticity coadjoint-orbit
quantization**, which changes the dynamical object and formalism:

An exact local execution already shows that this representation repairs the
kinematic freeze. On flat space set

```
q=exp(y-c t),
u=(q,c,0),
alpha=q,
beta=x-t q.                                                       (21)
```

Then `div u=0`, `partial_t u+(u dot grad)u=0`, and direct differentiation
gives

```
D_t alpha=D_t beta=0,
d alpha wedge d beta=-q dx wedge dy=d u^flat.                    (22)
```

Thus an advected pullback can represent the actual time-dependent Euler
**vorticity** while velocity is recovered separately. This is a positive
dynamical representation theorem, not just a proposed route. The example has
a constant background and is not a finite-energy isolated carrier, so it
licenses the changed formalism but not P2 or a particle interpretation.

1. Take the actual Euler momentum `m=rho u^flat` and vorticity two-form
   `Omega=d u^flat`. Represent `Omega` locally by Euler/Clebsch labels or
   globally by the full `SDiff(M)` coadjoint orbit; do not identify it with
   `i_u nu_g`.
2. On `R^3`, reconstruct the decaying divergence-free velocity by the
   Biot--Savart/Hodge inverse and use the nonlocal Hamiltonian
   `H[Omega]=(rho/2) integral |u[Omega]|^2`. This retains actual Euler
   evolution. Nonzero-helicity sectors require multichart labels or the full
   orbit rather than one global Clebsch pair.
3. After P2 supplies one actual stable localized carrier, embed its center,
   orientation, internal labels, and two-carrier exchange family in one Euler
   coadjoint orbit. Pull back the Kirillov--Kostant form
   `Omega_KKS(ad^*_xi m,ad^*_zeta m)=<m,[xi,zeta]>` and prove that the
   resulting collective form is nondegenerate modulo exact symmetries.
4. Compute every two-cycle period and test the prequantization condition
   `[Omega_KKS]/(2 pi hbar) in H^2(C,Z)`. This is the exact place where a
   circulation/density/size combination could select or fail to select an
   action quantum; it cannot be replaced by assigning (19).
5. Compute the image of the physical `2 pi` rotation and identical-carrier
   exchange loops after quotienting target relabelings. A nontrivial holonomy
   must be derived or explicitly licensed on this actual carrier space.
6. Derive the collective Hamiltonian `H(P)` and linearized principal symbol
   before comparison with (20). A Lorentzian limit earns only its controlled
   scale range; a Galilean/nonlocal result redirects P4 without erasing the
   topology result.

Attempt 0001 now supplies an exact exposing requirement for steps 2--3. Its
`material-balances.md` constructs a compactly supported smooth Euler initial
velocity whose decaying pressure has an exact quadrupolar `r^-3` tail, and two
states identical on a material tag have different tag acceleration. Therefore
an R2b collective reduction earns closure only by retaining the full ambient
velocity/pressure state or by proving an invariant/conditional elimination of
it. Center, velocity, spin, shape, and finitely many local moments alone are
insufficient on the unrestricted finite-energy Euler class. The full
coadjoint-orbit Hamiltonian plus Hodge reconstruction retains precisely that
nonlocal state; any finite carrier-moduli restriction must reproduce the 0001
pressure-tail test before its symplectic form or exchange phase can be called
physical.

This route is executable once LP2 supplies the same physical carrier family.
It does not assume a Lorentzian Skyrme action and directly tests action
normalization and exchange on Euler's own symplectic phase space.

## 8. Route verdict

`route_verdict: blocked_with_missing_construction`

The exact stationary Euler--strong-coupling Skyrme correspondence and the
configuration-space `Z_2` permission are established at their stated scope.
The explicit equations (21)--(22) additionally establish that the
failure-derived vorticity-pullback representation carries a genuine unsteady
Euler history at local/nonisolated scope.
The direct R2 promotion to a time-dependent physical quantum/relativistic
bridge is blocked by the absent action-preserving dynamical map: the naive
material-label map freezes, independent map variation is static, the FR sign
does not yet descend to actual Euler carrier exchange, and the continuous
Euler scale leaves `S_0/hbar` free. The named next construction is the R2b
same-carrier Euler coadjoint-orbit reduction after LP2, including KKS periods,
physical exchange-loop descent, and derived collective dispersion.

`evidence_scope: EXACT_STATIONARY_AND_TOPOLOGICAL_TRANSFER_WITH_DYNAMICAL_ACTION_GAP`

This verdict is route-scoped. P4 and the electron/neutrino parent remain
active; no finite set of route failures implies exhaustion.
