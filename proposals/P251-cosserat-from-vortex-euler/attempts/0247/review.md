# Independent review — complete periodic history/action/current join

Reviewer: `herdr optical-review pane w3:p3`, a separate Codex process from
the `/root` coordinator, 2026-09-05. I authored and implemented no part of
0241, 0246, or their source dependencies. My earlier 0242 and 0244 reviews
enter only through their frozen supplier licenses; this transaction does not
re-review either supplier. No expected conclusion was supplied.

## Decision and strongest supported statement

**Established as stated at the prepared periodic history/action/current
scope.** The new source correctly joins the reviewed acoustic and optical
suppliers realization by realization on one C016 background before positive
whole-field averaging. It retains the complete action and source cross forms,
uses the actual physical observation map, and proves on every fixed compact
time window that the physical fields satisfy the incompressible isotropic
micropolar equations and periodic constitutive virtual work through
`o(|K|^2)` in the declared finite time-derivative/operator norms.

The same preparation retains the literal material-tag plus continuous-ambient
momentum, pressure reaction, convective transport, intrinsic angular current,
initial current charge and `q_t` memory. Equality is correctly claimed for the
bulk current divergence/periodic virtual-work class, not for arbitrary local
tractions or for one preferred pointwise stress representative.

No false or absent load-bearing step was found inside the frozen 0247
boundary. The minimum repair is `none`. Geometry 0245 remains a distinct
pending input, and the compact Euclidean parent objective remains open.

## Frozen transaction and evidence boundary

The acceptance criteria are exactly `0247/README.md`. Review occurred on
branch `research/pr199-completion` at observed integration head
`b69839e642fe1dff494c25d7e5387cb694041625`. The reviewed source and its new
oracle are pinned below together with the exact input identities used by the
join:

| Artifact | SHA256 |
| --- | --- |
| `0246/README.md` | `df01810215c4e0ac027d8fcfd5659aa128a5588f63cbef8a55c6c0ca72a9b3aa` |
| `0246/joint-current-bridge.md` | `ed140e3ed0a09132defc66bfac709c0de271eef9a1ef05025c61acb5900b3ed5` |
| `0246/verify_current_bridge.py` | `84caf00d61c7a234d81c28eff4eb0e4c8616e363808eb950d104c037cec45dfc` |
| `0246/first.stdout` | `35f2768c7fdb204d8489e877e3c1090bd4a70089366559e957df0a47b51bf0e3` |
| `0246/first.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| `0241/joint-residual.md` | `554e42179f9c24a1b30df3df83919b9e21b1036f5c08cf5686203af386d9a3e0` |
| `0241/verify_joint.py` | `cbdd57503ff7ba057e0fcaf87cf33810ff3df1d3a75e19243e49dc3678e4ebc2` |
| `0241/joint-third.stdout` | `7f1cd2886a7fbcd9ef37937d78653563f07e548834bbdaf93bf42672daea0dcf` |
| `0232/material-torque-and-boundary.md` | `941176309dd6f83052d0c2321d60347bb6a56d57805d0c81716e2130b7fec23e` |
| `0235/current-improvement.md` | `2d1da22649c3292fae5936e861fe53ed536546e1e73b0fb618c8ce4d109b1db7` |
| `0242/review.md` | `2825e6045608f2bf235d9ccca23155150096cbcdc39761e54003adaf682ec251` |
| `0244/review.md` | `e87d8632c97a16eaf172164f96048f4dc050fe3976229c4841c815e60618d88e` |

The 0241 residual and current verifier are used only for the physical branch
substitution, `Q` superpotential and endpoint-memory identities consumed by
0246. The older suppliers and their reviews are unchanged inputs. No verdict
from the still-pending geometry transaction is assumed.

## One common preparation and error map

Both reviewed suppliers consist of finite linear initial-data maps on the
same stationary C016 Euler field, with the same laboratory inputs and
whole-state O(3)/time-reversal action. Taking their sum on each realization is
therefore one actual Euler/Lin preparation; averaging occurs only afterwards.
The optical two-fraction law stays positive, while signed controller
coefficients remain initial-field amplitudes rather than probabilities.
Linearity gives the combined physical history, but does not delete quadratic
action cross terms. Those are included in the complete finite joint jet sent
to the reviewed 0228 normalizer.

With `C=i[K cross]`, the physical observation is

    Y=(U,Phi)=T(A,B)+e,
    T=[[I,-j C/(2rho)],[C/2,I]].

Whole inversion symmetry has the needed content here: an axial-to-polar map
is odd in `K`, as is a polar-to-axial map, so there is no unaccounted even
second-order cross observation. The reviewed suppliers allow any fixed finite
time-derivative list, including the two derivatives required in the residual.
Their leading clock, coefficient, remote-control and initial-map errors can
therefore be placed in one `e,r_A,r_B=o(|K|^2)` operator estimate after the
common nested diagonal. This does not infer history from action matching.

The ordering is noncircular. For each history accuracy, both finite acoustic
parity controllers and finite optical current controllers are chosen first.
Their possibly large low-band costs become fixed constants. The enlarged
moment/current derivative list, observation inverse and all cross forms then
enter one finite exponent `D_N`; the later choices of moment order, remote
integration order, `|K|=h^(D_N+1)` and finally small `h` dominate the cubic,
remote, clock and inverse-map errors. No polynomial cost in the already fixed
acoustic low-frequency width is used.

## Physical state inverse and branch residual

At `K=0`, each transverse acoustic position/rate block has fundamental matrix
`[[1,t],[0,1]]`, while each optical block has
`[[cos(nu t),sin(nu t)/nu],[-nu sin(nu t),cos(nu t)]]`. Both determinants are
one. The spatial observation map contributes determinant

    det T=(1+j|K|^2/(4rho))^2

on the five physical position coordinates, with the longitudinal optical row
unchanged. Hence the full position/rate map has a bounded inverse on every
fixed time window near `K=0`. The controlled `C^1` observation error is a
small finite-dimensional operator perturbation, so invertibility persists.
This licenses an operator statement over all retained physical initial data,
not merely one selected trajectory.

Substitution into the canonical physical equations is exact before taking
orders:

    M0 Y_tt+K2 Y
      =M0 T r+(K2 T-M0 T D)z+M0 e_tt+K2 e.

The 0241 calculation shows that `K2 T-M0 T D` has no coefficient below
spatial order three, including both off-diagonal rows; the longitudinal
optical block is exact. Since `K2` has a zero-order optical mass term, the
source correctly requires `e` itself, not only `C e`, to be `o(|K|^2)`, and
also retains `e_tt`. The supplied estimates meet those requirements. Thus the
physical micropolar residual is `o(|K|^2)` in coefficient operator norm.

The parameter map

    mu=rho a,  alpha=j nu^2/4,
    gamma_T=j(c_T-alpha/rho),  gamma_L=j c_L

is the exact 0227/0241 branch map. Selecting `c_T>alpha/rho` and `c_L>0`
before solving the output inverse gives positive bulk spin curvatures without
turning them into claims about an unprepared spectrum.

## Spin target, current memory and literal balances

The current join uses the full mechanical output supplied by 0242, not an
axial detector identity. For the actual whole-law coefficient `q`, evaluated
with its tag fractions and correlations before substitution, define

    Q_ij=q(t) epsilon_ijk U_t,k,
    S_int=S_full-div Q,
    N_int=N_full-Q_t.

The acoustic first-spin relation and the optical second-order full-spin target
are chosen together so that `S_int=j Phi_t+o(|K|^2)`. In particular, the
optical target includes the actual `div Q` contribution from
`U_opt=-j curl(Phi)/(2rho)`; it is not set to `j Phi_t` before the current
representative is applied. The finite target also retains `q_t`, the lower
endpoint, independently supplied initial `G`, and

    integral_0^t q U_t
      =q(t)U(t)-q(0)U(0)-integral_0^t q_t U.

These are actual current rows and finite-window target functions, not a
time-independent microinertia or a displacement dipole with omitted advective
terms. The chosen derivative list controls differentiation of the
`o(|K|^2)` spin error.

For the literal hybrid currents, 0232 gives

    J_H,t+div T_H=div sigma_H,
    S_full,t+div C_H=div mu_H-ax sigma_H.

Setting `F_full=sigma_H-T_H` and `N_full=mu_H-C_H` retains every transport
term and the continuous ambient pressure reaction. `T_H` is symmetric:
tag terms are `M V tensor V`, ambient terms are `rho chi_0 u tensor u`, and
their linear variations remain symmetric. Therefore
`ax F_full=ax sigma_H`. The `Q` shifts cancel identically between
`S_int,t` and `div N_int`, yielding

    rho U_tt=div F_int,
    j Phi_tt=div N_int-ax F_int+o(|K|^2).

This step neither removes convective transport nor equates the pressure-bond
couple current to a canonical local stress pointwise.

## Periodic constitutive virtual work and inherited action

For the stated isotropic energy, direct differentiation gives the complete
canonical force and couple stresses and

    partial W/partial Phi=ax F_can.

Combining their Euler-Lagrange equations with the physical residual and the
literal balances proves

    P_T div(F_int-F_can)=o(|K|^2),
    div(N_int-N_can)-ax(F_int-F_can)=o(|K|^2).

The longitudinal force divergence is precisely the incompressibility pressure
multiplier. For transverse periodic/compactly supported virtual displacement
`v` and arbitrary virtual rotation `psi`, integration by parts gives the exact
identity

    integral [Delta F:grad v+Delta N:grad psi+ax Delta F.psi]
      =integral [-div Delta F.v+(-div Delta N+ax Delta F).psi]

after cancellation of periodic boundary faces. It is `o(|K|^2)` by the two
displayed divergence estimates. On a cut domain the omitted integration term
is explicitly `(Delta F)^T v+(Delta N)^T psi` on the boundary; it is generally
nonzero and remains part of the 0232/0235 localization and superpotential
bookkeeping. The source therefore earns periodic constitutive virtual-work
equivalence, not equality of arbitrary free-boundary tractions.

The one-dimensional curvature-representative freedom is also correctly
classified. Holding `gamma_T,gamma_L` fixed changes the compatible-gradient
energy by a divergence and the couple stress by a divergence-free
superpotential. It vanishes in the periodic/compact-variation bulk action but
changes a free-boundary functional unless that boundary term is transformed
with it.

Finally, 0228 is used at its actual scope. After all physical sources,
observation rows and current targets are fixed, its same-cell Kelvin controls
match the complete branch phase and Jacobi-energy jets, including intrinsic
and cross terms. The branch forms are `T* M0 T` and the corresponding pulled-
back stiffness. Applying the same invertible physical map to both returns
`M0,K2`; using only an energy match or transforming only one form would not.
The normalizer's remote outputs are included in the common error diagonal.
The time-dependent `Q` memory stays in the exact endpoint/boundary action and
does not become an autonomous bulk coefficient.

## Oracle assessment and findings

The strongest evidence is the analytic join above. The captured 0246 oracle
records seven of seven passing exact checks with process exit zero: both branch
Wronskians, the local angular derivative and complete couple stress, symmetry
of convective transport, and the full virtual-work identity including its
boundary flux. Its two negative controls expose deletion of the
antisymmetric-force torque and deletion of the cut-boundary term. The 0241
captured oracle separately exposes a wrong optical-current normalization and
omission of `q_t` memory while deriving the exact cubic branch mismatch.

Those checks prove their encoded identities; the source-level objective bridge
from actual histories to current divergence is supplied by equations (1)-(8),
not by the tally. No additional exposing check was warranted: the potentially
false joins—state invertibility, convective axial stress, angular
superpotential, boundary term and both-form pullback—are decided by explicit
algebra or already captured negative controls.

I found no source counterexample, missing common preparation, unsupported
factorization of whole-law averages, current-memory deletion, action/history
substitution, or periodic/free-boundary scope switch. No correction check is
needed.

## Precise supplier license and parent status

A later claim transaction may import 0246 as the following supplier and no
stronger one:

> Given the accepted stationary periodic C016 cell, the independently
> reviewed 0242 optical and 0244 hybrid-acoustic supplier licenses, the
> reviewed 0228 complete same-cell form normalizer, and the exact 0232/0235
> current representative, there exists one positive-law sequence of smooth
> finite-energy actual Euler/Lin preparations such that, on every fixed
> compact time window and finite derivative list, the physical hybrid
> displacement `U` and covariance rotation `Phi` have an invertible common
> initial-state chart and satisfy the incompressible isotropic micropolar
> equations through `o(|K|^2)`. In those same physical variables, the full
> inherited phase/Jacobi-energy action is normalized through degree two, the
> literal material-plus-ambient momentum and angular balances retain pressure,
> transport, initial charge and `q_t` memory, and their force/couple currents
> are equivalent to a canonical isotropic representative in periodic bulk
> virtual work through `o(|K|^2)`, with the complete cut-boundary difference
> explicit.

This license does not supply pointwise equality of literal and canonical
stresses, arbitrary free-boundary traction, an unrestricted all-`K` Euler
invariant subspace, acoustic-time uniformity, nonlinear finite-amplitude
stability, the pending 0245 geometry verdict, stationary compact Euclidean
ring/EPS density, an accepted claim, or full P251 completion.

- Verification: `symbolic_verified` for the exact branch, action/current and
  virtual-work algebra, supplemented by the analytic common-preparation and
  fixed-window error construction.
- Review: independently audited and supported at the complete periodic
  history/action/current scope.
- Compatibility: compatible with the accepted C016 fixed-cell sector and the
  frozen supplier licenses; no accepted statement changes here.
- Epistemic: established route evidence inside active P251, not canon.
- Route verdict: established as stated.
- Evidence scope: exact conditional periodic source join with analytic
  prepared asymptotics.
- Minimum repair: none.
- Correction check: not needed.
- Parent verdict: active; the Euclidean geometry/density objective remains
  open.

Signed: `herdr optical-review pane w3:p3`, separate independent non-author
reviewer of 0246, 2026-09-05.
