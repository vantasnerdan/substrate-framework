# Euler similarity, action selection, and causal propagation

## 1. Exact similarity orbit

Let `(u,p)` solve fixed-density incompressible Euler.  For `A,B>0`,

    u_AB(x,t)=A u(Bx,ABt),       p_AB(x,t)=A^2 p(Bx,ABt).       (1)

Substitution gives a common factor `A^2 B` in both terms of the momentum
equation and preserves `div u=0`.  With `y=Bx`, the exact weights are

| object | weight |
|---|---:|
| length, time | `B^-1`, `(AB)^-1` |
| velocity, vorticity, pressure | `A`, `AB`, `A^2` |
| mass, energy | `B^-3`, `A^2 B^-3` |
| linear momentum, vorticity impulse | `A B^-3`, `A B^-3` |
| angular momentum / physical KKS action | `A B^-4` |
| helicity, circulation | `A^2 B^-2`, `A B^-1` |
| material tag inertia | `B^-5` |

The action weight follows either from angular momentum or directly from the
orbit one-form `rho integral u dot xi dx`: a physical displacement scales as
`B^-1`.  An orientation-preserving dilation preserves knot, link, Hopf, and
degree labels.  It therefore moves dimensionful action continuously inside
one topological class.

For any carrier with nonzero action `J` and any target `J_*>0`, taking `B=1`
and `A=J_*/J` produces the same topology with action `J_*`.  Fixed
circulation only fixes `A/B`; then `J` still scales as `B^-3`.  Fixed helicity
fixes `A/B` up to a square root and has the same residual size freedom.  Two
independently prescribed dimensionful invariants can select one member of a
similarity family, but their numerical values are then initial/boundary data;
Euler supplies no mechanism that chooses the same values for a neighborhood
of initial conditions.

**Route A verdict.** Topology-only and conservation-only universal nonzero
action selection are refuted by the explicit same-topology similarity orbit.
This is route scoped.  A scale-breaking carrier constraint or a new discrete
action law remains a positive candidate.

## 2. Exact compact-source pressure response

Taking divergence of Euler gives

    Delta p=-rho partial_i partial_j(u_i u_j),
    p=rho(-Delta)^-1 partial_i partial_j(u_i u_j).         (2)

For a compact quadratic source, put

    M_ij=integral u_i(y)u_j(y)dy.

Two integrations by parts and the Newton expansion give

    p(x)=rho M_ij partial_i partial_j(1/(4 pi |x|))
          +O(|x|^-4)
        =rho M_ij(3x_i x_j-delta_ij |x|^2)/(4 pi |x|^5)
          +O(|x|^-4).                                    (3)

Choose the compact smooth swirl `u_0=f(s,z)(-y,x,0)`.  Axisymmetry gives
`M=diag(m,m,0)`, `m>0`, hence

    p(x)=rho m(|x|^2-3z^2)/(4 pi |x|^5)+O(|x|^-4).        (4)

On the positive axis,

    p(0,0,z)=-rho m/(2 pi z^3)+O(z^-4),
    partial_z p=3 rho m/(2 pi z^4)+O(z^-5).               (5)

Although `u_0` vanishes there, its exact Euler initial acceleration is
`partial_t u_0=-grad p/rho`, which is nonzero at arbitrarily large `z`.
The remainder cannot cancel the displayed leading coefficient for all large
`z`.  Thus the full velocity/pressure evolution has no strict finite domain
of dependence.  The statement is an initial-response theorem; it does not
claim a uniform long-time expansion.

For generic compact vorticity, the Hodge/Biot--Savart recovery is itself a
global integral with algebraic multipoles.  Special compact velocities can
cancel those velocity multipoles, as the exposing swirl does at `t=0`; the
pressure calculation shows that such cancellation does not restore a strict
Euler cone.

## 3. Restricted bands remain useful

An invariant or accurately controlled carrier band with dispersion
`omega(k)` and bounded `|grad_k omega|` has a finite group-velocity scale for
packets restricted to that band.  This is an effective propagation theorem
only after spectral projection, complement leakage, pressure tails, and the
time interval are bounded.  Equations (3)--(5) show why it cannot be promoted
to an exact full-field Lorentz cone.  A fixed ring harmonic such as the
current `ell` in 0052 is an internal oscillator, not yet a propagation band;
the latter requires a smooth spatial family `Pi(k)`, its dispersion, and a
wavepacket estimate on the physical array or large-radius limit.

For a three-dimensional array whose leading cell coupling decays as `r^-3`,
the absolute shell sum scales as

    integral_R^infinity r^2 r^-3 dr=integral_R^infinity dr/r,  (6)

so even absolute summability fails.  A decay `r^-p` is absolutely summable
for `p>3`, while the first spatial moment normally used in a finite
group-velocity/Lieb--Robinson estimate requires `p>4`.  Thus cancellation
only to `r^-4` is still logarithmically borderline for that stronger bound;
an `r^-5` leading multipole is a concrete neutral-cell target.  Every such
tail remains instantaneous in the full field.

**Route B verdict.** A restricted finite-group-velocity sector is a viable
positive construction.  Bare incompressible Euler cannot supply an exact
full-observable Lorentz cone because its pressure response is elliptic and
instantaneous.  P253 still needs the actual same-carrier band and tail bound.

## 4. Minimal foundation extensions and their ledgers

Three materially different additions survive the exact tests.

1. **Action/circulation law.** Postulate a dimensionful action `h_*` and an
   integral condition such as `2I/h_*=N`.  This selects discrete symplectic
   areas only after `h_*` and the permitted `N` are supplied.  It adds no
   causal speed or probability rule.
2. **Hyperbolic substrate.** Add compressibility or a relativistic hyperbolic
   field with characteristic speed `c_*`, then derive incompressible Euler as
   a limit.  This supplies a causal scale before the singular incompressible
   limit, but it selects no action quantum.
3. **Scale-setting environment.** Add a boundary, lattice, or constitutive
   energy with a length/energy minimum.  It may select a carrier size and an
   effective band.  Discreteness, action quantization, and universal speed
   still depend on the actual added structure.

A general classical hyperbolic action makes the independence exact.  Write

    S_(K,c)[phi]=K integral L(phi,partial_t phi,c grad phi,...) dx dt. (7)

Multiplying by `K>0` leaves the Euler--Lagrange equations and their
characteristic cone unchanged, while multiplying the canonical one-form,
symplectic form, and every action period by `K`.  Conversely changing the
ratios in `L` can change the characteristic speed without selecting the
overall prefactor.  Thus `c` cannot select `K`, and `K` cannot select `c`, in
the classical theory.  A single enlarged microscopic action could contain
both `h_*` and `c_*`, but that is one model carrying two independent
dimensionful inputs, not a derivation of either from bare Euler.  No candidate
above presently supplies measurement probabilities or an exchange character.

The least disruptive coherent extension candidate is compressible
Euler--Korteweg or a superfluid phase--density system: barotropic dynamics can
supply finite sound speed, capillarity a healing length, and a `U(1)` phase
sector.  Its overall action normalization is still independent; its state is
`(rho,m,phase/tag)` rather than the incompressible state, and the exact 0042
retained-history theorem would need a controlled limit or a new derivation.

**Route C verdict.** Explicit extensions can conditionally supply either
missing conjunct.  The minimal current ledger needs independent action-scale
and causal/hyperbolic structure.  Choosing one is a foundational change and
is not made by this attempt.

## 5. Consequence for the particle objective

The exact result is useful because it prevents a classical `CP1`, a KKS
period, or a finite carrier group velocity from being mistaken for the
electron/neutrino bridge.  A viable same-substrate particle route must now
join:

- a persistent localized carrier and physical positive doublet;
- a mechanism fixing its action class, rather than an initial-value label;
- analyzer and first-exit dynamics with retained unresolved state;
- an exchange-loop realization and selected character; and
- a controlled effective propagation band with explicit algebraic leakage,
  or a declared hyperbolic foundation extension.

Nothing here refutes such a joined construction.  It proves that topology and
bare incompressible Euler alone do not provide its action or exact causal
conjuncts.
