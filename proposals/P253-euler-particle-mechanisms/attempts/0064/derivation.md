# Exact localized-Euler charge and Coulomb-scaling audit

## 1. Compact-vorticity identities

Let `omega` be smooth, compactly supported and divergence free. Integration of
complete divergences gives

    integral omega_i=integral partial_j(x_i omega_j)=0,    (1)

and, for `M_ij=integral x_i omega_j`,

    M_ij+M_ji=integral partial_k(x_i x_j omega_k)=0.       (2)

With hydrodynamic impulse

    I_a=(1/2) epsilon_(aij) M_ij,                         (3)

(2) gives `M_ij=epsilon_(aij) I_a`. Expanding the exact whole-space
Biot--Savart potential and then taking its curl yields the first permitted
coefficient

    u(x)={3 n(I dot n)-I}/(4 pi r^3)+O(r^-4).             (4)

There is no `r^-1` or `r^-2` velocity term for a smooth compact-vorticity
carrier. The `r^-3` term is leading when `I!=0`; if it vanishes, decay is
faster. Impulse changes the anisotropic dipole coefficient; it does not act as
a scalar Coulomb monopole.

For two such carriers separated by `d n`, integration by parts in the kinetic
cross term gives

    E_12=rho integral u_1 dot u_2
        =rho I_(1,i) I_(2,j) partial_i partial_j
             [1/(4 pi d)]+O(d^-4)
        =rho{3(I_1 dot n)(I_2 dot n)-I_1 dot I_2}
             /(4 pi d^3)+O(d^-4).                        (5)

When the displayed contraction is nonzero, the translation force is
`-grad_d E_12=O(d^-4)`; cancellations give faster decay. It is orientation dependent
and has neither the `d^-1` energy nor `d^-2` force of a scalar electric charge.
This exact route verdict applies to localized compact-vorticity Biot--Savart
interaction, not to a separately constructed collective mediator.

## 2. Why a smooth whole-space flux monopole is unavailable

The only radial field with a Coulomb-sized `r^-2` amplitude is

    u_q(x)=q x/r^3.                                      (6)

Its flux through every centered sphere is `4 pi q`. A smooth divergence-free
field on all of `R^3` has zero flux through every bounding sphere, so `q=0`.
Distributionally, `div u_q=4 pi q delta_0`; moreover its kinetic energy
diverges at the origin. A puncture, source/sink, or additional continuity
field can carry this flux, but that is a declared foundation change rather
than a smooth incompressible Euler particle.

Exterior circulation around a compact vortex core does not evade (6): it has
zero spherical flux and enters the dipole/multipole field, not a scalar
monopole.

## 3. Pressure cannot supply the missing monopole

For compact velocity, or for localized velocity with `T_ij=u_i u_j`
in `L^1` and the finite first moments/boundary decay needed below, physical
pressure obeys

    p=rho G*partial_i partial_j T_ij,
    G(x)=1/(4 pi |x|).                                   (7)

Under those weighted hypotheses the source and its first moment vanish by two
integrations by parts, and the Taylor remainder is `o(r^-3)`. Moving the
derivatives onto `G` gives the first permitted pressure coefficient

    p(x)=rho partial_i partial_j G(x)
               integral u_i u_j dy+o(r^-3)
        =rho{3 n_i n_j-delta_ij}M^u_ij/(4 pi r^3)
             +o(r^-3).                                   (8)

Thus quadrupole order is the first permitted localized pressure term and its
gradient is first permitted at `r^-4`; the displayed coefficient can vanish,
giving faster decay. The density factor in (8) is physical. The stationary
compact-velocity supplier is stronger: outside its support Euler gives
`grad p=0`, so all exterior pressure multipoles vanish. In the virial identity
this reads

    integral u_i u_k=-(delta_ik/rho) integral p,          (9)

when the exterior constant is subtracted. Equation (9) uses physical pressure;
the density-free form uses kinematic pressure.

## 4. Accepted collective continuum and its precise limitation

C-CST-018 has one massless transverse acoustic sector after the gapped optical
angle is eliminated at zero frequency. On divergence-free displacement its
static principal operator is `-mu Delta`. A prescribed point body force has
the transverse Oseen Green tensor

    G^T_ij(x)=(delta_ij+n_i n_j)/(8 pi mu r),             (10)

so the prepared continuum contains a `1/r` response channel. This does not yet
make a localized Euler carrier electrically charged:

1. the source in (10) is a vector force, not an isotropic signed scalar;
2. an internal compact stress has source `f=div Sigma` and
   `integral f=0`, so its far displacement is a derivative of (10), at most
   `O(r^-2)` before additional moment cancellation;
3. the accepted theorem is prepared finite-window/second-variation response,
   not an autonomous invariant background-plus-defect sector;
4. the optical angle is gapped and therefore cannot by itself give an exact
   Coulomb Green function.

A positive collective-mediator route must construct a carrier functional that
sources a massless channel with opposite signs while preserving total
momentum through background recoil, and must derive reciprocal work/current
from the same action. Fixing a vector direction or inserting an external body
force does not supply isotropic electric charge.

## 5. Translation Ward identity removes the collective force monopole

The missing source cannot be supplied by an arbitrary localized defect while
retaining translation invariance. If the defect energy depends on the
massless displacement only through strain, its first variation has the form

    delta E_def=integral Sigma_ij partial_j(delta U_i),
    f_i=-delta E_def/delta U_i=partial_j Sigma_ij,         (11)

with compact `Sigma`. Hence

    integral f_i=0,       fhat_i(K)=i K_j Sigmahat_ij(K). (12)

The `1/|K|^2` transverse Green multiplier applied to (12) has only a
`1/|K|` singularity, corresponding to `U=O(r^-2)`. A second internal defect
couples through its own stress derivative, so the reciprocal pair energy is
generically `O(d^-3)`, not `O(d^-1)`.

Equations (11)--(12) are the translation Ward identity for the accepted
collective action. A nonzero force monopole requires explicit dependence on
`U`, an external pin/body force, or compensation by an enlarged background
state carrying the opposite total force. The latter can remain a candidate,
but its recoil and infinite/thermodynamic-volume limit are load-bearing parts
of the construction. Moreover the acoustic Green tensor is transverse and
vector-valued; a scalar isotropic sign law still needs a separately derived
source representation.

Global compensation does permit a collection of prescribed vector forces
`F_a` with `sum_a F_a=0` to retain common-shift invariance and exchange an
`O(d^-1)` field. It does not produce an isotropic scalar charge for each
carrier. For a local rotation-covariant signed scalar `q` with no internal
orientation, the general linear Fourier source has

    f_i(K)=i K_i F(|K|^2) q.                              (13)

The massless acoustic response is transverse, and

    P^T_ij(K) f_j(K)=0.                                   (14)

At `K=0`, an algebraic scalar-to-vector monopole `F_i=C_i q` would require
`R C=C` for every `R in SO(3)`, whose only solution is `C=0`. A nonzero source
therefore needs an internal vector or pseudovector and gives the anisotropic
contraction `q_1 q_2 a_(1,i)G^T_ij(n)a_(2,j)/r`. Haar averaging is a
statistical preparation, not a covariant coupling of one carrier. Thus neither
local stress coupling nor globally compensated vector forcing realizes the
required isotropic scalar electric source within the accepted field content.

## 6. Route verdicts and continuation

- Route A is refuted as a Coulomb mechanism by the first permitted exact
  `d^-3` cross-energy and `d^-4` force for compact vorticity, while preserving
  impulse as a real anisotropic interaction; coefficient cancellation gives
  faster decay.
- Route B is refuted for smooth source-free whole-space Euler by flux and
  finite-energy identities. Punctured/source models remain distinct candidates.
- Route C is refuted as a Coulomb mechanism by the exact pressure source
  moments; stationary compact velocity has no exterior pressure interaction.
- Route D is refuted for an isotropic scalar source in the accepted transverse
  collective field: the local Ward identity removes an internal force
  monopole, and the only rotation-covariant scalar-to-vector linear map is
  zero. Globally compensated vector-force pairs can retain a `1/r` response
  but are oriented relational sources, not localized scalar charges.

The continuation ladder leaves three explicit candidates:

1. an autonomous carrier/background action with globally compensated vector
   sources and exact recoil, tested for unavoidable orientation dependence;
2. a massless scalar/phase or compressible longitudinal sector, which is a
   separately authorized foundation change and must rebuild or extend the
   0042 state/observable map; and
3. a nonlocal or topological transverse source, whose locality, finite-energy,
   superposition, and extra-structure costs must be frozen before use.

The compact-vorticity hypothesis is load bearing for the power conclusion.
For a constant vector `a`, the smooth divergence-free field

    u_a(x)=(a cross x)/(1+|x|^2)^(3/2)                  (14a)

has zero spherical flux, belongs to `H^s(R^3)` for every finite `s`, and has
finite kinetic norm `||u_a||_2^2=pi^2|a|^2/2`, but it has an oriented
transverse `r^-2` tail and nonintegrable `r^-3` vorticity. Its Fourier
transform is

    uhat_a(k)=4*pi*i (k cross a) K_1(|k|)/|k|,          (14b)

so two translated copies satisfy

    E_ab(d)=2*pi*rho [a dot b+(a dot n)(b dot n)]/d
             +o(d^-1).                                 (14c)

This is admissible smooth finite-energy local Euler initial data and a genuine
escape from the compact-vorticity power law. It remains oriented and
anisotropic, and supplies no scalar sign, conserved charge/current, Gauss law,
or persistent carrier. Therefore it strengthens the continuation universe
without repairing the electric-charge mechanism.

For candidate 1 the joint quadratic action is

    A[U,X_a]=integral dt{(rho/2)||U_t||^2-(mu/2)||grad U||^2
                         +sum_a F_a dot U(X_a)},          (15)

with `div U=0` and `sum_a F_a=0`. Variation gives the same transverse source
used in (10); eliminating `U` gives the reciprocal pair term
`-F_1 dot G^T(X_1-X_2)F_2` (up to the declared action sign convention), and
its translation derivative is the force. This executes reciprocity for the
prescribed vector-source Green problem and confirms the orientation
dependence. It does not yet execute autonomous carrier recoil: (15) has no
carrier kinetic/KKS term, no constitutive derivation of `F_a` from a carrier,
and a point coupling has ultraviolet-divergent self-energy. The next positive
construction must replace `F_a delta_(X_a)` by a smooth compact form factor,
adjoin the actual carrier action/state variables, derive `F_a` from that same
coupling, and obtain total carrier-plus-medium momentum from translation
Noether symmetry. `sum_a F_a=0` then gives global compensation and a dipolar
far field for the neutral assembly while retaining the finite pairwise `1/d`
cross term. None of these missing carrier rows supplies scalar charge or
transfers the prepared continuum to an autonomous Euler defect.

No accepted velocity, pressure, transverse displacement, or gapped angle
field supplies candidate 2 or 3. No electron, magnetic moment, P5, or parent
conclusion is claimed here.
