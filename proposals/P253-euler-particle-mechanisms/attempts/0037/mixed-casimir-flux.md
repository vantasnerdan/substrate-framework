# Exact mixed-Casimir flux supplies the threshold derivative

The missing derivative in `scattering-block-reduction.md` is not a normal-
form assumption.  It is a local conservation law of the exact axisymmetric
Euler equations.  This section derives that law before using any spectral
projection.

## 1. Every mixed Casimir has an axial density and flux

In a frame translating at speed `c`, write

    Psi=psi-c r^2/2,
    xi_t+{Psi,xi}=0,
    zeta_t+{Psi,zeta}=2 xi xi_z/r^4.                 (1)

For an arbitrary smooth label function `D`, define

    G_D'(xi)=xi D(xi),
    q_D=zeta D(xi).                                  (2)

The material equation for `xi` and the derivation rule for the bracket give
the pointwise identity

    partial_t q_D+{Psi,q_D}
        =(2/r^4) partial_z G_D(xi).                  (3)

Moreover, for every smooth `q`,

    r{Psi,q}=partial_z(Psi_r q)-partial_r(Psi_z q).  (4)

Suppose the axis and exterior traces in (4) vanish.  This holds for the
smooth finite-excess perturbations here; the exact exterior is irrotational
and the axis fields have their regular Cartesian orders.  Radial integration
then yields

    partial_t M_D+partial_z F_D=0,                   (5)

    M_D(z)=integral_0^infinity r zeta D(xi) dr,
    F_D(z)=integral_0^infinity
       [Psi_r zeta D(xi)-2 G_D(xi)/r^3] dr.          (6)

The background contribution in `F_D` is subtracted when it is not absolutely
integrable; (5) is unchanged because that contribution is independent of
`z`.  Integrating (5) over the whole axis recovers the exact mixed Casimir
`integral zeta D(xi) dnu`, including the two-end interpretation in `0034`.
Thus (5) is the local form of the actual leaf invariant rather than an added
amplitude constraint.

## 2. Linearization about the exact solitary member

For the exact solitary state `(Psi_s,zeta_s,xi_s)` and perturbation
`(K eta,eta,chi)`, (5) linearizes to

    partial_t delta M_D+partial_z delta F_D=0,        (7)

    delta M_D=integral_0^infinity r
       [eta D(xi_s)+zeta_s D'(xi_s)chi] dr,           (8)

    delta F_D=integral_0^infinity [
       (K eta)_r zeta_s D(xi_s)
       +Psi_s,r {eta D(xi_s)+zeta_s D'(xi_s)chi}
       -2 xi_s D(xi_s)chi/r^3] dr.                  (9)

Equations (7)--(9) are identities for the full Hodge-coupled operator in
`0034`, with no radial-mode truncation.  Flat-label collars cause no extra
term: there `zeta_s=0`, an accessible label variation vanishes, and the
irrotational exterior contributes through `K eta` in the first term of (9).

The tempting choice `D=1` detects the critical mode, but it is not a bounded
coordinate in the physical energy topology.  An axis-concentrating smooth
streamfunction `psi_epsilon=r^2 chi(r/epsilon)g(z)`, with `g` a compact
derivative of zero integral, has energy and column graph norm tending to
zero while `M_1` tends to `2g`.  It lies in the dynamically accessible core
tangent and has zero integrated mixed-Casimir variations.  Thus the earlier
inference from nonzero `M_1` to a bounded oblique projection is withdrawn.
The exact identity (7) survives; only this coordinate choice fails.

Choose instead a nonnegative smooth `D` supported strictly inside a regular
label interval, away from the axis and the flat exterior.  In the column
limit, where `xi=L(r)` and `eta=-Delta_*psi/r^2`, integration by parts gives

    delta M_D=integral (psi_r/r) partial_r[D(L)] dr
       -partial_zz integral psi D(L)/r dr.            (10)

There is no radial boundary term.  On a fixed low-frequency window, (10) is
bounded by the physical streamfunction energy.  Its value on the positive
zero-frequency critical mode is

    m_D,0=integral Phi(r)f_0(r)D(L(r))/(c_0^2 r) dr>0  (11)

when the support of `D` meets the nontrivial transition region.  Hence this
regular-label density is a bounded coordinate on the one-dimensional
zero-frequency critical space.

Equation (7) now gives an exact output axial derivative for this bounded
density.  It does not yet give the physical spectral projection.  At nonzero
frequency the density and its companion flux mix the two column propagation
branches, and the fixed-cutoff coordinate has free off-diagonal terms.  The
required next construction is a bounded change from `(delta M_D,delta F_D)`
to the actual adjoint oscillator pair, followed by the scaled weighted
estimate of its solitary-wave coupling.  Conservation alone does not prove
that the `mu^-1` loss cancels.

## 3. What this closes and what remains

The exact gain is the local mixed-Casimir flux law and a bounded regular-label
coordinate with nonzero critical overlap.  The remaining achievement is to
construct the full adjoint spectral pair, prove its `Q` inverse and exterior
trace estimates uniformly in the physical energy graph norm, and derive the
scaled derivative estimate rather than infer it from (7).  Only those steps
can decide whether the one-pass estimate cancels the critical loss.  This
document does not claim that cancellation, the final operator bound, the
nonlinear bootstrap, or any nonaxisymmetric result.
