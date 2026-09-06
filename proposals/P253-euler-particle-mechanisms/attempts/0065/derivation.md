# Topological and nonlocal charge audit

## 1. A point puncture carries flux only by remaining a puncture

`R^3\{0}` retracts onto `S^2`, so

    H^1_dR(R^3\{0})=0,   H^2_dR(R^3\{0})=R.             (1)

For a velocity field, the closed two-form is `alpha=i_u vol`. The generator
with flux `4*pi*q` is represented by

    u_q=q x/|x|^3,    integral_(S_R) alpha=4*pi*q.       (2)

It is divergence free only away from the puncture and obeys

    div u_q=4*pi*q delta_0.                              (3)

Its exterior energy is finite, but its core energy is not:

    (rho/2) integral_(eps<r<R)|u_q|^2 dx
       =2*pi*rho*q^2(1/eps-1/R).                        (4)

This divergence is not special to the radial representative. Fixed flux and
Cauchy--Schwarz on every sphere give

    integral_(S_r)|u|^2 dS >= |integral_(S_r)u dot n dS|^2/(4*pi*r^2)
                            =4*pi*q^2/r^2.               (4a)

Radial integration yields the right side of (4) as a lower bound for every
field with that flux on the punctured annulus; `u_q` saturates it.

Filling the puncture by a smooth source-free incompressible core makes every
bounding-sphere flux zero by the divergence theorem, so it destroys the
nonzero class in (1). A puncture can therefore label signed flux, but its
source, core energy and force law are additional physics. Route A is refuted
as a smooth finite-energy whole-space Euler charge and preserved as an exact
singular/source candidate.

## 2. A compact loop carries circulation, not a point monopole

For a smooth embedded closed loop `K`, `H^1(R^3\K)=R`; a meridian detects the
circulation class. Locally around a thin core of radius `a`,

    u=Gamma/(2*pi*s) e_phi+O(1),                         (5)

so the ideal line limit has logarithmic core energy. A finite-radius vortex
tube makes the energy finite and realizes circulation by smooth vorticity.
However a closed compact tube has `integral omega=0`. Its far velocity is the
impulse dipole already derived in 0064,

    u(x)={3n(I dot n)-I}/(4*pi*r^3)+O(r^-4).             (6)

The nontrivial meridional `H^1` class controls near-core circulation. It does
not become the `H^2` flux of a point complement and does not generate a
`1/r` pair energy. Route B is established as a genuine topological Euler
label and refuted as a scalar Coulomb mechanism.

## 3. Helicity is an invariant, not an Euler Gauss source

On the decaying whole-space Hodge sector,

    H(omega)=integral B[omega] dot omega,                (7)

where `B=curl(-Delta)^-1` is self-adjoint on the paired divergence-free
domain. Hence

    delta H=2 integral B[omega] dot delta omega.         (8)

This exposes both useful and limiting facts. Helicity can encode linkage and
is conserved on smooth Euler evolution, but it is already nonlocal through
the Hodge inverse. Euler variation makes it a coadjoint invariant/Casimir
row; (8) does not produce an independent massless scalar equation or a Gauss
law whose source is `H`. Its value also changes continuously under Euler
similarity unless a separately defined map sector imposes an integer degree.

The obstruction is representation theoretic, not merely local. Let a
translation-invariant linear source map take a rotation scalar `q` to a vector
`f`. For every nonzero Fourier wavevector, rotational covariance about `k`
forces

    f_i(k)=i k_i F(|k|^2) q(k),                          (9)

for every nonzero `k`, even when `F` is a singular nonlocal multiplier on the
declared tempered/finite-energy domain. Therefore

    P^T_ij(k)f_j(k)=0.                                  (10)

For a pseudoscalar, `q k` is an axial longitudinal vector under `O(3)`. If the
accepted body-force target is polar, parity forbids it; if an axial target is
allowed, the transverse projector still kills it. The alternative construction
`epsilon_ijk k_j k_k` vanishes. At an ordinary `k=0` value no nonzero
rotationally invariant vector exists. Contact distributions supported at
`k=0`, and nonlinear or configuration-dependent topological maps, are not
classified by (9)--(10). Thus a scalar or Hopf sign cannot linearly source the
accepted transverse massless mode in the covered homogeneous multiplier class
without an additional orientation, tensor, broken background, or
nontranslation-invariant structure.

Route C is refuted as an automatic helicity-to-electric-charge map. It leaves
helicity and Hopf data available as internal labels for a later construction
that derives a distinct mediator and joint action.

## 4. What an explicitly nonlocal construction costs

There are two minimal ways to evade (9)--(10):

1. Introduce an internal vector `a` and use `f=P^T(a q)`. The response has a
   `1/r` transverse Green function, but the pair law contains
   `a_1 dot G^T(n)a_2` and is anisotropic. Global force compensation and
   background recoil remain necessary.
2. Introduce a scalar field `phi` with

       A_phi=integral dt dx {phi_t^2/(2c^2)-|grad phi|^2/2}
             +sum_a integral dt q_a phi(X_a).           (11)

   Its static equation has the desired scalar Green function and reciprocal
   Coulomb energy. It also adds a new propagating degree of freedom, speed,
   action normalization and point/core renormalization. None follows from the
   incompressible Euler action or the accepted transverse sector.

A nonlinear or configuration-dependent topological transverse map need not be
covered by (9). It remains a live candidate, but must declare its zero mode,
broken symmetry or nonlocal branch/frame and then prove finite energy,
superposition, recoil, causality and reciprocity. Contact distributions at
`k=0` likewise remain unclassified until their physical core/domain is given.
Merely naming an inverse derivative does not accomplish those steps.

Route D is established as an exact assumption ledger: an anisotropic vector
route is compatible with the accepted transverse Green function, while an
isotropic Coulomb route requires a new scalar/longitudinal or equivalent
singular Gauss sector. This does not authorize that foundation change.

## 5. Verdict and continuation

The puncture, loop-cohomology, helicity/Hopf-as-automatic-source and
translation-invariant linear isotropic multiplier routes do not produce a
smooth finite-energy scalar Coulomb charge from bare incompressible Euler.
This is a coverage result for those frozen classes, not for contact zero modes,
nonlinear/configuration-dependent topological maps, or changed foundations.
It preserves three positive objects: signed puncture flux with a class-wide
singular energy lower bound, smooth vortex-loop circulation with dipolar
interaction, and helicity/Hopf internal labels.

The strongest next in-substrate electron route is an autonomous oriented
carrier/background source with exact recoil, followed by a test whether a
same-carrier internal doublet can turn its tensor interaction into the
required observed scalar law without statistical averaging. Nonlinear or
configuration-dependent topological maps and contact zero modes remain a
second in-substrate candidate with explicit domain costs. The materially
different alternative is the separately governed scalar/longitudinal
foundation extension (11), which must extend the exact retained-state/history
construction of 0042 and cannot be called derived incompressible Euler.
