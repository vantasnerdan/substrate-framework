# Exact quantum-bridge calculations and their physical boundary

## 1. Route A: what the Koopman state does and does not encode

Let `Gamma` be an invariant classical Euler-carrier phase space, let `mu` be
an Euler-invariant probability measure, and let `Phi_t` be its flow. On

    H_K=L2(Gamma,mu),
    (U_t Psi)(x)=Psi(Phi_(-t)(x)),                      (1)

`U_t` is exactly linear and unitary. A real physical classical observable
`f` acts as `(M_f Psi)(x)=f(x)Psi(x)`, so

    [M_f,M_g]=0,
    <Psi,M_f Psi>=integral_Gamma f |Psi|^2 dmu.        (2)

This is a complete positive classical ensemble representation. Multiplying
`Psi` by an arbitrary state-dependent phase leaves every expectation in (2)
unchanged. Differential generators on `H_K` can fail to commute with
multiplication operators, but that algebraic fact does not turn them into
physically implemented Euler measurement operations.

Route A is therefore **established as a classical Hilbert/state/dynamics
representation and refuted as an automatic P4 quantum completion**. The
mechanism is explicit: the physical observable algebra inherited from
classical functions is commutative, the Hilbert phase has no declared Euler
observable, and composition remains classical product probability. A Born
rule for noncommuting alternatives is not derived by rewriting a density as
`|Psi|^2`.

## 2. Route B1: exact sphere quantization and continuous action scaling

For the physical Euler rotation suborbit from `0005`,

    Omega_j=j sin(theta)dtheta wedge dphi,
    integral_(S2) Omega_j=4*pi*j.                      (3)

Given a universal positive action `hbar`, a prequantum line exists exactly
when

    N=(1/(2*pi*hbar)) integral Omega_j=2*j/hbar in Z. (4)

For `N>=0`, the invariant Kähler polarization on `S2=CP1` yields the
holomorphic section space of `O(N)`, of complex dimension `N+1`, carrying the
spin-`N/2` irreducible `SU(2)` representation. Hence `N=1` gives an exact
two-state spin-one-half representation and a physical `2*pi` rotation acts
as `-1` in its lifted representation.

This is a stronger result than prequantization alone, but it does not select
the sector. Under the exact Euler rescaling

    u_(A,B)(x,t)=A u(Bx,A B t),                        (5)

with fixed material density, the angular momentum/action scale is

    j_(A,B)=A B^(-4) j.                               (6)

The dimensionless carrier topology is unchanged while (6) ranges
continuously. With externally fixed `hbar`, (4) selects discrete members of a
continuous Euler family; setting `hbar=2j` merely defines units so that
`N=1`. No Euler equation, conserved topology, or variational minimum in the
current inputs fixes this equality.

Route B1 is **established conditionally as a polarized Euler-orbit
quantization and blocked as physical action/class selection**. The missing
construction is a same-substrate mechanism fixing `rho U L^4` and choosing
the odd minimal integral class without using the desired answer as input.

## 3. Route B2: two-carrier exchange and the character-selection problem

For two noncoincident labeled centers in `R3`, remove the center of mass and
write the relative vector as `r in R3\{0}`. Quotienting by label exchange
identifies `r~-r`; the unordered center configuration retracts to `RP2`, so

    pi_1(C_2(R3))=Z2.                                 (7)

Its flat complex line systems have exactly two unitary characters,

    chi_+(generator)=+1,  chi_-(generator)=-1.        (8)

The center topology permits both bosonic and fermionic exchange. A product of
two identical internal line bundles also admits two swap equivariant
structures differing by (8). Neither the one-carrier Chern number (4) nor the
classical Euler equations choose between them.

The Finkelstein--Rubinstein route is more structured. If a persistent Euler
carrier supplied a based map `n:S3->S2` of Hopf charge `Q`, if its accessible
phase space embedded into the corresponding full map component with the
nontrivial fundamental-group class preserved, and if the odd deck character
were chosen, then the source theorem gives

    phase(2*pi rotation)=phase(exchange)=(-1)^Q.      (9)

Thus odd `Q` supports a consistent spinorial and fermionic quantization. But
each italicized physical bridge is open for Euler, and the source itself
chooses the deck-odd wavefunctions. Equation (9) relates rotation and exchange
after quantization; it does not derive why nature selects `chi_-`.

Route B2 is **established as the exact center/Hopf topology and conditional
exchange descent, but blocked as Euler-selected fermionic statistics**. Its
next positive construction is an LP2/LP3 invariant noncollision carrier
sector with a proved map into the full Hopf configuration component, followed
by a dynamical or locality principle selecting the deck character.

## 4. Route C: positive-frequency linear modes

Let a stable linearized carrier sector have real symplectic form `Omega`, a
positive conserved quadratic energy `g`, and an invertible generator `A` that
is skew-adjoint in `g`. When the polar decomposition defines

    J=A(-A^2)^(-1/2),  J^2=-1,                         (10)

the compatible complex inner product constructed from `g,Omega,J` gives an
exact one-particle Hilbert completion and the classical linear flow is
unitary. This is the precise positive-frequency structure that a stable
carrier could supply dynamically. The hyperbolic returned sector established
in `0032/0038/0039` does not satisfy this premise and is an adverse carrier
selection result.

Canonical quantization of (10) still adds the CCR representation and an
action unit; its Fock functor is bosonic. A CAR/fermionic Fock space does not
follow from a real classical symplectic wave phase space. Nor does the
classical norm derive projective measurement probabilities. Route C is
**established as a conditional dynamically selected complex one-particle
space and refuted as an automatic fermion/Born construction**. A positive
carrier metric remains useful, but P4 additionally needs Route B2's topology
and an action/measurement selection principle.

## 5. Route D: exact Galilean boundary and effective continuation

Incompressible Euler is invariant under the Galilean transformation

    u_U(x,t)=u(x-U t,t)+U,
    p_U(x,t)=p(x-U t,t).                               (11)

Around a quiescent state, Leray projection gives `delta u_t=0`; around a
uniform background `U`, transverse Fourier modes obey

    omega(k)=U dot k                                  (12)

for both physical polarizations. There is no pair of bare acoustic branches
`omega^2=c_*^2|k|^2`, and the same-field pressure remains an elliptic
constraint. Therefore bare constant-density incompressible Euler does not
carry an exact finite-speed Lorentz representation on its full state space.

This does not exclude an emergent restricted band on a structured Euler
background. Such a route must derive a two-sided conical band, a common
finite speed, the collective boost brackets, and controlled errors on the
same carrier sector. It necessarily declares the background/rest frame and
its regime. Route D is **refuted for an exact bare-Euler Lorentz bridge and
blocked positively on a same-substrate effective-band construction**.

## 6. Joined verdict

Euler now supplies an exact classical Koopman representation, a physical KKS
action, a conditionally polarized finite-dimensional spin representation, and
topologies on which fermionic exchange is mathematically consistent. These
are more than analogy. They also expose the irreducible missing physics:

- continuous Euler scaling does not select `hbar` or the `N=1` orbit;
- topology permits two exchange characters and does not select the odd one;
- a Koopman norm is classical probability unless a physical noncommutative
  observable and measurement rule is constructed;
- a positive classical wave sector canonically leads to CCR/bosonic
  quantization unless joined to additional topology; and
- bare incompressible Euler has Galilean convective branches, not a finite
  invariant speed.

Consequently LP4 is not earned by these routes. The strongest failure-derived
positive candidate is now sharply typed: construct a persistent odd-Hopf
Euler carrier whose accessible phase-space inclusion preserves the `Z2`
rotation/exchange loop; find a same-substrate action-selection principle; and
derive a positive two-sided band whose complex structure, observable algebra,
and effective boost relations live on that identical carrier. Any proposed
additional stochastic, topological, or microscopic rule must be named as a
new substrate hypothesis and tested rather than hidden inside the word
quantization.
