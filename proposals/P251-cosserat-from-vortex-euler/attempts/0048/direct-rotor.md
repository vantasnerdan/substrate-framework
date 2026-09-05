# Direct EPS common and physical relative-angle action

## 1. The actual EPS field and its relative orbit

Choose the actual R³ Beltrami field in the construction of EPS
[1210.6271](https://arxiv.org/abs/1210.6271), Theorems 1.1 and 8.3. In the
proof of Theorem 8.3, pages 56–57, a **finite** spherical-Bessel/harmonic
sum is acted on by `(curl curl+lambda curl)/(2 lambda²)`. This supplies the
knotted-tube field, not just a comparison field. Its derivatives decay
O(1/r). More particularly, its derivative under a rotation about the origin
also decays O(1/r): angular differentiation preserves finite spherical
harmonic sums, and the curl polynomial commutes with physical rotations.
Changing the rotation center adds a constant translation derivative, still
O(1/r). Merely multiplying the theorem's bound on grad u by r would not
prove this stronger rotational bound; the explicit source construction does.

Pick a point x0 on an EPS core vortex line with omega0(x0) nonzero and a
unit axis e perpendicular to it. Let K=e cross (x-x0). The actual rotation
tangent is

    vK=e cross u0-(K.grad)u0,  curl vK=lambda vK.

It is nonzero, since vK(x0)=e cross u0(x0), and O(1/r). It need not belong
to L². Therefore this construction uses the following specified **relative
orbit**, not the finite-energy Euler phase space with a hidden norm change.

Its states are compact volume-preserving vorticity rearrangements of omega0,
followed by a global rotation about e. The corresponding velocities are the
decaying Biot–Savart solutions, equivalently `R_B*(u0+v_c)`, where the
vorticity change of v_c is compact. A compact rearrangement has zero integral
vorticity change (differentiate along its compact generating isotopy), so
`v_c=O(r^-3)`, `curl v_c=0` outside a compact set. The same holds for compact
tangent velocities `v_xi=P(xi cross omega0)`.

Define relative energy by centered rotation-invariant balls:

    Erel(u)=lim_R (rho/2) integral_{|x-x0|<R} (|u|²-|u0|²).

For a pure rotation the integrals cancel **exactly at every R**. For a compact
rearrangement the cross integrand is O(r^-4) and the squared perturbation
O(r^-6), hence absolutely integrable. After rotation the ball is unchanged,
so Erel is exactly independent of the common angle B. This is a fixed
renormalization prescribed by the symmetry, not an adjustable counterterm.

For compact generators use the usual Euler orbit Hessian and KKS form

    H(xi,eta)=rho integral [v_xi.v_eta-v_xi.curl v_eta/lambda],
    Omega(xi,eta)=rho integral omega0.(xi cross eta).

Extend them by one global rotation direction K. `H(K,K)=0` both because the
pure relative energy is exactly constant and because its Beltrami Hessian
integrand vanishes pointwise. For compact eta the mixed integrals are finite.
The curl integration boundary term is O(R^-2), since vK=O(R^-1) and
v_eta=O(R^-3), so

    H(K,eta)=rho integral (vK-curl vK/lambda).v_eta=0.

Likewise Omega(K,eta) is a finite compact-support integral. One global axis
per microscopic realization suffices; no pairing between two divergent
global rotation generators is invoked. The Lie bracket of K with a compact
generator is compact, so the usual KKS closedness/Jacobi calculation has
only these finite pairings. It defines the relative symplectic form on the
declared rotation-plus-compact orbit; its local finite-dimensional action
exists wherever the selected KKS matrix is nonsingular.

The relative angular impulse along e is also fixed by centered balls:

    Jrel=-rho/2 lim_R integral_{|x-x0|<R}
                    |x-x0|² e.(omega-omega0).

Its pure-rotation contribution cancels at every R; its remaining vorticity
change is compact. Integration by parts, using
`curl(|x-x0|² e)=-2 K`, gives the exact moment map

    dJrel(eta)=Omega(K,eta)=:l(eta).

Thus the common angle and compact momentum response are defined on one
Euler relative orbit. No time-dependent rigid velocity `Bdot K` is appended
to u0, and no infinite material locked inertia is replaced by a finite one.
This remains a constrained quadratic Euler-orbit ensemble, as in 0045; it
does not assert that this finite-dimensional family contains every Euler
trajectory or supplies a global finite-L² material displacement mode.

## 2. Compact physical jets, momentum, and positive high-frequency cages

As in 0045, express the finite norm inequalities in one fixed reference
length unit (so the shorthand |k|>=1 is dimensionless). The final action
coefficients are evaluated by the physical integrals, not by assigning that
reference length a constitutive meaning.

Choose a second nearby core point x1, still with e not parallel to its
vorticity direction. In two small disjoint balls define smooth compact
divergence-free jets by

    xi_Rj=curl[-chi_j |x-xj|² e/2],  chi_j=1 near xj,
    Q_R=xi_R0-xi_R1.

Near the tracked core points these are the rigid spatial rotations
`e cross (x-x0)` and `-e cross (x-x1)`. The physical material orientation
jets therefore have angles `theta_plus=B+q`, `theta_minus=B-q`. At each
point the vorticity direction changes by the corresponding cross product
with e; these are observable rotations, not relabeling amplitudes.

Choose a nonnegative smooth cutoff chi in a third small region, away from
both physical jets, with curl vK nonzero, and define

    eta0=curl(chi curl vK),
    l0=l(eta0)=-rho integral chi |curl vK|² != 0.

The displayed identity follows from compact integration by parts and
`l(eta)=-rho integral vK.eta`. Here K cross omega0-vK is a global gradient:
its curl vanishes and R³ is simply connected, so it pairs to zero with a
compact divergence-free eta. No L² Leray projection of the growing global
generator force is presumed. Analyticity and vK not identically zero let
the support be chosen away from the two jets. Choose all regions in a
small neighborhood where a coordinate component omega_z has one positive
lower bound. This carrier axis z need not equal the rotation axis e.

Add a compact negative-helicity circular cage A_k in a body region disjoint
from eta0, and
a circular pair C1_k,C2_k in a fourth region disjoint from the body and both
jets. The exact formula is that of 0045:
`C_i=-curl(phi p_i)/k`, `p1=(cos kz,sin kz,0)`,
`p2=(-sin kz,cos kz,0)`. Choose signed k with k/lambda>0. Put

    r0=eta0+A_k,  q0=Q_R+C1_k,  s0=C2_k.

The body and internal supports are disjoint, so
`Omega(r0,q0)=Omega(r0,s0)=0` **exactly**, independent of projection tails.
The internal B0k=Omega(q0,s0) is positive at a finite threshold by 0045.

Here is a joint finite bound, not an assumption of energy orthogonality.
The three principal projected cage fields have disjoint supports except for
the two orthogonal circular internal polarizations. Hence their leading H
is `rho(1+|k|/|lambda|) diag(A_b,A_i,A_i)`, with both A_b,A_i positive
integrals of `(phi omega_z)²`. The full projection errors obey the D/|k|
and E estimates of 0045, summed over the three cage amplitudes. Their H
remainder is bounded by an explicit constant times rho|t|² by the same
Cauchy–Schwarz expansion. The two fixed compact attachments eta0,Q_R have
finite H and finite norms `||v||+||curl v||/|lambda|`. Moving curl to these
fixed fields bounds every attachment/cage cross by a k-independent constant.
Thus there is a finite C_H, computed by adding those displayed norm bounds,
such that the COMPLETE compact 3 by 3 matrix satisfies

    H_raw >= rho[(1+|k|/|lambda|) A_min-C_H] I,
    A_min=min(A_b,A_i)>0.

Choose a finite |k| above `|lambda| C_H/A_min` and the KKS threshold. All
three compact energy directions are then strictly positive. Disjoint support
has removed only principal energy crosses; the exact projected crosses are
retained in H_raw.

Since K and omega are smooth on every compact support, integrating each
oscillatory moment integral by parts twice gives finite constants L_b,L_s
with

    |l(A_k)|<=L_b/|k|²,  |l(s0)|<=L_s/|k|².

For the 1/k cutoff terms one integration by parts suffices. These constants
are sums of L1 norms of first and second z derivatives of the corresponding
compact amplitudes. Hence a further finite threshold ensures
`|l(r0)|>=|l0|/2>0`. Also l(q0) is uniformly bounded. The selected construction
below improves this to an exactly fixed common moment b=l0.

An implementation refinement suggested by `/root` imposes the moment
conditions using the FIXED eta0, so every high-frequency energy change stays
bounded. This is the selected formula, replacing a projection along r0:

    r=[1-l(A_k)/l0] eta0+A_k,
    Q=q0-[l(q0)/l0] eta0,
    S=s0-[l(s0)/l0] eta0,
    b=l(r)=l0.

The corrections are supported away from the jets. In addition
Omega(r,eta0)=0 by disjoint support of A_k and eta0. Thus they give

    Omega(K,Q)=Omega(K,S)=0,
    Omega(r,Q)=Omega(r,S)=0,
    Omega(Q,S)=B0k=:c != 0.

All attachments are bounded multiples of FIXED compact generators. Thus
the same finite H estimate applies to the selected three generators. For
an explicit choice of its constant, sum uniform L² bounds R0 and R1 for
their attachment velocities and curls; the moment estimates above provide
the bounded coefficients. Put A_max=max(A_b,A_i), M0=sqrt(A_max), and use
summed cage projection constants D,E and a transverse-amplitude derivative
bound G as in 0045. One sufficient constant is

    C0=2 M0 D+D²+(D M0+D G+M0 E+D E)/|lambda|,
    C_H=C0+R0²+R0 R1/|lambda|
        +2(R0+R1/|lambda|)(M0+D).

It bounds the full symmetric matrix remainder in operator norm by rho C_H.
The required uniform attachment coefficients can, for |k|>=1, be bounded
by `1+L_b/|l0|`, `(|l(Q_R)|+L_q)/|l0|`, and `L_s/|l0|`, respectively;
L_q is obtained by the same compact oscillatory moment bound as L_s.
There is consequently an explicit finite threshold making their full H
positive. We have constructed the four independent physical
directions `(K,r,Q,S)` with canonical KKS blocks b and c, a null common
energy direction K, and a positive FULL compact 3 by 3 energy matrix.

## 3. Full action, including the gyroscopic connection

Use physical configuration coordinates `(B,q)` and compact shape coordinates
`(y,s)` along `(r,S)`. Write the positive compact Hessian in this order as

    H_compact = [[P11,P12,g1],
                 [P12,P22,g2],
                 [g1, g2, h]],
    D=diag(b,c),  P>0,  Kq=h-g^T P^-1 g>0.

The exact quadratic orbit action is

    L=(y,s) D (Bdot,qdot)^T
      -[(y,s) P (y,s)^T+2 q g^T(y,s)^T+h q²]/2.

Eliminating BOTH shapes gives, with V=(Bdot,qdot)^T,

    (y,s)^T=P^-1(D V-g q),
    L2=V^T M V/2-q V^T n-Kq q²/2,
    M=D P^-1 D>0,  n=D P^-1 g.

The term `-n2 q qdot` is a total derivative. The term `-n1 q Bdot` is a
genuine gyroscopic coupling and is retained for a single circulation sign.
Dropping it by completing a square would change the physical problem.

Now take the declared time-reversal pair `{u0,-u0}` with equal weights and
the SAME geometric generators and coherently tied physical B/q fields.
The conjugate shape coordinates are INDEPENDENT in the two fluid
realizations: the full ensemble action is

    L_ensemble=[L_plus(B,q,y_plus,s_plus)
                +L_minus(B,q,y_minus,s_minus)]/2.

Vary and eliminate each pair separately, giving
`(y_plus,s_plus)=P^-1(D V-g q)` and
`(y_minus,s_minus)=P^-1(-D V-g q)`. Tying the conjugate momenta as well as
the physical angles would cancel the first-order KKS term before reduction
and would not yield the claimed inertia; that is a different ensemble.

Both are stationary constant-lambda Beltrami fields with identical pressure
and unoriented EPS tube geometry. The compact H is unchanged; Omega, b,
and c reverse sign. Ratios used to define Q,S remain unchanged. Therefore
P,g,h,M,Kq are identical and n reverses sign. Averaging the two **already
reduced** actions leaves

    L_pair=V^T M V/2-Kq q²/2.

The equality of positive and negative circulation weights is the declared
independently testable zero-handedness premise. It is not averaging away a
positive fluctuation or assuming zero microscopic energy: every positive
quadratic coefficient remains intact. The conjugate shapes in the two
realizations take their respective Euler-derived stationary values.

## 4. Physical absolute-angle field and finite positivity of its map

Let `M=[[m00,m01],[m01,m11]]`. The actual cage-section angle is beta=B-q.
Define

    d=m00+m01,  e_mass=m00+2 m01+m11>0,
    a=e_mass/d,
    Psi=beta+a q.

Whenever d is nonzero, Psi has unit weight under a common physical rotation
and the exact action becomes

    T=J_Psi Psidot²/2+J_beta betadot²/2,
    V=K_Psi(Psi-beta)²/2,
    J_Psi=d²/e_mass>0,
    J_beta=det(M)/e_mass>0,
    K_Psi=Kq/a²>0.

The nonzero d is earned by this construction, not inferred from arbitrary
positive definiteness. The fixed-eta0 construction above yields
`|P12|<=rho C_H` and
`P22>=rho[(1+|k|/|lambda|) A_i-C_H]`.
Here b=l0 is fixed and c has the explicit upper bound
`c_max=rho[(1+|lambda|/2) B0+Tstar]` for |k|>=1, using the 0045 KKS
formula with `B0=integral omega_z phi²>0` and its cutoff bound Tstar.
From the exact inverse of P,

    m01/m00=-(c/b) P12/P22.

The additional explicit finite threshold

    |k| > (|lambda|/A_i) C_H [1+2 c_max/|l0|]

ensures `|c P12|<|b| P22/2`, hence `d>m00/2>0`. Thus a is finite and
positive. No target mass or stiffness is used in that selection.

All coefficients are finite same-Euler integral functionals of the selected
actual EPS field, cutoffs, and finite carrier. Their dimensions are the
physical energy/KKS dimensions before any declared cell-volume normalization.
The global common direction has not introduced an assigned body mass.

## Route result

**Established as stated (individual review pending).** On the explicitly
defined EPS relative orbit, one actual common spatial rotation plus compact
core/cage rearrangements gives a positive physical relative-angle action
and a positive absolute-angle inertia from the SAME Euler KKS action. The
time-reversal-paired coherent ensemble removes its computed odd gyroscopic
connection and leaves positive cage gradient inertia as well. Both physical
angle jets and the exact finite integral coefficients are identified.

This closes the common/EPS same-object construction in the relative-orbit,
paired-ensemble scope. The parent still performs its spatial-gradient and
affine translation joining on these exact same fields and coefficients;
the present result does not silently supply those independent calculations.
