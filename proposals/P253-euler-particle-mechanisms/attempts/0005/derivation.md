# Exact Euler coadjoint-orbit action and the surviving spinorial bridge

This attempt uses the actual constant-density incompressible Euler phase space. It does not
identify a map-space tangent, an accepted abstract algebra, or a Clebsch cover with a particle.
The exact advance over attempt 0003 is that the transported object is now the **vorticity
two-form** `B=d(u^flat)`, while the velocity is recovered nonlocally from `B` and the complete
boundary/harmonic data. Thus time-dependent Euler histories, including their ambient pressure
closure, are retained.

## 1. Convention, measure, and domain map

Let `(M,g,nu)` be an oriented Riemannian three-manifold. On a domain with boundary, all
Euler velocities and variations are tangent to the boundary. On `R^3` they have the decay
needed by the displayed integrals. The density `rho>0` is constant and `G=SDiff(M,nu)`.
Its Lie algebra `g` is the divergence-free vector fields, with the ordinary Jacobi bracket.
The regular dual and pairing are

```
g* = Omega^1(M)/d Omega^0(M),
m = [rho u^flat],
<m,xi> = integral_M rho u^flat(xi) nu.                    (1)
```

Exact one-forms pair to zero by incompressibility and the boundary/decay condition. The factor
`rho` is not optional: `m` is physical momentum one-form density. When `xi` is an infinitesimal
displacement, (1) has dimensions of action. The kinetic Hamiltonian is

```
H(m) = (rho/2) integral_M |u|^2 nu.                       (2)
```

Define the coadjoint operator by

```
<ad*_xi m,eta> = -<m,[xi,eta]>.
```

Then `ad*_xi[m]=[L_xi(rho u^flat)]`. A physical tangent to the right-reduced orbit is
`X_xi(m)=-ad*_xi m`. The KKS form, including its sign and density, is

```
Omega_KKS,m(X_xi,X_eta)
    = -<m,[xi,eta]>
    = -rho integral_M u^flat([xi,eta]) nu.                (3)
```

For `delta H/delta m=u`, one has

```
dH(X_eta) = <u,-ad*_eta m>
           = -<m,[u,eta]>
           = Omega_KKS(X_u,X_eta).
```

Thus `i_(X_H) Omega_KKS=dH` gives

```
partial_t m = -ad*_u m,
partial_t[rho u^flat] + [L_u(rho u^flat)] = 0.            (4)
```

For a representative this is

```
partial_t u^flat + L_u u^flat = -d Pi,
partial_t B + L_u B = 0,       B=d(u^flat),               (5)
```

which is incompressible Euler after absorbing `|u|^2/2` into `Pi`. Equation (5) is an actual
unsteady transfer. It avoids the 0003 freeze because it does **not** set the velocity-flux
two-form `i_u nu` equal to the transported vorticity form `d(u^flat)`.

### Hodge/Biot-Savart closure

On a closed manifold, let `G_2` be the Green operator for the Hodge Laplacian on two-forms.
For exact closed `B`, the co-closed velocity representative is

```
u^flat = delta G_2 B + h,                                 (6)
```

where `h=P_H u^flat` is harmonic **instantaneous reconstruction data** and `P_H` is the
fixed-metric harmonic projector. Its coefficients in a fixed harmonic basis are not asserted to
be constants of motion. Indeed `d delta G_2 B=B`. Hence

```
H[B,h] = (rho/2) integral_M |delta G_2 B+h|^2 nu.         (7)
```

For general nonzero `H^1(M)`, the full Euler--Poincare/orbit-A equation (5), not the local
Clebsch curl row below, supplies the missing harmonic evolution. Projecting (5) gives the exact
coupled reconstruction system

```
partial_t B+L_u B=0,
partial_t h=-P_H(L_u u^flat)=-P_H(i_u B),
u^flat=delta G_2 B+h.                                    (7a)
```

The second equality uses Cartan's formula; the exact part is killed by `P_H`. Equivalently,
Kelvin's theorem preserves `oint_gamma(t) u^flat` on material loops. Those circulations
constrain the combination `delta G_2 B+h`; they do not generally freeze `h` in a fixed metric
basis. Membership in the orbit `m=Ad*_g m_0` carries precisely these global relations. For
example, on a flat torus `B=0` and spatially constant `u=h(t)`, (7a) gives `partial_t h=0`, the
row that is absent if one retains only the curl equation.

On `R^3` with decay, `h=0`; writing `B=i_omega nu` gives the full-space inverse

```
u(x) = (1/(4 pi)) integral_R3
         omega(y) cross (x-y)/|x-y|^3 dy.                 (8)
```

On a bounded multiply connected domain, (6) is replaced by the Hodge-Morrey right inverse
`K_bc B` for the stated normal boundary condition, plus the harmonic field fixed by all
circulations. The Green operator, boundary condition, and circulation vector are therefore part
of the state; `B` alone is not a complete Hamiltonian input. This is exactly where the ambient
pressure/nonlocality exposed by 0001 lives rather than being discarded by a local carrier
closure.

## 2. Two exact actions

### 2.1 Full Euler-Poincare action

For a volume-preserving material map `eta_t`, `u=dot eta o eta^{-1}`, and material measure
`nu(a)`, Hamilton's action is

```
S_EP[eta] = integral_dt (rho/2) integral_M |dot eta_t(a)|^2 nu(a).  (9)
```

Eulerian variations are constrained by

```
delta u = partial_t xi + [u,xi],   div xi=0,
```

with endpoint-vanishing `xi`. Substitution into `delta S_EP` gives

```
delta S_EP = -integral_dt <partial_t m+ad*_u m,xi>,       (10)
```

so (9) yields (4). This is the complete time-dependent variational principle, not the static
strong-coupling Skyrme functional.

### 2.2 First-order action on one physical coadjoint orbit (representation A)

Fix `m_0` and use an inverse-reconstruction variable `g`, writing locally

```
m=Ad*_g m_0,    delta g g^{-1}=-xi,
```

so that `delta m=X_xi=-ad*_xi m`. The one-form
`theta=<m_0,g^{-1}dg>` obeys `d theta=Omega_KKS` on these physical tangents in the
convention (3). Since a first-order action uses a potential
`Theta` with `Omega=-dTheta`, the correctly signed local orbit action is

```
S_A[g] = integral_dt [ -<m_0,g^{-1}dot g>
                       -H(Ad*_g m_0) ].                   (11)
```

Its stationary curves obey (4). A change of local section by the stabilizer changes the first
term by a boundary term only when the stabilizer character is consistently specified. A global
one-form in (11) cannot exist when `Omega_KKS` has a nonzero period. The global classical
object is the closed two-form (3), with local action charts.

## 3. Local/scoped advected labels and exact descent (representation B)

Assume here either `H^1(M)=0` with the stated boundary/decay conditions (in particular the
decaying Euclidean sector), or a specific global cocycle/circulation sector whose extra
variables and admissible variations have already been proved. On a spatial chart use Clebsch
fields

```
u^flat = d theta + alpha d beta,
B = d alpha wedge d beta.                                 (12)
```

For constant density the signed action is

```
S_B = integral_dt [ -rho integral_M(theta_t+alpha beta_t) nu
                    -(rho/2) integral_M |u|^2 nu ].       (13)
```

The `theta`, `alpha`, and `beta` variations give respectively

```
div u=0,
(partial_t+u dot grad) beta=0,
(partial_t+u dot grad) alpha=0.                           (14)
```

Consequently `B` satisfies (5). Locally, taking `d` of the momentum equation loses the pressure
exact form; for general nonzero `H^1(M)` it also forgets the harmonic row (7a). Thus (13) is a
local canonical cover of the time-dependent Euler orbit. The `H^1=0`/decay hypothesis removes
the missing harmonic row, but it does not by itself prove that one regular Clebsch chart covers
every field; global use still requires proved regular chart coverage or multiple/enlarged labels.

Its canonical one-form and symplectic form are

```
Theta_B = -rho integral_M alpha delta beta nu,
Omega_B = -delta Theta_B
        = rho integral_M delta alpha wedge_field delta beta nu.   (15)
```

For a physical displacement `xi`, advected labels vary as
`delta_xi alpha=-xi(alpha)`, `delta_xi beta=-xi(beta)`. Integration by parts gives

```
<m,[xi,eta]>
 = rho integral alpha [xi,eta](beta) nu
 = -rho integral [xi(alpha)eta(beta)-eta(alpha)xi(beta)] nu,

Omega_B(delta_xi,delta_eta) = -<m,[xi,eta]>.              (16)
```

Equation (16) is exactly (3). Hence, on its regular coverage, the Clebsch map

```
(theta,alpha,beta) |--> [rho(d theta+alpha d beta)]       (17)
```

is a symplectic realization, and quotienting its local relabeling kernel produces the physical
KKS pullback in that coverage. This is not an analogy to an accepted abstract API; it is the
physical Euler pairing. It is not a proof that these three bulk labels complete every orbit on
an arbitrary manifold with nonzero `H^1`.

On overlaps, canonical target changes satisfy

```
alpha_j d beta_j-alpha_i d beta_i=dS_ij,
theta_j=theta_i-S_ij.                                    (18)
```

They leave `u^flat` and the local spacetime canonical density invariant up to the required total
derivative. Triple-overlap constants encode transition data, but the displayed bulk action does
not by itself vary those constants or derive (7a). On general `H^1`, the complete action used in
this attempt is therefore `S_EP` or representation A, with `(B,h)` reconstructed by (7a); a
Clebsch completion requires a separately supplied and varied harmonic/circulation row. A single
regular, single-valued Clebsch chart forces

```
u^flat wedge du^flat=d(theta d alpha wedge d beta),       (19)
```

so its helicity integral vanishes on a closed domain or under vanishing boundary flux. Nonzero
helicity therefore requires singularities, multiple charts, or a larger label set. Independently,
even a helicity-zero field may require multiple **orbit** charts: a nonzero KKS period forbids a
global potential on the reduced orbit. The label gauge bundle can be global while no global
section of it exists.

## 4. Exact physical rotation suborbit from the 0001 compact swirl

Take the already constructed admissible Euler initial datum

```
y=x-D, r=|y|, u_n(x)=f(r) n cross y,
```

where `f` is a smooth nonnegative radial bump supported in `r<a`, and `|n|=1`. Define

```
I_1=integral_0^a r^4 f(r) dr,
I_2=integral_0^a r^4 f(r)^2 dr.
```

Direct angular integration gives

```
H = (4 pi rho/3) I_2,
L = integral_R3 y cross (rho u_n) dx
  = (8 pi rho/3) I_1 n = j n,    j>0.                    (20)
```

The genuine global action used for the homogeneous-space statement is the ordinary rigid
`SO(3)` action inside the larger volume-preserving group. It rotates this compact field to
another compact field on the same coadjoint orbit. Its stabilizer **inside the rigid-rotation
subgroup** is exactly `SO(2)` about `n`; hence the physical rotation submanifold is

```
O_rot = SO(3)/SO(2) = S^2.                               (21)
```

Compactly supported fields `xi_a=chi(r) a cross y`, with `chi=1` on a ball containing the
support, remain useful tangent lifts: they agree with rigid rotation wherever `m` is nonzero,
so their KKS bracket pairing equals the rigid pairing. They are not an `so(3)` algebra in the
transition annulus, because

```
[xi_a,xi_b]=-chi(r)^2 r_(a cross b),
```

which is not `xi_(a cross b)` there. Likewise, a time-dependent compact lift chosen to equal the
instantaneous rigid generator on a ball containing the moving support lifts any rigid rotation
path and proves compact-fluid orbit membership. The cutoff objects are therefore tangent/path
lifts, not a substitute compactly supported `SO(3)` subgroup.

For the rotation generator `r_a(y)=a cross y`, the ordinary vector-field bracket is
`[r_a,r_b]=-r_(a cross b)`. Restricting (3) gives

```
Omega_rot(X_a,X_b)
 = -<m,[r_a,r_b]>
 = <m,r_(a cross b)>
 = (a cross b) dot L.                                    (22)
```

Writing `L=j n(theta,phi)` and orienting the sphere outward,

```
Omega_rot = j sin(theta) dtheta wedge dphi,
integral_S2 Omega_rot = 4 pi j
                      = (32 pi^2 rho/3) I_1.              (23)
```

This is a nondegenerate physical KKS two-form with action dimensions and a completely fixed
density/sign convention. The full coadjoint stabilizer is not classified here and is much larger
than the rigid `SO(2)` stabilizer. The submanifold omits translations, profile/amplitude changes,
general vorticity rearrangements, ambient modes, and every non-rigid deformation. Its
Hamiltonian restriction is constant, and it is not asserted to be invariant under the full Euler
Hamiltonian. The compact swirl is smooth finite-energy initial data, not an R1 stable particle.

### Local potentials, period, and rotation phases

North/south action charts are

```
A_N=j(1-cos theta)dphi,
A_S=-j(1+cos theta)dphi,
dA_N=dA_S=Omega_rot,
A_N-A_S=2j dphi.                                         (24)
```

Importing a universal action constant `hbar`, a prequantum line bundle exists only if

```
(1/(2 pi hbar)) integral_S2 Omega_rot = 2j/hbar = N in Z. (25)
```

The overlap transition is `exp(i 2j phi/hbar)`. For the physical loop obtained by rotating an
equatorial orientation through `2pi`, the geometric connection holonomy is

```
exp[(i/hbar)oint A_N] = exp(i 2pi j/hbar) = (-1)^N.       (26)
```

The loop is contractible in `S^2`; (26) is curvature holonomy, not a fundamental-group sign,
and changes for a general loop with its enclosed symplectic area. A rotation about the field's
own axis projects to a constant base loop. Its connection holonomy is one, but the **equivariant
SU(2) lift** of the rotation action acts on the fibre after a physical `2pi` rotation by `(-1)^N`.
These are distinct statements. For odd `N`, geometric quantization of this finite-dimensional
suborbit can carry a half-integer rotation representation, conditional on the imported line
bundle, polarization, and quantum interpretation. A line bundle on this restriction is not yet a
proved extension to the full Euler orbit or to a stable-carrier invariant phase space.

Ritter's review explicitly notes that prequantization alone produces a Hilbert space that is too
large and requires polarization before it is quantization. Equations (25)-(26) therefore do not
derive the Born rule, select `N=1`, produce an electron Hilbert space, or prove exchange
statistics.

## 5. Scaling and why the period is useful but not self-quantizing

For a dimensionless profile `bar u` define

```
u_(U,L)(x)=U bar u((x-D)/L).
```

Then

```
H_(U,L)=rho U^2 L^3 bar H,
T_(U,L)=L/U,
j_(U,L)=rho U L^4 bar j,
integral Omega_rot=4 pi rho U L^4 bar j.                  (27)
```

Thus `H T`, the orbit action, and every KKS period have the same classical scale
`rho U L^4`. At fixed spatial profile, the Euler scaling
`u_c(x,t)=c u(x,c t)` gives `H_c=c^2 H`, `T_c=T/c`, and period/action `c` times the
original. Amplitude scaling moves between coadjoint orbits and continuously changes (25).
Consequently Euler plus its classical scaling does not select `hbar` or a discrete `N`; with an
externally fixed `hbar`, (25) selects only a discrete subset of a continuous classical family.

## 6. No-swirl R1 carrier: exact quotient and narrow zero test

For the centered axisymmetric no-swirl R1 field used in 0002/0007,

```
u=u_r(r,z)e_r+u_z(r,z)e_z,
omega=omega_theta(r,z)e_theta.
```

Therefore `u dot omega=0` pointwise, and
`x cross u=(z u_r-r u_z)e_theta` integrates to zero on every centered axisymmetric tag.
The full angular-momentum integral is zero only when its convergence is separately established.
For a directed ring of nonzero impulse and no extra stabilizer, the physical orientation space is
again `SO(3)/SO(2)=S^2`, not `SO(3)`, and `pi_1(S^2)=0`. The map in the homotopy exact
sequence

```
pi_1(SO(2))=Z --> pi_1(SO(3))=Z_2 --> pi_1(S^2)=0
```

is reduction modulo two, so the framed `Z_2` loop is killed by quotienting the axial stabilizer.
Where the global pairing converges, `j=0` also makes the pure rotational restriction (22) zero;
orientation alone is not a symplectic spin orbit. Nonzero hydrodynamic impulse can still pair
translations with orientations on a larger Euclidean/fluid orbit, so this zero does not annihilate
the full KKS form. It refutes only the direct orientation/FR identification for the no-swirl
carrier. Intrinsic swirl, twist, and other internal collective sectors remain active routes.

## 7. Strongest bridge and the concrete missing physical construction

The following statement is established as an exact classical result:

> A physical incompressible-Euler coadjoint orbit has the correctly normalized nonlocal
> Hamiltonian, KKS action, and local/scoped Clebsch realization. For general nonzero `H^1`,
> the full orbit action carries the harmonic/circulation row. The compact swirl contains an
> explicit `S^2` rotation suborbit whose KKS period is `4 pi` times its physical angular
> momentum. Conditional prequantization gives Chern number `N=2j/hbar` and an odd-`N`
> spinorial `2pi` lift.

That is a useful quantum bridge because it identifies an actual Euler symplectic class and the
classical quantity that would have to be quantized. It is not yet an electron mechanism. The
next physical construction is now specific:

1. Construct or import-and-audit a **nonlinearly stable localized Euler carrier with intrinsic
   swirl/twist and nonzero finite `j`**, and prove that its rotational/internal symplectic sector
   is invariant or adiabatically controlled. The no-swirl R1 supplier cannot fill this slot.
2. Derive a scale-selection mechanism that fixes `rho U L^4` and selects an odd integral class,
   rather than choosing `hbar` or `N=1` after the fact.
3. Supply LP2: a stable two-carrier embedding with a collision/exclusion domain. Only there can
   an exchange loop be tested for survival in the full physical Euler quotient and related to the
   spinorial lift. Neither `pi_1(S^2)=0` nor (26) alone gives Fermi statistics.
4. Quantization still needs a selected polarization/Hilbert space, Hamiltonian spectrum,
   composition rule, an extension from the rotation suborbit to the controlled carrier phase
   space, and physical probability/measurement postulates. These are declared inputs or new
   mechanisms, not consequences of prequantizability.

For relativity, the exact orbit action remains Galilean. A viable emergent sector must derive on
the same stable carrier an energy-momentum law and boost generators satisfying, at its stated
accuracy,

```
H(P)^2=M_*^2 c_*^4+c_*^2 |P|^2,
{K_i,H}=P_i,
{K_i,P_j}=delta_ij H/c_*^2,
{K_i,K_j}=-epsilon_ijk J_k/c_*^2.                         (28)
```

Bare Euler does not supply (28), but microscopic Galilean symmetry is not by itself a no-go for
an emergent restricted sector. The executable test is to pull `H`, `P`, `J`, and the KKS form
back to the eventual stable intrinsic-swirl/LP2 manifold and compute its dispersion and Poisson
algebra. This is a route gap, not a physical impossibility.

## 8. Route verdict

R2b-A (physical orbit) is **established as stated** for the exact Euler Hamiltonian/action,
general-`H^1` harmonic/circulation evolution, and the compact-swirl rotation period. R2b-B
(multichart labels) is **established locally and by gauge quotient** on its regular Clebsch
coverage when `H^1=0`/decay applies or the global sector has been separately proved. The
displayed bulk labels do not establish a complete arbitrary-`H^1` action; such a sector uses
R2b-A or must add and vary its harmonic/circulation variables. The direct no-swirl orientation/FR route is
**refuted with mechanism named**: its stabilizer quotient is `S^2`, its centered spin vanishes,
and its rotational KKS density is zero where defined. The particle implication remains active,
with stable intrinsic angular momentum, scale selection, LP2 exchange descent, quantum
completion, and emergent boost algebra as the named constructions.

## 9. P253/0011 bounded correction receipt

The independent review in `attempts/0011/review.md` audited the pre-correction derivation hash
`0cde70fc4fe089654d42af7095d15380bab26a5f8b4cb8c3d1b81846913ff294` and required two
bounded representation repairs. This revision:

1. restores (7a), the general-`H^1` harmonic/circulation row supplied by the already-established
   Euler--Poincare and representation-A actions, and scopes representation B to local/decaying
   or separately proved global sectors; and
2. identifies the cutoff rotations as tangent/path lifts while using genuine rigid `SO(3)` in
   the larger volume-preserving group for the global `S^2` orbit.

The coadjoint sign/density, Hodge inverse, two full actions, local Clebsch pullback, swirl energy
and angular momentum, KKS form and period, chart jump, conditional prequantum result, scaling,
and no-swirl quotient findings are unchanged. The six recorded symbolic checks address only
those unchanged algebraic quantities and were deliberately not rerun.
