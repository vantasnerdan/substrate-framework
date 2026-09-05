# Exact canonical phase pullback with the compact Euler orbit pair

This is the positive same-action joining theorem. It is an exact
restriction of the original quadratic Euler material phase action under
the original phase Cauchy--Born exclusion of nonaffine relaxation. The
stronger free-trajectory reconstruction and its observable correction are
kept separately in 0095. No new invariant-manifold condition is imposed on
the conditional action proved here.

## 1. Canonical material variables and the internal graph

Use the real L2 pairing; complex Bloch formulas use its Hermitian
polarization and real action. Define

    B0 xi = (u0.grad)xi,           C0 xi = (xi.grad)u0,
    Vop xi = P(xi cross omega0),
    A xi = Vop xi-curl(Vop xi)/lambda,
    K(xi,eta) = integral xi.Hess(p0).eta-rho <B0xi,B0eta>.

Here B0 is the RAW material transport operator, skew under the retained
global boundary/ensemble pairing. It is not a projection onto selected
coordinates. The same canonical material Jacobi Hamiltonian is

    H_J(eta,pi) = ||pi||^2/(2rho)-<pi,B0eta>
                                  + integral eta.Hess(p0).eta/2,
    Theta_J = <pi,d eta>,
    pi = rho D_t eta   on an unrestricted material solution.

The fixed-Kelvin internal tangent has canonical graph

    T xi = (xi, rho(A+B0)xi)
         = (xi, rho[Vop xi+C0xi]).

The equality uses `curl(xi cross u0)=B0xi-C0xi`. In the compact sector
of 0085, `Vop xi=xi cross omega0` exactly, so all terms are actual compact
fields. The graph is a subspace of the SAME canonical phase space;
the momentum is not an additional assigned fluid rotor momentum.

### Independent signed symplectic calculation

The canonical symplectic convention is
`Omega_J((xi,pi_xi),(eta,pi_eta))=<xi,pi_eta>-<eta,pi_xi>`.
Consequently its pullback is

    rho integral [xi.Vop eta-eta.Vop xi
                     +xi.C0eta-eta.C0xi]
      = rho integral [2 omega0.(xi cross eta)
                                      -omega0.(xi cross eta)]
      = rho integral omega0.(xi cross eta).

The first equality to `2 omega0.(xi cross eta)` uses solenoidality and
the self-adjoint Leray projector (or the exact compact force identity).
The last C0 difference is a pointwise curl identity. This proves the
KKS sign separately; it does not infer symplecticity from an energy norm.

For a transparent nonzero witness on a periodic cell, take
`omega=(cos z,sin z,0)`, `u=-omega`, `lambda=-1`,
`Q=e_z`, `S=(sin z,-cos z,0)`. Then `VQ=(-sin z,cos z,0)`,
`VS=e_z`, `pi_Q=0`, `pi_S=rho e_z`; both the canonical pairing and KKS
cell-average are exactly rho. This witness tests the signed local/cell
algebra, not the compact support construction or a zero-mean torus
circulation reduction. The latter is supplied by 0085 on its actual
full-space compact sector.

### Hamiltonian calculation

Substitution into H_J cancels its mixed B0 terms and gives

    H_J(Txi) = rho ||Axi||^2/2 + K(xi,xi)/2 = H_orbit(xi,xi)/2.

The final equality is the exact same-field identity
`H_orbit-K=rho||Axi||^2` of 0091. Together with the independent symplectic
calculation it proves that the internal H and KKS used in 0085 really
are the pullback of ONE original material canonical action.

## 2. A joint macro/internal embedding, with every cross term

Let `U=Gx` be a declared divergence-free slow material extension and
`V=Gw` its retained macro momentum-velocity column. G includes the ambient
fluid; it is not a tube-volume-fraction normalization. Write internal
displacement `Ez` for all retained compact core, cage and reaction columns.
Use the explicit joint canonical phase embedding

    eta = U+Ez,
    pi = rho [B0(U+Ez)+V+AEz].

Modulo one displayed integration by parts in the time differential, its
complete one-form is

    Theta = rho<V,dU> + rho<B0U,dU>
          + rho<(A+B0)Ez,E dz>
          + rho<V,E dz> + rho<(A+2B0)Ez,dU>.

Indeed the omitted exact differential is `rho d<B0U,Ez>`; B0 is skew.
The complete Hamiltonian is

    H = rho||V||^2/2 + rho<V,AEz> + H_orbit(Ez,Ez)/2
                                  + K(U,Ez) + K(U,U)/2.

No cross has been replaced by an isolated-cell inverse, and neither
`<V,E dz>` nor `<V,AEz>` is dropped at an affine macro jet. The bare
material shear is K(U,U), whose stationary isotropic average is negative
as computed in `slow-translation-lift.md`. Positive compact orbit cages
tied geometrically to STF strain can repair this SAME H through 0096;
one does not have to borrow the different coadjoint affine modulus.

### Full operator form and exact macro momentum elimination

Set `y=(x,z)`, `T=(G,E)` and `J=(0,AE)`. Then

    Theta = rho <(B0T+J)y+Gw, T dy>,
    H = y^*[T^*KT+rho J^*J]y/2
                 +rho Re<w,G^*Jy>+rho<w,G^*G w>/2.

When `G^*G` is invertible on the declared mean space (or after explicitly
quotienting its redundant coordinates), variation of w gives exactly

    w = (G^*G)^-1 G^*(T ydot-Jy).

Therefore, with `P_G=G(G^*G)^-1G^*`, the reduced Lagrangian is

    L = rho ||P_G(T ydot-Jy)||^2/2
        +rho Re<(B0T+J)y,T ydot>
        -K(Ty,Ty)/2-rho||Jy||^2/2.

This formula retains the entire noncommuting macro/internal mass,
gyroscopic and stiffness structure. It is the input to 0102, before any
reaction elimination, symmetry pairing, or normal-form expansion.

For a constant coherent G, compact divergence freedom gives
`G^*E=G^*AE=G^*B0E=0`. The mean mass is consequently the original
`rho||Udot||^2/2`, and the internal action is precisely the 0085 action.
This recovers the exact Galilean/centroid calculation without appending
mass to an already reduced mean kinetic energy.

## 3. Affine moments and actual observations

At an affine macro jet `U=a+h r`, two different moments occur:

    C(xi) = rho integral r cross xi,
    L(xi) = rho integral r cross Vop xi.

They have different dimensions: C is the moment of a displacement
generator and L is the moment of its induced velocity. For compact
divergence-free v, `integral r cross curl v=2 integral v=0`, so

    rho integral r cross Axi = L(xi).

The rotational part of `<V,E dz>` therefore contains C, whereas that of
`<V,AEz>` contains L. Discarding the first while retaining the second
would produce a false absolute-spin inertia. The six-row independent
moment construction is being developed in 0103; the short modular check
here merely independently corroborates its rank.

For a macro displacement tangent G with macro momentum graph rho B0G,
the cross symplectic row is exactly

    Omega(G,xi) = rho integral G.Vop xi
                      -2rho integral xi.sym(grad G).u0.

Thus a rigid rotational G has the actual L row, translations have zero
row, and STF affine gradients have an additional transport row. This
last row remains in the operator formula above; vanishing of a velocity
STF first moment does not remove it.

The actual Eulerian velocity encoded by the canonical phase point is

    delta u_phase = pi/rho-C0eta
                  = V+curl(U cross u0)+Vop Ez.

Its tag displacement is `delta chi_a=-(U+Ez).grad chi_a=-U.grad chi_a`.
The internal compact fields vanish in the boundary collar, as do their
induced velocities. At the uniform coherent jet their parcel mean and
ambient velocity are identically zero, so tube centroid motion is the
same V common to the ambient phase. At nonuniform jets the exact
centroid/pressure and boundary observation formulas of the companion
material bridge are used rather than identifying V with Udot prematurely.

For example, the actual linear parcel spin observation of a general
phase point is

    delta S_a = rho integral_Da r cross delta u_phase
          +rho integral_boundary(Da) (eta.n) r cross u0,

with the centroid correction zero because the stationary parcel has
zero base mean velocity. The compact internal contribution is L(Ez).
The macro boundary term and within-parcel affine velocity spin remain.

The velocity of the material path represented by eta, in contrast, is
`eta_t+curl(eta cross u0)`. Its difference from the phase velocity is

    R = Udot-V+E zdot-AEz.

At a uniform coherent jet this is compact and has zero centroid moment,
but its angular moment is generally

    rho integral r cross R = C(E zdot)-L(Ez).

This formula is the concrete physical observation which the full moment
normalization/reaction calculation must test. It is not erased by calling
canonical p a material velocity. If a stronger freely reconstructed Euler
interpretation is required, 0095 retains the full complement and current.
For the original conditional phase Cauchy--Born action, the exact
pullback, complete macro mass and all observable/correction rows above
are established without an unrestricted invariant ansatz assumption.
