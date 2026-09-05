# Translation solvability and exact axisymmetric virial identities

The input is0253's equation, with epsilon=1/R and fixed straight profile Phi:

    E_e(phi) = -Delta phi + e/(1+e x) phi_x
              -lambda^2 phi + (1+e x)^2 g(phi) = 0.       (1)

All integrations below are over the full meridional plane, with flat compact
support away from the cylindrical axis when R is finite. Pressure is gauged
to zero outside the support. They are exact analytic identities, not small
numerical residuals. G(s)=integral_0^s g(t)dt.

## 1. The actual first translation row

At e=0, -Delta Phi-lambda^2 Phi+g(Phi)=0. An ordinary differentiable
continuation phi_e=Phi+e h+o(e) would solve

    L h = -f,   L=-Delta-lambda^2+g'(Phi),
    f=Phi_x+2x g(Phi).                                    (2)

Inside the positive region, differentiating the straight equation gives
L Phi_x=0. For a regular continuation class in which the boundary pairing
vanishes, the necessary translation condition is

    0=<Phi_x,f>
      = integral Phi_x^2 - 2 integral G(Phi).              (3)

The second equality follows by integrating 2x partial_x G(Phi). The flat
edge makes that boundary term zero. This is a necessary condition for the
specified regular perturbation class; it does not assert Fredholm theory for
L with its singular g'(Phi), or exclude other finite-radius solutions.

The straight two-dimensional scaling identity is

    integral [G(Phi)-lambda^2 Phi^2/2] = 0.                (4)

It follows by differentiating the compactly supported energy
integral (|grad Phi|^2-lambda^2 Phi^2)/2+G(Phi) under
Phi(x,z)->Phi(x/k,z/k). The gradient integral is invariant in dimension two;
the potential integral scales by k^2. Therefore (3) is exactly

    Q(Phi)=integral Phi_x^2-lambda^2 integral Phi^2=0.      (5)

For radial Phi this is half the poloidal gradient integral minus the axial
flow integral. No arbitrary pressure constant is available to change it:
G(0)=0 is fixed by the exterior pressure gauge. The z-translation pairing
vanishes by parity for a radial profile. Moving its center by O(e) adds a
kernel column to h and does not change (3).

## 2. A genuine profile parameter is available

0253 permits g to be designed from a positive strictly decreasing radial
profile with an exact Bessel core c J0(lambda rho), and an exact flat
exponential edge. It imposes positivity of g only on the sufficiently small
label interval, not on the entire transition. Under this freedom, Q is
not sign constrained.

Keep the core through a fixed rho0 with Phi(rho0)=c0>0. A short transition
of width delta taking most of c0 to a small positive value has, by
Cauchy-Schwarz, radial gradient integral bounded below by a positive
constant times c0^2/delta. Smooth endpoint matching and a short exponential
edge can be arranged without changing that bound. Its L2 integral stays
bounded as delta tends to zero. Consequently Q is positive for sufficiently
short transitions.

Alternatively, match the core derivative through a fixed collar to a broad
decreasing tail Phi_L(rho)=c1 H(rho/L), with 0<c1<c0, on a macroscopic
annulus of radii proportional to L. Choose H positive, strictly decreasing
and smooth there, and flat at its outer edge. The core/collar contributions
are bounded independently of L; the broad tail gradient integral is O(1)
and its L2 integral is bounded below by a positive constant times L^2.
The inner join can have arbitrarily small negative slope followed by the
scaled tail; strict monotonicity imposes no uniform slope lower bound.
Thus Q is negative for sufficiently large finite L. An exponential edge
can be chosen with the same logarithmic exponent as0253; allowing its
positive amplitude and width is an explicitly declared profile parameter.

The space of such monotone positive profiles with a common exact inner core
can be connected continuously in H1 and L2, using smooth monotone splices
and a continuously moving outer boundary. Hence an intermediate profile
has Q=0. This is an existence argument for the scalar compatibility row,
not for a finite-R PDE solution. A proof of the later inverse needs its own
transversality or additional parameter argument; continuity and a sign
change alone do not prove a simple zero.

At a balanced profile, multiplying the straight equation by Phi gives

    integral g(Phi) Phi = -lambda^2 integral Phi^2 < 0.   (6)

The transition must therefore contain negative g. A globally nonnegative
absorption law would preclude this particular large-R regular continuation.
Allowing a signed transition is already part of0253's inverse-profile design;
it preserves its positive small-label logarithmic absorption and zero-g core.

## 3. Finite-radius identities independent of perturbation theory

For a genuine compact axisymmetric solution, use the physical streamfunction
psi and stream swirl F(psi), with

    -Delta* psi = F(psi) F'(psi)-r^2 B'(psi),
    p=B(psi)-[|grad psi|^2+F(psi)^2]/(2r^2).              (7)

Its variational functional is

    J = Er+Ez+EF+EB,
    Er=integral psi_r^2/(2r),  Ez=integral psi_z^2/(2r),
    EF=-integral F(psi)^2/(2r), EB=integral r B(psi).      (8)

Varying J and multiplying its Euler equation by r recovers (7), including
the minus sign in Delta*=partial_rr+partial_zz-r^-1 partial_r.
Axis separation and flat boundary justify both following variations.

Under psi(r,z)->psi(r/k,z), the four terms scale respectively as
(k^-2,1,1,k^2). Under psi(r,z)->psi(r,z/k), they scale as
(k,k^-1,k,k). Stationarity gives exactly

    EB=Er,             Ez=2Er+EF.                       (9)

In physical velocity components these include

    integral u_z^2 dV
       = [integral u_r^2 dV + integral u_theta^2 dV]/2. (10)

This is also the Cartesian steady-Euler stress identity: integrating
partial_j[x_k(u_i u_j+p delta_ij)] gives
integral(u_i u_k+p delta_ik)=0 for compact velocity and exterior-zero
pressure. Its diagonal rows require isotropic integrated kinetic stress.
The cylindrical Jacobian and the average of the azimuthal basis vectors
produce the factor two in (10).

Substituting psi=R Phi into the large-R leading terms of (9) yields (4)
and (5). Thus the translation condition has an independent physical origin;
it is not an artifact of the chosen formal expansion.

## Result and continuation

The exact translation/virial identities and a tunable scalar profile balance
are established as stated. The regular continuation of a fixed unbalanced
profile fails this necessary row; the candidate family remains viable through
signed transition-profile design. Next achievement: realize the balanced
profile in a finite-R degenerate free-boundary inverse, with its actual
translation gauge, kernel/cokernel, physical smoothness and nonzero twist.
No statement here supplies that inverse, an arbitrary finite-R exclusion,
a positive density supplier, or the full Euler response/action/current join.

## 4. Frozen0255 convex-family refinement

The more specific0255 family can retain one common outer radius and exact
outer profile. First choose a sufficiently large finite radius and a broad
negative-Q profile, arranging its final amplitude to be small before the
fixed exponential edge. At that same radius choose a steep early cutoff,
followed by a tiny strictly decreasing positive tail to the identical outer
neighborhood. Keeping that outer amplitude small leaves room for a tail
above it. Making the early transition thin produces positive Q while the
rest of the L2 integral remains bounded at this fixed radius. Smooth
monotone joining preserves the exact common inner and outer neighborhoods.

Write Phi_eta=Phi_0+eta H, H=Phi_1-Phi_0. Every convex combination is positive
inside, strictly decreasing away from the center, and has precisely those
common neighborhoods. Both endpoint signs are strict. Moreover

    Q(Phi_eta)=Q(Phi_0)
       +2 eta integral (Phi_0,x H_x-lambda^2 Phi_0 H)
       +eta^2 integral (H_x^2-lambda^2 H^2).             (11)

This polynomial has degree at most two. Opposite endpoint signs force an
interior root at which its derivative is nonzero: a multiple quadratic root
has even multiplicity and cannot change sign. Thus this declared family
supplies the simple scalar border, strengthening the arbitrary-continuous-
path observation in section2. It still supplies no infinite-dimensional
inverse or boundary regularity theorem. The support and equation are fixed
within each finite family; their initial large choice precedes eta selection.
