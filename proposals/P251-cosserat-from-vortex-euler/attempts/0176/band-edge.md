# Actual fixed-marker initial band edge and its time-dependent phase

## 1. Preserved Euler and physical observation

Use the actual n=8 packet and positive marker mechanism of0164, keeping
the full Kelvin, pressure, tag, spin and displacement-current definitions
of0147/0158. Put E=2n+2, c=sqrt(2)Omega delta_star and P=p/p_star.
At a simple positive node z of F(x)=exp(-x/2)L_n^1(x), choose physical
annuli at x=z+/-h with fixed h>0, widths sufficiently small before the
high-carrier limit. Their weights are positive and have ratio encoded by
r=(w_plus+w_minus)/(w_plus-w_minus)>1. A common positive scaling has
no effect on the phase; it remains in the actual mass and spin.

The exact leading physical numerator, not an eigenfrequency renamed
as an angle, is

    A(P,t)=sum w_i F(x_i P^(3/2)) exp(-i c P x_i t/2),
    gamma=sigma(P)+Im(A_t/A),
    sigma(P)=-2Omega+E c P^(-1/2).                     (1)

The full finite packet and actual material transport supply the higher
pressure/reference/packet errors specified in0164. Their constants may
depend on the fixed h and marker condition number. Physical x_i are
fixed during differentiation; a p-dependent moving ring is not allowed.

## 2. Exact nodal quotient and the constructive initial root

For the leading simple-node Taylor field put S=P^(3/2), d=z(S-1),
w_plus=(r+1)/2, w_minus=(r-1)/2. Removing the common nonzero F'(z),
the two-ring numerator is exactly

    A=e^(-i c P z t/2)[(r d+hS)cos theta
                      -i(d+r hS)sin theta],
    theta=c P h t/2.                                  (2)

At t=0 its signed mean is

    mu(P)=z+h(d+r hS)/(r d+hS),
    mu(1)=z+r h,
    mu_P(1)=(3z/2)(1-r²),
    mu_PP(1)=(9/2)r(r²-1)z²/h+O(1).                 (3)

Thus the finite-h nodal initial-band-edge equation is

    (3z/2)r²-h r-[E+5z/2]=0.                         (4)

It has the explicit positive root

    r_h=[h+sqrt(h²+6z(E+5z/2))]/(3z)>1.

As h->0 it gives r_h²=(5z+2E)/(3z)+O(h), the proposed seed.
At that root gamma_P(1,0)=0 exactly in the nodal problem. Its
r-derivative is c(3z r-h)/2>0, so it is a transverse, not fitted,
zero of a physical response. Also

    partial_P² gamma²
      =9Omega c r(r²-1)z²/h+O(Omega c)+O(c² C_h)>0.    (5)

At small fixed h, the full smooth F changes normalized equation(4)
by O(h). Replacing rings by positive smooth annuli of sufficiently
small width changes it smoothly. The simple zero and strict curvature
therefore persist by the finite-dimensional IFT. Fix these actual
annuli first. High-order Kelvin/pressure preparation and the finite
packet then perturb the normalized root by O(delta C_h) plus its
declared packet/transfer errors. A final small adjustment of the same
positive ratio solves the actual initial physical-clock derivative
gamma_p(p_star,0)=0; no frequency is supplied as an input. The Gaussian
packet's nonzero carrier derivative factor only rescales this equation.

The old reference/spin/G equations are not assumed unchanged. Their
analytic radial rows remain independent on the open positive annular
supports; solve them with their new actual amplitude and beta targets.
The radial-weight root is a further observable equation, while the
reference marking has its own finite moment controls. Keeping the
combined Jacobian in this triangular order licenses the simultaneous
IFT: first the unmarked response and ratio, then the independent small
reference/spin/G marking. Actual pressure corrections are kept in both
blocks. Every finite choice has positive material density and positive
action mass; no zero-width or zero-density limit supplies an inertia.

## 3. Initial zero is not a finite-time band edge

Define R(P)=(d+r hS)/(r d+hS) and
D=cos² theta+R(P)² sin² theta. Equation(2) gives the ACTUAL nodal
clock and amplitude connection at every time where A is nonzero:

    gamma(P,t)=sigma(P)-c P z/2-c P h R(P)/(2D),
    ell(P,t)=c P h (R(P)²-1)sin theta cos theta/(2D).  (6)

At P=1, R=r and R_P=(3z/2)(1-r²)/h. Expanding(6) at the EXACT
finite-h root(4), not only its h->0 approximation, yields

    gamma_P(1,t)
      =-3c³ h²(r²-1)(E+2z)t²/8+O(t⁴).              (7)

The coefficient is strictly negative. It is not a numerical residual
and is not removed by the initial root. Independently

    ell(1,t)=c² h²(r²-1)t/4+O(t³).                   (8)

Consequently the integrated mechanical angular current retains

    H=j(t)theta+
       integral j(gamma_t/gamma+ell)theta+integral e_spin,

after the actual G0=j0 theta0 match. It is not made a configuration-only
current by solving(4). The complete0172 Wronskian variance uses the
actual all-time physical columns, so an initial gamma_P zero cannot
be substituted for its time-dependent D1,D2 rows.

Equations(7),(8) establish a concrete mechanism for this route's
nonautonomy: finite positive material annuli sample different advected
axial phases with signed nodal response. The stronger assertion that
the two-annulus initial root itself yields exact all-window autonomy is
refuted in its leading nodal model. This is not a refutation of the
physical initial band edge, whose construction and positive curvature
are established, or of the parent Euler-to-Cosserat objective.

The full nonlinear F changes bounded h terms, but does not justify
declaring all higher physical pressure/time terms zero. The finite-order
construction can give a controlled improvement of the group error;
it is not an exact autonomy proof. The failure-derived next repair is
registered in moment-repair-contract.md: control signed phase variance
with a third positive annulus and expose its next actual cumulant.
