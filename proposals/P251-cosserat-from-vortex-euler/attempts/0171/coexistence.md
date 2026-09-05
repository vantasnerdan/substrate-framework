# One fixed generalized-force-free field: acoustic and optical histories

## 1. Object and ordering

Fix the triangular reciprocal lattice, Psi, lambda, rho and an axial
reference period L_s, all positive. Use exactly the field of0163:

    psi=Psi sum_{j=1}^3 cos(b_j.x), v=J grad psi,
    W=F(psi)=sqrt(C+lambda² psi²), u=(v,W).

Choose C sufficiently large, but FINITE, for the positive fixed-profile
acoustic theorem in0163. Keep this C fixed in every limit below. The
pressure per unit mass is -(v²+lambda² psi²)/2; curl u=f u with
f=-lambda² psi/W. Its smooth periodic norms are finite fixed constants.
In particular no constant-curl inverse or constant-factor EPS theorem
is an import for this field.

Let Omega=3Psi lambda²/2 and psi0=3Psi. At the core maximum,

    v=-Omega Jx+O(r³), zeta0=-2Omega,
    W0=sqrt(C+9lambda²Psi²),
    a=F'(psi0)=3lambda²Psi/W0>0,
    b=F''(psi0)=lambda² C/W0³>0.                    (1)

The positive optical carrier is k=p=2pi N/L_s, N a sufficiently large
positive integer. Define

    ell⁴=2/(a p³), delta²=2a/p, cD=Omega delta.       (2)

Thus the acoustic limit k_ac ->0 and the optical limit p ->infinity
are compatible on the SAME fixed field. One may use axial Bloch fibers
or repeated axial cells with commensurate carriers. Neither limit
changes C or the normal periodic geometry. The acoustic assertion is
on its acoustic interval; the optical assertion is on fixed
0<=Omega t<=T. Their common observation interval is the latter, unless
a separate long-optical-time estimate is supplied.

The result below is positive coexistence of actual prepared physical
histories and their phase actions, not an autonomous coupled isotropic
constitutive law or a knotted Euclidean ambient realization.

## 2. The generalized profile enters only after the needed orders

The exact triangular Taylor jet is

    psi=psi0-Omega r²/2+Omega lambda² r⁴/32
        -Psi lambda⁶ r⁶(10+cos(6theta))/7680+O(r⁸).

Consequently

    W-W0=-a Omega r²/2
         +(a Omega lambda²/32+b Omega²/8)r⁴+O(r⁶). (3)

The coefficient involving b is retained; it is not set to zero by a
constant-curl assumption. From(2),

    p=2a/delta², ell²=delta³/(2a²),
    lambda² ell²=O(delta³), p ell⁴=O(delta⁴).

At fixed C the cubic radial swirl first enters the scaled velocity
operator at delta³, the ENTIRE axial quartic at delta⁴, the nonradial
normal C6 term at delta⁶ and its axial Doppler term at delta⁷. Constants
can depend on the chosen a,b,lambda,Psi; these are fixed before p.
Small a is not silently a uniform estimate in C.

For comparison only, the local affine axial profile
W_lin=W0+a(psi-psi0) has exactly the same terms through delta². It is
the axial reflection, up to a uniform Galilean velocity, of the local
profile with slope -a. Axial reflection reverses k, horizontal
vorticity and axial displacement together; it leaves axial physical
spin and the scalar KKS pairing invariant. This is an exact Euclidean
covariance statement, not a constant-curl replacement of lambda by a.

## 3. Full pressure, Kelvin preparation and the resonant oscillator

With the uniform laboratory phase p W0 restored when observing labels,
the normal velocity V, axial velocity B and pressure P obey

    D V+(V.grad)v+grad P=0,
    D B+V.grad W+ipP=0, div V+ipB=0,
    D=partial_t+v.grad+ip(W-W0).                      (4)

Eliminate B=i div V/p by the divergence of the first equation and the
second equation. At the quadratic core this gives, without discarding
axial shear,

    (p²-Delta)P=2Omega curl V+2ip grad W dot V.        (5)

The inverse in(5) is on the FULL normal triangular torus. Expanding
its exact massive resolvent on scaled Gaussian sources gives

    L=L0+delta L1+delta² L2+O(delta³),
    L0=partial_theta+J,
    L1 V=i R² V/2-grad curl V,
    L2 V=i grad(R dot V)-(1/2)grad Delta curl V.       (6)

For example the axial-shear pressure contribution has precisely the
same sign as the trapping Doppler contribution, since p grad W is
-p a Omega r at principal order. A choice k=-p gives the inverted
quadratic potential. The two coefficients balanced in(6) are
Omega/(p² ell²)=p a Omega ell²/2=Omega delta/2.

Put e_+=(1,i)/sqrt(2), m=2,n=2 and

    F_R=R exp(-R²/2)L_2^1(R²) exp(i theta).

The resonant scalar equation is
(-Delta_R+R²)F_R=12F_R; the relative Eulerian frequency is
-6Omega delta+O(Omega delta²). The first nonresonant column is
delta e_- D_+²F_R/8. At the next order the nonresonant L0 equations
have nonzero integer denominators. Circular angular grading removes
the opposite resonance; the pressure cancellation of the two retained
terms in L2 leaves the plus-l radial resonant modulation, not an
asserted isolated eigenvalue. This explicit construction is the same
normal differential equation as0166 through delta²; (3)--(6) supply
its new-field license and its first omitted terms.

Prepare a divergence-free displacement, with axial completion
xi_z=i div xi_h/p, by the finite-order Kelvin expansion. At its
principal order

    v_pert=P(xi cross omega0), xi_h=i V/(2Omega).       (7)

Indeed xi cross (-2Omega e_z)=2Omega Jxi and
J e_+=-i e_+. Higher coefficients come from this invertible circular
principal block and the pressure return in(5). Choose the exact
initial velocity as P(xi cross omega0), not a truncated velocity that
only approximately satisfies Kelvin. Transport it by full Euler and
reconstruct xi by full Lin. This only requires a finite-order local
preparation: it does not invert curl on the variable-factor background
or presume a finite-dimensional invariant Euler manifold. Cutting the
mode off at fixed physical radius makes the preparation error
exponentially small; any prescribed finite number of formal orders
can be retained before this cutoff.

## 4. Complete action and literal material observations

For velocity amplitude V_*, the exact phase KKS number is the complete
cell integral

    beta=rho integral omega0 dot(Re xi cross Im xi).

Its leading density from(7) is -rho |V|²/(4Omega), so

    beta=-3pi rho V_*² L_s ell²/(4Omega)[1+O(delta)]<0. (8)

Here integral |F_R|² d²R=3pi. The horizontal vorticity generated by
F' is O(Omega a ell), and xi_z/xi_h=O(1/(p ell)); their KKS product is
O(delta²). It is retained in exact beta, not inferred to vanish.
The generalized b term enters still later. All pressure images and
exterior-to-core fluid enter this action integral. A tag's volume is
not substituted for L_s or for the action domain.

The tag is a smooth nonnegative label density supported in a single
elliptic tube, with an axial label cutoff of width c/p_star and a
signed radial marking b_tag(x) cos(2theta+p_star s), x=r²/ell_star².
The plus sign is the axial reflection of the negative-carrier label;
p_star and the actual physical tag are fixed when differentiating a
nearby carrier. A small modulation keeps the total density nonnegative.
Observe the ORDINARY centered Euclidean quadrupole

    Q=rho integral_labels [(x-X_D)+i(y-Y_D)]² w0,
    theta_obs=(1/2) delta arg Q.

The response to a rigid physical rotation is one. Axial transport
s(t)=s0+integral W(r(t))dt is used inside the carrier phase. Therefore
the uniform pW0 contribution cancels by actual material travel, not
by choosing a rotating physical clock. The label's normal reference
turns at -2Omega at principal order.

For the untruncated reference radial weight chi=x, the exact Gaussian
Laguerre Laplace ratio is

    integral x³ e^-x/2 L_2^1(x) /
    (2 integral x² e^-x/2 L_2^1(x))=19/3.

Subtracting the Eulerian oscillator value 6 gives the observed phase

    gamma=2Omega+Omega delta/3+O(Omega delta²),
    p² partial_p² gamma²=Omega² delta+O(Omega² delta²)>0. (9)

These are derivatives with the SAME finite physical tag fixed. To
realize its ideal coefficient, use the three-bump construction of0166
on the densities x D0^j[(x/2-19/3)e^-x/2 L_2^1(x)], j=0,1,2,
D0=(3/2)x partial_x. Their Gaussian-stripped degrees are 3,4,5;
three interior bumps give an invertible matrix. The tail is
exponentially small, so a sufficiently distant finite cutoff followed
by this exact three-row correction remains nonnegative. A single
value fit would not license the two derivatives in(9).

Full centered material spin, rather than Eulerian velocity spin, is

    S=rho integral_labels [xi cross u0
             +r cross(v_pert+(xi.grad)u0)]_z w0,       (10)

with centroid terms retained. On the radial comparison it is
rho integral r[v_theta+(2O+rO')xi_r]w0, O=-Omega+O(r²).
The leading expression cancels POINTWISE when(7) is used. The first
spin is the pressure torque at order delta. Its radial polynomial is

    P_2=L_2^1-2(L_2^1)'=x²/2-5x+9.                   (11)

The plus-l order-delta² resonant velocity remainder has this same
zero principal spin, hence spin error order delta³. The additional
generalized-profile term is order delta⁴ in velocity and cannot
spoil the relative order-delta² spin estimate.

The six pressure rows P_2,DP_2,D²P_2,xP_2,D(xP_2),D²(xP_2), with
D=c0+(3/2)x partial_x acting also on e^-x/2, together with reference
rows 1,x have a nonzero eight-row minor. The attached verifier derives
that minor, not a determinant for a different mode. The exact moment
map, including full beta, Kelvin preparation, axial filter and
second-order pressure, is an O(delta) perturbation after dividing by
its explicit spin/volume scales. The finite-dimensional implicit
function theorem therefore supplies a marking with

    S=eta Pi_theta + O_scale(delta²), eta=1 or eta=1/2,

uniformly on the fixed optical interval and through two relative
carrier derivatives. eta=1/2 is a prescribed row normalization, not
by itself a constructed standing-carrier superposition. The reference
moment is O(delta/(p_star L_s))=O(delta³). The reference x constraint
cancels the cubic radial clock term; the next normal reference error
is O(delta⁶), hence O(delta³) after this division. The axial nonlinear
profile does not alter the unperturbed normal reference geometry.

For exact observed amplitude c_theta and phase gamma, the pulled-back
moving phase action is

    M=-beta/(gamma |c_theta|²)>0,
    L=M/2[(theta_dot-partial_t log|c_theta| theta)²
                                                   -gamma² theta²].

All time and carrier connections remain. In particular (9) is a
positive physical optical response, not permission to erase the
connections or call an unobserved eigenfrequency the tag clock.

## 5. Uniform full-pressure estimate at fixed C

Let E=psi0-psi. It is an exact invariant of v and of the full u.
The weight exp(sigma sqrt(2E/Omega)/ell), flattened smoothly inside
the core, has zero transport derivative and logarithmic gradient at
most C/ell. Massive periodic pressure images decay like exp(-p d).
Weight conjugation leaves exp[-(p-C/ell)d]; p ell ->infinity.
For grad(p²-Delta)^-1 div use its unweighted L2 bound first. In the
commutator, the weight ratio minus one removes one degree of the
near-diagonal order-zero singularity, giving O(1/(p ell)); the local
delta term commutes. Thus the Schur argument applies to the integrable
commutator, not to the bare Calderon--Zygmund kernel.

The Euler energy exponent depends on fixed periodic norms of v,W
and on T, never on p. Nested invariant weights control polynomial
moments and two relative carrier derivatives. The differentiated
local axial factor is p(W-W0)=O(Omega delta R²)+O(delta⁴R⁴);
the exterior has exponential smallness. Duhamel gives the stated
velocity and typed spin errors. Lin transport moves the actual tag.
The finite-order pressure expansion, reference cancellations and
moment solve thus leave O(Omega delta²) clock error and relative
O(delta²) spin error with two carrier jets, below the positive
Omega² delta scale in(9). Nonlinear histories follow on the fixed
interval by taking perturbation amplitude small last.

## 6. Exact joint linear phase interface and its next achievement

0163 gives the same field's prepared actual axial acoustic histories,
physical common velocity, positive c_W² and phase mass rho on
T_ac/|k_ac|. At the linear level Euler is independent of s. Therefore
distinct axial harmonics remain exactly invariant and the COMPLETE
fluid canonical one-form, symplectic form and Jacobi Hamiltonian
pairing between k_ac and p vanish by their axial integral whenever
k_ac != +/-p. This includes pressure and displacement reconstruction;
both preserve axial harmonics. On a common repeated cell this is
ordinary Fourier orthogonality. On Bloch fibers it is the identical
orthogonality of distinct periodic axial harmonics. No pressure
interaction or harmonic mean is dropped to obtain this statement.

Thus the exact linear phase action is the direct sum of the two
actual transported phase actions, on their common time interval.
This establishes coexistence on one field; it does NOT manufacture
a nonzero coupled coarse Fourier transfer between different carriers.
A localized physical tag can observe both phases. Its observation
matrix contains the acoustic tag response as well as c_theta; it is
not diagonal merely because the full action is.

In particular the displacement angular moment
rho integral_tag r cross xi has leading radial row e^-x/2 L_2^1(x),
different from(11). Its symmetric displacement moment, axial centroid
and mean/ambient current also remain their actual integrals. Spin
matching alone does not identify any of these with M theta. Joining
these complete observation/current rows to a shared macroscopic
modulation is the next coupled-action construction. Generic wavevector
isotropy, a whole-Euclidean stationary knotted environment, and a
uniform long-optical-time statement are not consequences of axial
coexistence. The periodic elliptic tubes used here are genuine
invariant material tubes; their geometry is not renamed an EPS knot.

Route verdict: established as stated, with0163's fixed-profile
acoustic theorem and0166's explicit finite-order optical/moment method
as declared unchanged dependencies. Parent objective remains active.
