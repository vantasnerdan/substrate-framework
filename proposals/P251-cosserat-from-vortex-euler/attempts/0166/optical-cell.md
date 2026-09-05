# Actual material quadrupole and spin on the triangular constant-curl cell

## Statement and observation

Fix Psi,lambda,L_s,rho>0 and a fixed interval0<=Omega t<=T. On the exact
periodic field of `core-geometry.md`, take k=-p=-2pi N/L_s with N a
sufficiently large positive integer. Construct two smooth fixed-Kelvin
initial displacements, their actual Euler/Lin histories, and a smooth
nonnegative material label in one elliptic core. Its observable is the
ordinary centered Euclidean quadrupole angle

    Q(t)=rho integral_labels [(x-X_D)+i(y-Y_D)]² w0 da,
    theta=(1/2) delta arg Q.                              (1)

Its response to a physical rigid rotation is exactly one. No metric
calibration or chosen Floquet logarithm is used. It is a registered
material shape, not an asserted absolute vortex director. Its phase
action has positive scalar mass and the observed phase-clock curvature

    gamma=2Omega+(Omega/3)delta+O(Omega delta²),
    p² partial_p² gamma²=Omega² delta+O(Omega² delta²)>0,
    delta=sqrt(2lambda/p).                               (2)

A fixed radial marking makes literal axial material spin eta times its
canonical momentum to relative O(delta²), uniformly on this interval
and through two relative carrier jets. eta=1 or the predeclared eta=1/2
row normalization may be selected. No pair construction is inferred.
This does not assert a positive
acoustic branch of the same field; that remains the separate0163 route.
Relative error means the operator norm of the two-phase observation row
divided by its nonzero characteristic action/angle/spin scale, not division
by a scalar oscillation at its zeros. The two relative carrier derivatives
are p partial_p and its square, with the physical tag fixed at p_star.

## 1. Complete pressure and Lin reconstruction

With axial factor exp(iks), remove only the translation phase kW0.
The normal velocity V, axial velocity b and full pressure P satisfy

    D_k V+(V.grad)v+grad P=0,
    D_k b+V.grad W+ikP=0, div V+ikb=0,
    D_k=partial_t+v.grad+ik(W-W0).                        (3)

All pressure inverses are those of the full triangular normal torus.
Retain b=i div V/k. The physical displacement obeys exact Lin
reconstruction D_t xi-(xi.grad)u=v_pert. Its initial velocity is
prepared by v_pert=P_Leray(xi cross omega0), and Euler transports that
specified Kelvin data. A restricted orbit ODE is not substituted.

At v=-Omega Jr, W=W0+lambda Omega R²/2, divergence of the two normal
equations and the axial equation in(3) give

    (k²-Delta)P=2Omega curl V+2ik grad W dot V.             (4)

Both the strain and axial-shear contributions are present. With k=-p,
ell4=2/(lambda p³), r=ell R and time normalized by Omega, the actual
pressure resolvent expansion gives

    L=L0+delta L1+delta² L2+O(delta³),
    L0=partial_theta+J,
    L1 V=iR²V/2-grad curl V,
    L2 V=i grad(R dot V)-(1/2)grad Delta curl V.           (5)

Its remainder is the full inverse applied to the next differentiated
Gaussian source, estimated below. The actual cubic radial velocity
enters at delta³, and the first triangular C6 term at delta6.

Take e_+=(1,i)/sqrt(2), m=2,n=2,l=m-1 and

    F=R^l exp(-R²/2)L_2^1(R²) exp(il theta).

The resonant equation is [-Delta_R+R²]F=2(2n+m)F, yielding relative
Eulerian frequency omega=(2-m)Omega-(2n+m)Omega delta+O(Omega delta²).
The first nonresonant velocity is delta e_- D_+²F/8, at angular l+2.
The Kelvin inverse's pressure projection is included in xi and KKS.
At second order solve the nonresonant L0 equations with their nonzero
integer denominators. The opposite resonance at l-2 is absent by
circular angular grading. The general pressure cancellation derived
in0155 reduces here to DeltaG=0. Retain the resonant plus-l radial
source of order delta²; it is not an isolated-spectrum assertion.
Every other residual is O(delta³).

## 2. Literal spin and its exposed cancellation

For any smooth material tag the centered spin variation is

    S=rho integral_labels [xi cross u0
                +r cross(v_pert+xi.grad u0)]_s w0 da,

with actual centroid terms included where they contribute. This formula
does not confuse Eulerian and material time derivatives of xi. For the
radial comparison field it is exactly

    S=rho integral_labels R_phys
                     [v_theta+(2O+R_phys O')xi_r] w0 da.  (6)

At O=-Omega and xi=i e_+F/(2Omega), the integrand vanishes POINTWISE.
The first spin is the order-delta pressure torque. The order-delta²
plus-l residual has the same zero principal spin, so its observed spin
starts at delta³. Hence it gives relative spin error O(delta²), not
O(delta). This typed cancellation is essential: an unspecified small
velocity residual would not imply the spin bound. The actual O' term
and triangular corrections are retained at their derived orders.

## 3. Full-cell KKS and the physical clock

Give the leading velocity a declared amplitude V_*. The phase columns
are Re xi, Im xi. Integrate KKS over the COMPLETE action cell:

    beta=rho integral_cell omega0 dot (Re xi cross Im xi),
    beta=-3pi rho V_*² L_s ell²/(4Omega) [1+O(delta)]<0.   (7)

Here integral_R² |F|²=3pi. The exact beta, not this approximation,
enters the moment solve. All outside-core Leray velocity and periodic
pressure images are retained. The tag is not the action domain.

Choose radial label density chi(x), x=R², an axial label cutoff of width
c/p_star, and modulation b(x) cos(2theta-p_star s), with p_star fixed at
the selected integer. The complete density is nonnegative. Its small
reference Q rotates at -2Omega at principal order. Its mode numerator
also includes actual axial material travel. For reference weight chi=x,
the exact Laguerre Laplace ratio is <x>/2=6+1/3, yielding(2) rather
than the unobserved eigenfrequency.

To keep that coefficient AND its carrier curvature with a compact
smooth tag, start with chi0=x times a sufficiently distant scaled cutoff.
Define three linear functionals

    R=exp(-x/2)L_2^1(x)(x/2-19/3), D0=(3/2)x partial_x,
    F_j(h)=integral h x D0^j R dx, j=0,1,2,
    M_jb=F_j(q_b), chi=chi0-sum_b[M^(-1)F(chi0)]_b q_b.  (8)

The three analytic densities have distinct polynomial degrees after
their Gaussian is removed, so three interior bumps give invertible M.
The tail is exponentially small, so one finite cutoff preserves chi>=0
and a nonzero angle denominator. A common rescaling may bound chi by1.
This fixes the phase coefficient and its two carrier jets, not only its
value. `three-jet-cutoff.md` proves the complete construction. Without(8),
retain the actual nearby positive cutoff coefficient and its derivatives.

The MODE is instead cut off at a fixed small PHYSICAL core radius.
Its scaled radius tends to infinity and the residual is polynomial
times exp(-c/ell²), including carrier jets. The tag cutoff in(8) does
not pin or localize the Euler mode.

## 4. Fixed-tag time/carrier solve and action

The actual pressure-spin polynomial is

    P=L_2^1-2(L_2^1)'=x²/2-5x+9.                         (9)

Fit its value and first slow-time coefficient, with their first two
relative carrier derivatives. The six rows are P,DP,D²P,xP,D(xP),
D²(xP), where D=c0+(3/2)x partial_x acts also on the Gaussian.
Changing the scalar prefactor c0 only performs triangular row operations.
Add reference rows1,x. The verifier computes a NONZERO eight-row minor
for n2,m2; no minor from the different0155 mode is borrowed.

The axial tag occupies O(1/(p_star L_s)) of the full cell. Its required
dimensionless reference moment is accordingly O(delta/(p_star L_s)),
which is O(delta³) at fixed geometry. The reference x row cancels the
first cubic-core phase. Remaining radial/C6 reference error is O(delta6),
hence O(delta³) after dividing by the small reference moment. Reference
normal centroid is exactly zero by inversion; the complete centered
linear observable is nevertheless retained. No other current row is
set to zero by convention.

After dividing by the explicit leading spin and axial-volume factors,
the exact map is an O(delta) perturbation of the invertible eight-row
map. It includes second-order pressure, exact Kelvin initial data,
beta and the fixed axial filter. The finite-dimensional implicit-
function theorem fixes b for S=eta Pi_theta. A common reduction of
marking and reference amplitude preserves nonnegative density. The
first-time match leaves relative O(delta²) error on fixed Omega t;
the two carrier jets are fitted rather than inferred from t=0.

For the actual measured angle row c_theta and its phase rate gamma,
the same moving phase action is

    M=-beta/(gamma |c_theta|²)>0,
    L=M/2[(theta_dot-(partial_t log|c_theta|)theta)²
                                                   -gamma² theta²].

Time and parameter connections stay in the action. Equations(2),(7)
give its positive sign and natural-scale curvature; they do not turn
this moving observation into an autonomous continuum equation.

## 5. Actual full-pressure finite-time estimate

Let E=3Psi-psi, a genuine invariant of the normal flow, and let
R_E=sqrt(2E/Omega), with physical length units. Use the weight
exp(sigma R_E/ell), smoothly flattened as a function of E inside
the elliptic cell. Its transport derivative is zero and its logarithmic
gradient is bounded by C/ell. The periodic pressure Green function has
exponentially decaying modified-Helmholtz images. Weight conjugation
replaces exp(-p distance) by at most exp[-(p-C/ell)distance]. Since
p ell tends to infinity, weighted pressure and commutator bounds are
uniform and retain EVERY image.
For the order-zero pressure return grad(k²-Delta)^(-1)div, apply its
unweighted L2 bound first. In the conjugated commutator the factor
w(x)/w(y)-1 cancels one order of the two-dimensional1/r² kernel
singularity. Its near-diagonal integral is O(1/(p ell)); the remaining
tail is bounded by the displayed exponential. Thus Schur's estimate is
used on the integrable commutator, not incorrectly on a bare singular
Calderon--Zygmund kernel. The local delta part has zero commutator.

The weighted Euler energy exponent is C_{Psi,lambda,T}, independent of p.
Finitely many nested weights control polynomial R moments and two
relative carrier derivatives. The latter generate the centered local
factor p(W-W0)=O(Omega delta R²) and controlled pressure derivatives;
outside the weighted core there is an exponential tail. All global
background norms are fixed bounded periodic norms, not unknown constants
attached to an improving local approximation.

Duhamel for the exact Euler group controls the typed residual in(5),
and exact Lin transport moves the actual labels. The cancellation in(6)
and the fitted time/carrier rows give relative O(delta²) spin/action
error and O(Omega delta²) clock error with two relative carrier jets.
This is smaller than the positive Omega² delta curvature. It is a
finite-time statement, not an invariant nonlinear mode subspace.
Actual nonlinear histories follow on the fixed interval by taking
the perturbation amplitude sufficiently small last.

## 6. Same-field and current boundaries

Integer carriers are smooth finite-action fields on the SAME triangular
three-torus as0163. Continuous neighboring carriers are Bloch fibers
with cell action density. A whole-space packet needs its own envelope
normalization and material observations. This proof does not identify
those different objects or a distant knotted tube with this core.

The displacement angular moment L_xi=rho integral_tag r cross xi is
not assumed to be M theta. Its leading radial functional is
exp(-x/2)L_2^1(x), different from the matched pressure polynomial(9).
Its symmetric displacement moment and axial centroid row remain actual
integrals for the joint current consumer. Spin matching alone does not
manufacture their equality. Acoustic sign, generic-wavevector response,
absolute frame and complete current closure remain separately active.
