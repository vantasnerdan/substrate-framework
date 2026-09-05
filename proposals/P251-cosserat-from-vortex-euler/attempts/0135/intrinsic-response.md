# Intrinsic response: full reaction and a material surface-mode construction

This attempt changes the intrinsic Euler operator, not a trial gradient
modulus. Its two complementary results are an exact dynamical reduction
interface for the smooth polygon and an actual Rankine surface mode with
positive intrinsic spatial curvature. The latter admits a positive
mechanical/canonical moment match at a fixed nonzero axial carrier. The
smooth force-free continuation is a next construction, not an imported
conclusion. Accepted v0.175.0 claims remain unchanged.

## 1. Retain the optical-frequency reaction

For real mode coordinates use
`L=-z^T Omega zdot/2-z^T H z/2`, with `Omega^T=-Omega`.
The equation is `Omega zdot+H z=0`; hence the `exp(-i nu t)` pencil is
`D(nu,k)=H(k)-i nu Omega(k)`. Partition it into the retained centroid
sector c and the full remaining Euler sector f. Whenever `D_ff` has an
inverse on the declared quotient, exact elimination gives

```
E = D_cc-D_cf D_ff^-1 D_fc.
```

Every mixed symplectic row occurs inside D; replacing the inverse by a
static H inverse is a different calculation. In an even-k channel,
write `D=D0+k² D2+o(k²)` at fixed nu. With `R=D_ff,0^-1`,

```
E2 = D_cc,2-D_cf,2 R D_fc,0-D_cf,0 R D_fc,2
     +D_cf,0 R D_ff,2 R D_fc,0.
nu2 = -w* E2 v/(w* (partial_nu E0) v).
```

The last expression assumes a simple optical pole and the left/right
null vectors at that pole. Without an even channel, retain the first
jet and its eigenvector correction as well. The same resolvent gives
the full physical observation `O_c-O_f R D_fc`; the observation is not
necessarily the uncorrected centroid coordinate.

Attempt0041 supplies the actual self/mutual Biot--Savart second
variation on its affine twist/radius family. In the two alternating
triangle coordinates its positive logarithmic terms are
`Cq,Cx = 3 rho Gamma² d² log(d/epsilon)/(8 pi)+O(1)`.
With *that family's* canonical pair the frequency coefficient is
`(h Cq+K Cx)/P0²`. This is a valid contribution to D_cc,2.
Attempt0036 also supplies a genuine planar radial-core gap after
translations and its stated two-dimensional Casimir directions are
separated. It does not supply an inverse for all three-dimensional
axial/stretch directions of D_ff at the optical frequency. In
particular, k-dependent divergence completion is singular on some
coordinates that become radial Casimir changes at k=0.

Thus the polygon route's missing construction is precisely the actual
three-dimensional optical-frequency inverse and its mixed second jet,
or a bound proving its contribution smaller than the positive
logarithmic term. The named positive logarithm is not discarded, nor
is the missing inverse replaced by a fixed-core assumption. A two-mode
exact example in the verifier shows why a positive retained gradient
alone cannot decide the result: eliminating a dynamical oscillator
gives `nu2=(a-c²/(L-nu0²))/(2 nu0)`, which can have either sign.

## 2. Actual Rankine surface branch, with the material clock

Use the existing `rankine_modes` Euler equations: radius a, density
rho, core `u=Omega r e_theta`, exterior `u=Omega a² e_theta/r`,
`Omega>0`. Vorticity is a finite-radius patch; velocity is continuous.
This is exact Euler, not a smooth EPS field. Perturbations have phase
`exp(i m theta+i k z-i omega t)`, m>=2. Put
`x=|k|a`, `s=(omega-m Omega)/Omega`, and `sigma=Omega s`.
Exterior perturbations are irrotational; match actual pressure and
normal displacement, including the exterior Doppler factor. The
canonical boundary determinant is

```
F=(4-s²) K+s² J-2 m s,
J=lambda J_m'(lambda)/J_m(lambda),
K=x K_m'(x)/K_m(x), lambda=x sqrt(4-s²)/|s|.
```

Near the simple surface root s=-1,

```
J=m-3x²/[2(m+1)]+o(x²),
K=-m-x²/[2(m-1)]+o(x²),
F_s(-1,0)=-6m,
s=-1-x²/[2(m²-1)]+o(x²).
```

For m=2 the next remainder can contain x^4 log x. In particular

```
omega²=Omega² (m-1)²-Omega² x²/(m+1)+o(x²),
sigma²=Omega²+Omega² x²/(m²-1)+o(x²).
```

The first sign is not a physical material-angle verdict. A transported
m-lobed marker rotates with the core at Omega. Its central complex
shape moment has background factor `exp(i m Omega t)`; the ratio of
its variation to its background removes precisely that factor. Its
angle perturbation has the intrinsic frequency |sigma|. This is an
actual transported observation, not a laboratory frequency relabeling.
For m=2 the physical squared-frequency coefficient is therefore
`Omega² a²/3>0`.

## 3. Complete KKS/action, including the exterior pressure response

Let P be the real radial pressure divided by density in the core,
proportional to `J_m(lambda r/a)`, normalized to approach `p r^m`
with p>0 as x->0. The full core Euler equations give

```
v_r=i(sigma P'-2 Omega m P/r)/(4 Omega²-sigma²),
v_theta=(2 Omega P'-sigma m P/r)/(4 Omega²-sigma²),
v_z=k P/sigma.
```

The actual rotating-material displacement is `eta=i v/sigma`.
Write `eta_r=A(r)`, `eta_theta=i B(r)` before the common phase. Then

```
A=(2 Omega m P/r-sigma P')/[sigma(4 Omega²-sigma²)],
B=(2 Omega P'-sigma m P/r)/[sigma(4 Omega²-sigma²)],
beta=4 pi rho Omega Lz integral_0^a A(r) B(r) r dr.
```

Lz is the length used to normalize the axial mode; at k=0 one may use
action per unit length. This is the KKS integral
`rho integral omega0 dot(Re eta cross Im eta)`. Exterior base
vorticity is zero, so it contributes no missing volume term to this
integral. The normal-displacement and pressure matching already include
its full velocity/pressure response. A divergence-free generator can
be extended through a finite exterior collar; its exterior relabeling
does not change the KKS integral. The actual exterior material
displacement need not be periodic in labels because of differential
rotation. No periodic exterior-label claim is used here.

At x=0, `A=B=m p r^(m-1)/[sigma(2 Omega+sigma)]`, so beta is
strictly positive. Continuity keeps beta positive for small nonzero x.
The two real mode columns have `Omega_KKS=-beta J2`. Their actual
laboratory evolution has generator `omega J2`; the Hamiltonian Hessian
is consequently `-beta omega I2`. The physical SO(2) generator is
`m J2`, so the rotating-material Hamiltonian is
`H_rot''=-beta sigma I2>0`. This obtains the energy from the actual
Euler mode and symplectic form. It does not import the impermeable-wall
kinetic-energy formula, which would omit the surface pressure terms.

If the measured material angle has row length c, its exact scalar
action after eliminating the conjugate quadrature is

```
M_can=-beta/(sigma c²)>0,
L_theta=M_can [theta_dot²-sigma² theta²]/2.
```

Its mass also has a spatial jet. The normal-form propagation
combination is `K2-sigma0² M2=M0 partial_(k²)(sigma²)>0`, not an
assumption that M2 vanishes. Thus the actual surface-mode action does
have positive intrinsic optical curvature. Mechanical spin equality
is a separate measurable condition constructed next.

## 4. A nonnegative transported tag with two signed radial moment controls

Inside the uniform core choose the actual material mass fraction

```
w0(r,theta,z)=mu chi(r) chi_z(z) [1+eps b(r) cos(m theta)],
0<=chi,chi_z<=1, |eps b|<1, 0<mu<=1/(1+|eps b|_infinity).
```

The signed function b changes lobe orientation with radius; it does
not make the material density negative. Take an even axial window with
`Z0=integral chi_z>0`, `Zc=integral chi_z cos(kz)!=0`, and m=2
if an ordinary covariance orientation is desired. Angular symmetry
centers this tag and its first variation. Define

```
Q = integral chi b r^(m+1) dr,
T = integral chi b r P dr,
G = integral chi r^m(P'+m P/r) dr.
```

The complete material calculation of 0131 applies inside the core:

```
Q_m^0 = rho mu pi eps Z0 Q exp(i m Omega t),
delta theta = c sin(sigma t),
c=Zc G/[eps Z0 Q sigma(2 Omega+sigma)],
delta S_z = rho mu m pi eps Zc T cos(sigma t)/sigma.
```

Spin includes both position and velocity variations. Pointwise its
coefficient is `rho r(2 Omega eta_r+v_theta)=rho m P/sigma`;
its time derivative is the actual pressure torque on the tag. Omitting
the deformation term reverses relevant conclusions. Therefore

```
j_tag = rho mu m pi eps² Z0 Q T (2 Omega+sigma)/(sigma G).
```

At k=0, P=p r^m implies T=p Q for *every* signed b. Hence
`j_tag<0` whenever Q is nonzero. Radial lobe reversals cannot repair
this endpoint; the ordinary positive-radial marker has the same sign.

At every sufficiently small fixed k!=0, the rows Q and T are linearly
independent on smooth radial bumps. Indeed `P/r^m` is not constant on
any radial interval, as follows directly from its Bessel ODE with
nonzero lambda. Choose two narrow smooth bumps in disjoint radial
intervals on which their weighted mean values of P/r^m differ. Their
2x2 moment matrix is invertible. It can therefore prescribe a nonzero
Q and the exact ratio

```
T/Q = -beta Z0 sigma²(2 Omega+sigma)
      /[rho mu m pi Zc² G].
```

Substitution gives **j_tag=M_can>0**, with the same exact Euler mode,
material angle, full pressure torque and KKS/action. Scale both bump
coefficients together to enforce |eps b|<1 without changing the ratio.
All coefficients are actual profile integrals. This is observation
design, not a fitted frequency or a postulated inertia. The tag mass
fraction scales the measured moment and is explicitly included.

The conditioning matters. With the stated pressure normalization,

```
P/r^m=p[1-lambda² r²/(2(m+1)a²)+O(lambda^4)].
```

The two-row determinant is O(k²). A bounded family b(k) that maintains
the opposite-sign ratio must have Q=O(k²); its angle chart degenerates
at k=0. Thus this construction establishes a fixed nonzero-carrier
match, not a nonsingular matched local rotor all the way to k=0.
A finite carrier can be chosen before a later macroscopic sideband
expansion. Such an expansion must retain both counterpropagating
carrier branches and their group response; identifying a standing
carrier with one autonomous two-field band is not done here.

## 5. Result and activated continuation

The full-frequency elimination identity and the finite-core logarithmic
contribution are established at their stated scope. The polygon's
actual three-dimensional mode transfer remains blocked on its named
frequency-dependent inverse, not refuted.

The Rankine route establishes an actual finite-radius Euler surface
mode with positive material-angle action and positive intrinsic
spatial curvature, and constructs a nonnegative material tag whose
physical spin matches the canonical momentum at any sufficiently small
fixed nonzero carrier. Its k=0 matched-tag continuation is refuted for
the declared marker family by exact proportionality of the radial rows.
This is stronger than a positive trial Hessian and weaker than the
requested smooth stationary EPS/coarse coupled mode.

The next source-compatible route is already under construction by the
parent in 0136: smooth the vorticity gradient in a generalized force-free
column, retain the radial mode and its pressure/moment rows, and transfer
the isolated pole and second jet at a fixed nonzero carrier. Separate
work must then assemble coherent physical translation and axial
sidebands. The exact smooth polygon remains a distinct interaction
candidate. No source here supplies a full three-dimensional EPS Bloch
pole theorem, and no accepted claim is changed.
