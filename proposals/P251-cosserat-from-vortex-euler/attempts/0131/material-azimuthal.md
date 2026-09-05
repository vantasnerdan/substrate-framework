# Higher azimuthal Euler waves: physical orientation, spin and action

This attempt tests the proposed positive-curvature higher-m route on
uniform rotation. It preserves exact positive field and observation
constructions while distinguishing them from an autonomous positive
rotor. No accepted C-CST statement changes.

## 1. Full Euler/Bessel field and its actual material displacement

Let `u_*=Omega r e_phi`, `Omega>0`, with constant density rho. Take
`N>0`, transverse radial wavenumber K, integer `m>=3`, and

```
c=N/sqrt(N²+K²),  sigma=±2 Omega c,  omega=m Omega+sigma,
D=4 Omega²-sigma².
```

Here omega is the laboratory field frequency; sigma is intrinsic to
rotating coordinates. With pressure divided by density
`pi=P(r) exp(i m phi+i N z-i omega t)`, the exact velocity profiles are

```
v_r   = i[sigma P'-2 Omega m P/r]/D,
v_phi = [2 Omega P'-sigma m P/r]/D,
v_z   = N P/sigma.
```

The common phase is implicit on these profiles. The complete radial,
azimuthal and axial Euler equations are
`-i sigma v+2 Omega e_z cross v=-grad pi`; incompressibility is precisely

```
P''+P'/r+(K²-m²/r²)P=0.
```

Thus `P=A J_m(Kr)` gives a regular actual Euler field, with pressure,
not an input oscillator. In rotating coordinates the Lin displacement
is `eta=i v/sigma`; in laboratory coordinates it is its physical vector
push-forward. It satisfies the full displacement/velocity relation.
Direct curl differentiation gives

```
curl v=-(2 Omega N/sigma)v.
```

The lower branch is therefore positive-helicity Beltrami with eigenvalue
`sqrt(N²+K²)`. The upper branch has the opposite eigenvalue.

There are two explicit domain choices, with different licenses:

- The whole-space Bessel wave is a generalized mode, not a finite-energy
  localized field. Its K-dependent dispersion does not alone supply a
  compact-cell Bloch band.
- An axial-periodic cylinder of radius R with an impermeable wall is a
  finite-energy exact mode when
  `sigma P'(R)-2 Omega m P(R)/R=0`. The lower condition is
  `c P'(R)+m P(R)/R=0`. Such modes exist: choose a radius just before
  a simple positive zero of J_m, where `-R P'/P>m`, set
  `c=-mP/(R P')` in (0,1), and then `N=Kc/sqrt(1-c²)`.
  This is declared wall geometry, not EPS localization. On a fixed
  cylinder the radial eigenvalue condition also matters; treating K as
  a freely variable Bloch wavenumber would be an additional assumption.

For a nonzero lower-branch endpoint at K=0, scale
`A~K^(2-m)`. Then pressure tends to zero, while the velocity tends to
a multiple of the exact polynomial mode

```
v=(1,i,0)(x+i y)^(m-1) exp[i N z-i(m-2)Omega t],
xi=-i v/(2 Omega),  pi=0.
```

Both full laboratory Euler and Lin equations hold. It has no axis
core-vector tilt for m>=3. It grows transversely and is only a local
comparison without completion. The m=3 Cartesian residual is explicitly
checked in `verify.py`.

## 2. Actual finite-amplitude Euler fields also exist

The real part of each mode remains Beltrami at every time. Hence
`v.grad v=grad(|v|²/2)`, and for every amplitude a,

```
u=u_*+a v,
p/rho=Omega² r²/2+a pi-a² |v|²/2
```

solves nonlinear Euler exactly on its declared domain. The impermeable
mode satisfies the same wall boundary at finite amplitude. This is a
meaningful positive result stronger than a linear tangent. Its field
pattern is an SO(2) relative equilibrium rotating at `omega/m`.

Varying a in this additive family is **not** thereby an isovortical
family on one finite Euler coadjoint orbit. In particular its total
axial angular momentum has no perturbation: the m-nonzero azimuthal
integral of the velocity is zero. This fact is not inconsistent with
the quadratic Noether charge of the isovortical completion in section 5;
that completion includes a second-order mean field. The two constructions
cannot exchange momenta without that completion.

## 3. A genuine m-lobed material orientation

Choose a centered positive material marker at reference time,

```
w_0(r,phi,z)=mu chi_r(r) chi_z(z)[1+d cos(m phi)],
0<d<1,  0<mu<=1/(1+d).
```

The marker can occupy an annulus with an even finite axial window; its
weight denotes the fraction of tagged material, not a change in the
Euler density. Both radial/axial edges may be smoothed. It has a nonzero
complex central shape moment

```
Q_m(t)=integral_tag rho (x+i y)^m,
theta(t)=arg Q_m(t)/m  (locally, modulo 2 pi/m).
```

Transverse centroid variations vanish for m>=3 by angular orthogonality,
so centralization introduces no omitted linear term here. Under a real
spatial rotation by alpha, `Q_m -> exp(i m alpha)Q_m`, and
`theta -> theta+alpha`. This fixes the rotation generator normalization
before any action calculation: theta is an actual material-shape angle.

The reference marker is not stationary in the laboratory:
`Q_m^0(t)=exp(i m Omega t)Q_m^0(0)`. Evaluating the actual Lin
displacement on the transported material yields

```
delta Q_m(t)=exp(i m Omega t)
             rho m integral w_0 zeta^(m-1) eta_+(t),
delta theta=(1/m) Im[delta Q_m/Q_m^0].
```

The common background rotation cancels. Thus the linear *physical*
orientation perturbation has frequency `|sigma|`, not `|omega|`.
The absolute angle is `Omega t+delta theta(t)`; subtracting its known
reference trajectory is not selecting a Floquet winding. Replacing the
material marker by a stationary lab mask would be a different experiment.

For explicit measured coefficients take real P and even chi_z, and put

```
Z_0=integral chi_z,       Z_c=integral chi_z cos(Nz),
R_m=integral chi_r r^(m+1) dr,
R_P=integral chi_r r P dr,
R_g=integral chi_r r^m(P'+mP/r) dr.
```

Choose `Z_c!=0`. The reference moment and variation are exactly

```
Q_0=rho mu pi d R_m Z_0,
delta Q_m/Q_m^0 = (C/Q_0) exp(i sigma t),
C=rho mu m pi Z_c R_g/[sigma(2 Omega+sigma)].
```

The radial formula follows from the actual displacement components
`eta_r=a(r)`, `eta_phi=i b(r)` and
`a+b=(P'+mP/r)/[sigma(2 Omega+sigma)]`.
The zero and 2m angular channels give the displayed coefficient; the
2m channel vanishes for the declared marker. For a sharp annulus
`R_g=[r^m P(r)]_a^b`, an independently useful exact simplification.
Therefore `delta theta=C sin(sigma t)/(m Q_0)` is a nonzero physical
angle whenever R_g is nonzero. The Fourier selection and deformation
identity are checked, rather than inferred from a harmonic's name.

## 4. Complete mechanical spin and its pressure torque

Keep both displacement and velocity terms in the material angular
momentum. Its axial first variation, pointwise before integration, is

```
delta s_z=rho r[2 Omega eta_r+v_phi]
         =rho m P/sigma.
```

Consequently the actual measured tag spin is

```
delta S_z=S_0 cos(sigma t),
S_0=rho mu m pi d Z_c R_P/sigma,
delta S_z_dot=-rho integral_tag partial_phi pi.
```

This last identity is the actual pressure torque on the material tag.
Using only `r v_phi` would fail it. For this prepared phase the measured
spin-rate coefficient is

```
j_tag=S_0 m Q_0/(sigma C)
     =rho mu m pi d² Z_0 R_P R_m (2 Omega+sigma)/(sigma R_g).
```

It is derived, not assigned or renamed from a canonical momentum.
On the lower branch it is negative for small core annuli, where
`R_P R_g>0`. A positive coefficient is also constructible: near an
impermeable lower-branch wall,

```
P'+mP/r=-(1-c)mP/(c r),
```

so a sufficiently thin interior annulus has `R_P R_g<0` and hence
`j_tag>0`. P cannot vanish at that wall together with P', by uniqueness
of the radial ODE. Thus actual material shape angles with nonzero,
even positive, mechanical spin-rate response do exist; rejecting all
higher-m physical observations would discard a useful true result.
Their canonical inertia is still a separate whole-field integral.

At the normalized K=0 polynomial endpoint, `P=O(K²)` and the measured
oscillatory spin is `O(K²)`, while the shape-angle rate has a finite
nonzero limit. At exactly K=0, `2 Omega xi_r+v_phi=0` pointwise:
every tag has zero axial spin variation, independently of its lobes.
There is no nonzero axial spin inertia for this pressure-free endpoint.

Most importantly, wherever the physical material orientation is smooth,

```
sigma²=4 Omega² N²/(N²+K²),
partial_(K²) sigma²=-4 Omega² N²/(N²+K²)²<0.
```

A clever marker changes its measured moment coefficient, not this
intrinsic frequency. Its positive spin variant therefore does not
restore the requested positive spatial curvature.

## 5. Full KKS and the actual SO(2) action

For a finite admissible mode (for example the cylinder above), write
the complex velocity V and displacement Xi=i V/sigma. Use the real and
imaginary columns as the two tangent directions, with amplitude q.
Axial-period/angular orthogonality makes their norms equal and their
inner product zero. Define the actual positive rotating kinetic Hessian

```
h=(rho/2) integral |V|² >0,
beta=2 rho Omega integral (Re Xi cross Im Xi)_z.
```

Taking the imaginary part of the integrated Euler equation multiplied
by V* (the pressure boundary flux vanishes) gives

```
integral (Re V cross Im V)_z
       =-sigma/(4 Omega) integral |V|²,
beta=-h/sigma.
```

This evaluates the sign and normalization of the complete KKS integral;
it does not insert an inertia. With `J=[[0,-1],[1,0]]`, the phase form
is `Omega_KKS=-beta J`. The intrinsic generator is `sigma J` and
`H_rot''=-Omega_KKS(sigma J)=h I`.

The physical rotation action is the vector push-forward
`V(x)->R_alpha V(R_-alpha x)`. For azimuth m its coefficient generator
is **m J**, not J. Its genuine Euler Noether momentum on an equivariant
isovortical chart therefore has quadratic Hessian

```
J_z''=-Omega_KKS(m J)=-m beta I.
```

This is the quadratic restriction of the physical rotation momentum
map, whose linear variation at the axisymmetric base is zero. The
full laboratory quadratic action retains its angular-momentum term:

```
H_lab''=H_rot''+Omega J_z''=-beta omega I.
```

For the lower m>=3 branch omega>0 while beta>0; the laboratory Hessian
is negative. A fixed nonsingular lab scalar observation `theta=c_row.q`
has the exact pulled-back inertia

```
M_lab=-beta/[omega |c_row|²]
     =h/[sigma omega |c_row|²].
```

It is negative on precisely that proposed positive-curvature lower
branch. A transported material angle has a time-dependent lab row;
including that row's connection recovers the intrinsic generator and
positive `h/[sigma² |c_row|²]`, not the negative laboratory one. This
is another reason not to combine the lab frequency with the material
angle's action.

There is an exact sign identity covering both helicities and every m,
away from a zero-frequency/singular chart:

```
M_lab partial_(K²)(omega²)
      =-h/[|c_row|²(N²+K²)] <0.
```

Thus no azimuthal Doppler choice in this uniform-rotation family has
both positive scalar inertia and positive squared-frequency curvature.
When inertia itself has a spatial jet, this coefficient is the
propagation combination `K_2-omega_0² M_2`, not an unjustified
identification with the raw potential coefficient K_2. The sign
identity is independent of the magnitude of the KKS norm or scalar row.
It is asserted for both helicities in `action_verify.py`.

Finally, a finite nonzero wave-pattern orientation is genuinely an
SO(2) angle `alpha=psi/m`, with `q=R(cos psi,sin psi)`. The quadratic
polar action is

```
L=-m I_wave alpha_dot+omega I_wave,
I_wave=beta R²/2.
```

It has a cyclic phase and `alpha_dot=omega/m`; its optical restoring
curvature and its second velocity derivative are zero. It is a relative
equilibrium angle, not a small optical displacement about a statically
oriented core. The charge `-m I_wave` is the fixed-orbit quadratic
momentum described above, not the unchanged global momentum of section
2's additive finite-amplitude family and not the linear tagged spin.

## 6. Decisions and the next materially different construction

Candidate A's positive-curvature route has verdict **refuted with
mechanism**: material transport removes the Doppler shift. Its strongest
meaningful positive evidence is an actual material shape angle with full
pressure torque and a nonzero spin row, including a positive-spin annular
choice. Candidate B's positive optical-rotor route has verdict **refuted
with mechanism**: the true pattern phase is cyclic, and the quadratic
lab scalar version has the wrong inertia sign whenever curvature is
positive. Its strongest positive evidence is an exact nonlinear Euler
relative equilibrium and a correctly normalized physical field angle.
The general sign identity covers both helicities and all m in this
uniform-rotation azimuthal route. `evidence_scope: REPRESENTATION_SCOPED`
refers to these declared uniform-rotation domains, not the parent class.

Candidate C, the tilted-vorticity m=1 repair, retains the separately
proved 0125/0128 physical observation and localization interface. This
attempt does not duplicate its review or close its remaining periodic
Euler problem. The failure-derived next candidate must alter the actual
intrinsic spatial operator, not merely its laboratory Doppler label:
for example coherent nonuniform/elliptic fixed-curl cells with the full
physical observation, or genuine intercell mode mixing carrying its
stress/shape branches. A concrete new background already exists here:
linearize about the finite-amplitude relative equilibrium of section 2,
which is nonuniform and stationary in its pattern frame. That is an
actual new operator, unlike relabeling a uniform-rotation frequency;
its SO(2) neutral mode and pressure/shape currents remain part of the
construction. The full Euler/Bloch construction of 0125 is the
appropriate interface. No proof here rules out those routes, and the
parent remains active.

`verify.py` records 17 exact field/material checks in `first-run.txt`.
`action_verify.py` checks KKS, physical rotation normalization, both
helicity sign identities, phase degeneracy, a moving-tag mutation and
an actual nonlinear Cartesian mode in `action-first-run.txt` (11 checks).
`moment_verify.py` derives the measured-inertia and KKS normalizations,
including the positive wall-marker sign and exposing missing-vorticity
mutation, in `moment-first-run.txt` (6 checks). No numerics, canonical
edits or accepted claim changes were required.
