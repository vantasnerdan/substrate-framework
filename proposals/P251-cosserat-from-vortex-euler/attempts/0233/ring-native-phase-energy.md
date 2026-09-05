# Actual ring-native phase/energy return and the first carrier forms

The first physical-history calculation does not automatically normalize
its inherited action. This body derives the additional actual Kelvin
sectors and the normalization matrix on the same literal-curl core.
It imports no C016 phase-control theorem.

## 1. Actual first carrier forms of the original lift

For Xi_k=(Jgrad S-ik lambda grad f,lambda S), S=Delta_1 f,
let Lambda=pi rho L with L the actual axial measure (2pi R in the
ring limit). Compactness removes the radial boundary terms. With m=-1
the complete KKS and Jacobi energy MATRIX coefficients are

    beta0=Lambda lambda^2 integral r Omega A^2,
    beta1=-2Lambda lambda^2 integral r W A^2,
    H0=Lambda lambda^2 integral r Omega^2 A^2,
    H1=-4Lambda lambda^2 integral r Omega W A^2.           (1)

The phase here is the full Kelvin cotangent, equivalently the KKS form
on the actual fixed-circulation leaf. The energy is

    H[w]=rho/2 integral[|w|^2-lambda^-1 w.curl w].         (2)

Although lambda is constant only on the generator support, (2) is
valid: the integration involving Xi cross u is supported there, and
P(Xi cross omega)=w before integrating against curl w. Pressure tails
are retained in w. The scalar first pressure correction contributes
only the exact radial boundary derivative in beta1 and H1. It is not
set to zero pointwise.

In the uniform-core limit (H/beta)'=-2Wc, whereas removal of the
unwanted rows of the original physical tag gives its Doppler slope
-Wc. This is a concrete first-order action/observation discrepancy,
not permission to rename the observed clock. Nor is it resolved by
assuming a signed average cancels quadratic energy. The general axial
lift and the actual return below address these distinct equations.

## 2. Two actual Kelvin signatures

In the same column take, for any real c,

    Xi_c=(Jgrad S,c lambda S),
    b=-lambda Omega S_theta,
    d=lambda^2 P_perp(r Omega S e_r),
    w_c=b e_zeta+(1-c)d.                                 (3)

This is the actual pressure projection of Xi_c cross omega. The full
real-pair phase and energy matrices are exactly

    beta_c=(2c-1) beta_1,
    H_c=(2c-1) H_1+(1-c)^2 D,
    D=rho times the real-column integral |d|^2 >=0.       (4)

Here beta_1,H_1 denote the c=1 zero-carrier forms, not derivatives
in (1). The exact identity behind (4) is

    lambda^-1 integral d.curl(b e_zeta)=integral |b|^2.

It follows before any radial approximation by projecting against a
solenoidal curl and using the radial component of curl(b e_zeta).
Both horizontal and axial velocities, and their helicity cross term,
are retained.

If the support is in r<=b0, the projection contraction gives

    D/H_1 <=lambda^2 b0^2/m^2.                           (5)

Thus c=1 has positive phase and energy; c=0 has negative phase and
strictly negative energy for lambda b0<1. These are genuine Euler/Kelvin
signatures, not negative probabilities or an assigned rotor mass.

The more general h lift in the companion body has

    beta_h=2Lambda lambda integral r Omega A h
                         -Lambda lambda^2 integral r Omega A^2,
    H_h=2Lambda lambda integral r Omega^2 A h
                         -Lambda lambda^2 integral r Omega^2 A^2
                         +rho integral_real |d_h|^2,
    d_h=lambda P_perp[r Omega(lambda S-h)e_r].            (6)

At uniform Omega0, H_h=Omega0 beta_h+rho integral_real |d_h|^2.
The extra positive horizontal energy is real; canceling the physical
secular row alone does not cancel it. An off-tag H correction can make
beta_h positive without changing the physical rows: use an off-tag
component of A and choose a compact H variation whose pairing with
Delta_1^2 A has the desired nonzero sign. The actual total phase and
energy are then measured by (6), not guessed from the observed clock.

## 3. Compact annular controls and two different action ratios

Put each control in an annulus disjoint from the observed tag and the
other control supports. The c=0 pressure in (3) has both an interior
and an exterior harmonic tail. Their two coefficients vanish when

    integral Omega A dr=0,
    integral r^2 Omega A dr=0.                           (7)

Indeed div(rOmega A e_r)=2Omega A+rOmega'A+rOmega A'.
Its regular m=1 interior and exterior Green moments reduce by parts
to the two rows (7), with nonzero fixed factors. Then d is compact in
that annulus. At the reference column all phase/energy cross forms
between disjoint controls vanish exactly. Additional finite moments
for time/parameter accuracy remain compatible with a nonzero profile;
the strict sign bound (5) holds independently of those profile choices.

Write nu_+=H_+/beta_+ and nu_-=H_-/beta_->0. The second ratio obeys
nu_-=nu_base-D/beta_+, while nu_+ is the positive weighted mean of
Omega on its annulus. Since Omega decreases quadratically near the
core, two well-separated fixed scaled annuli give either ordering of
nu_+ and nu_-. For nu_+<nu_-, put the negative branch sufficiently
closer to the axis that the bound D/beta_+ from (5) is smaller than
the computed decrease to the positive outer annulus. A fixed radius
ratio, followed by lambda b0 sufficiently small, suffices. The reverse
ordering follows with the positive branch inside the negative branch.
The gap is of order (lambda a)^2 with its actual positive lower bound;
it is not an assumed spectral separation of the full Euler operator.

Whiten the absolute energy matrices of one positive and one negative
pair. They become +I,-I, with phase coefficients 1/nu_+ and -1/nu_-.
Equal squared amplitudes give zero energy and phase

    a[1/nu_+-1/nu_-]J.                                  (8)

The two available orderings give both signs in (8). Conversely squared
amplitudes a and a nu_-/nu_+ give zero phase and energy

    a[1-nu_-/nu_+] I.                                   (9)

Again both signs are available. The actual two-by-two normalization
matrix is diagonal in these zero-energy and zero-phase coordinates,
with the two explicitly nonzero gap coefficients (8)-(9). This gives
independent control of phase and scalar inherited energy, not a generic
rank assertion. To cross zero smoothly in a parameter, use a fixed
positive baseline in both opposite-sign copies and vary their squared
amplitudes by the corresponding signed amount. Square roots then stay
smooth on a specified small parameter window. All ensemble weights
can remain positive.

## 4. Transfer and the physical-output boundary

Complete every annular generator by its compact potential and the exact
Piola map, as in 0222. Compute its actual global velocity P(Xi cross
omega_R), full phase and (2). The strict sign margins and the finite
normalization minor persist when R is selected after their actual
conditioning constants. Small cross forms at finite R are included in
the same finite-dimensional implicit equations; they are not discarded
by pretending that velocities of disjoint generators are disjoint.
The ordinary IFT with the nonzero minor of (8)-(9) retunes the full
actual phase/energy. Its inverse-gap amplitude cost is retained.

The controls have zero initial tag displacement because their generators
are off-tag. Their later physical output is controlled, not presumed
zero. The compact pressure recursion of 0226 applies on annuli with
both interior and exterior Green moments at every retained order;
Laurent/logarithmic harmonic weights replace the single exterior
moment where necessary. Polynomial coefficient operations and a finite
number of inverse Laplacians consume only finitely many such moments.
This constructs compact off-tag coefficient histories to the selected
finite order, and the actual full Euler/Lin residual gives the remainder.

Choose that finite order above the actual amplitude/conditioning powers
before taking the small inner-core and curvature parameters. This
provides a controlled physical-output return, not an exact invisible
Euler invariant manifold. The first physical-history matrix and its
normalization are the present scope. A simultaneous infinite-order
limit, a common-K stationary ensemble, and the acoustic/current join
still require their own completed estimates.
