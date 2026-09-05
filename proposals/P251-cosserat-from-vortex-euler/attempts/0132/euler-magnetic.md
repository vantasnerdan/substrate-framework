# Exact Euler response, and what transfers from the magnetic analogy

## Source and typed comparison

Primary source: H.K. Moffatt, JFM166 (1986),359–378,
[doi:10.1017/S0022112086000198](https://doi.org/10.1017/S0022112086000198),
author PDF at the URL in README. The source distinguishes magnetic
field transport from Euler vorticity transport (§3). Its elastic-wave
argument is magnetic (§6). The proposed Euler mean-field evolution in
§7 explicitly leaves microscopic feedback unresolved; its general
instability extension is a conjecture. Neither is used here as a theorem
about the active Euler continuum.

PDF retrieved from the public HTTP author archive after HTTPS certificate
and web-fetch failures; SHA256
`c518d49675c128d7327666535cc3eed5d3affa390df512456b89c8c9d704ffc7`.
The PDF is a local reading artifact, not redistributed in the repository.

## Independent exact Euler calculation

Choose reciprocal-length unit alpha=1, density rho=1 and the actual
steady ABC field

    u=(U3 cos z+U2 sin y, U1 cos x+U3 sin z, U2 cos y+U1 sin x),
    curl u=-u, p=-|u|²/2.

For nonresonant rational 0<|k|<1/3, let eta=(cos kz,sin kz,0).
The common periodic-cell average includes every product mode. The
constant-density Kelvin velocity variation is v=-P(eta cross u), with
the full mean-preserving pressure projector P. Put
b=curl(eta cross u), so curl v=-b. The coadjoint energy Hessian and the
magnetic energy Hessian at the SAME background are respectively

    H_E=<|v|²-v.b>, H_M=<|b|²-v.b>.

Thus H_E+H_M=<|v-b|²>; this identity supplies no individual positive
sign. Exact finite Fourier convolution in the existing accepted API gives

    H_E=-(U1²+U2²) k³/[2(1+k²)],
    H_M=(U1²+U2²) k²/2.                              (1)

Restoring units multiplies both by rho and the displacement amplitude
squared; H_E carries alpha k³/(alpha²+k²), H_M carries k², with U_i
the dimensional velocity amplitudes. This is a comparison of functionals,
not an insertion of a magnetic field or magnetic force into Euler.

The complete linear Euler acceleration at this initial Kelvin velocity
is L v=-P div(u tensor v+v tensor u). Extracting its slow Fourier modes,
write A=U1²+U2² and D=U1²-U2². Then

    (L v)_slow=k²/[2(1+k²)]
       ((A k+D) cos kz,(A k-D) sin kz,0).             (2)

Consequently <eta.(L v)_slow>=-H_E exactly. In the isotropic-amplitude
case U1²=U2² this is U1² k³ eta/(1+k²). It is not the magnetic
restoring acceleration, and it is not by itself an equation for eta(t).

The anisotropic expression is independently obtainable without Fourier
machinery. With only U1=1 nonzero, set Afield=eta cross u and

    phi=-(1-k) sin(kz) cos x/(1+k²), v=-Afield+grad phi.

Direct differentiation gives div v=0. Direct x integration gives

    <u_z v_x>_x=-k(k+1) sin(kz)/[2(1+k²)],
    <u_y v_z+u_z v_y>_x=k(k-1) cos(kz)/[2(1+k²)],

which reproduces (2). The first run's anticipated anisotropic force
omitted these terms; the canonical calculation and this independent
Cartesian derivation repair that attempted expectation. The isotropic
specialization and (1) are unaffected. No claim about unrelated results
in the primary paper follows from that repair.

## Positive scope and next construction

Counterhelical k<0 gives a strictly positive coadjoint Hessian on the
two-dimensional real macro plane when A>0. Nevertheless its KKS matrix
is zero: the two quadrature displacements cross to a constant axial
vector, while the background vorticity has zero mean. This positive
quadratic energy alone is not a normalized physical oscillator. A mean
velocity datum and the induced microscopic response must be retained.

Likewise the signed cubic term is not an isotropic positive quadratic
shear modulus. Cancellation of this term in a mirrored ensemble cannot
manufacture the missing even-order modulus. That statement concerns this
specific bare-displacement calculation, not all possible corrected cells.

The exact Euler energy, initial stress response and distinct magnetic
functional are established. The attempted direct magnetic transfer is
refuted by the field map and (1). A closed positive actual mean/spin
constitutive law remains an active construction:0133 changes the actual
background through a finite-amplitude nonuniform relative equilibrium;
0134 solves alternative stationary Euler pressure cells;0135 explores
actual nonuniform finite-core spatial spin interactions. Their Euler
response, not a borrowed magnetic modulus, decides the next route.
