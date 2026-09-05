# Actual field-changing Kelvin lift of the positive fixed-tag clock

The full coadjoint calculation repairs the zero-Euler-velocity limitation
of0216. The physical tag and its clock are unchanged, but the actual
vorticity now moves with the material displacement. The inherited energy
also changes by a derived factor, rather than being copied from a label
variation.

## Complete two-scalar Kelvin representation

Use the actual C016 cell u=(psi,J grad psi), H=-Delta_perp,
psi=cos b+A cos a, A=1/100. Its full curl is -u. Write an X-independent
solenoidal displacement as

    xi=(h,J grad s), T=(J grad psi).grad.

The complete Kelvin velocity is P(xi cross omega), not an independently
chosen Euler field. Its exact components are

    w_X=Ts,
    w_perp=P_perp[(s+h)grad psi],
    -curl_perp w_perp = T(s+h)=zeta.                      (1)

The omitted intermediate term is -grad(psi s), killed by the actual
periodic Leray projection. When the planar mean vanishes,
w_perp=J grad H^-1 zeta. A nonzero planar mean is an actual constant
velocity and is retained unless the initial finite first-shell moments
are explicitly zero. The passive choice below makes the entire planar
field exactly zero and needs no mean-mode exception.

Full Euler/Lin evolution gives

    h_t=-Th,
    s_t=H^-1 T(s+h)-Ts,
    zeta_t=-T(1-H^-1)zeta.                              (2)

The planar mean-free subspace can be fixed by P1(s+h)=0, with its
stationary/forced first-shell rows treated separately. In particular
r=zeta-Ts=Th can be NONZERO on the full three-dimensional Kelvin leaf.
The r=0 condition in0216's planar continuation implicitly fixed h=0;
it is not a restriction on all admissible three-dimensional generators.
Allowing h is the material representation change that unlocks the
positive route here.

## What the zero-axial-generator route actually establishes

For h=0 and zero planar mean the actual full velocity is
(zeta,J grad H^-1 zeta), with zeta=Ts. The canonical full coadjoint forms
give

    Omega_12=rho integral s1 Ts2,
    Hphysical=-rho/2 integral zeta(1-H^-1)zeta.           (3)

Thus the norm ||(1-H^-1)^(1/2)zeta||^2 is nonnegative, but the physical
Jacobi energy of this sector has the opposite sign. The whole first
Fourier shell is its actual nullspace. Dropping that shell before the
Euler reconstruction loses its forced physical output.

The positive-weight representation is useful for the generator, not a
license for positive angular action. For clarity, the primary source
[Shvydkoy--Latushkin, Theorem3](https://arxiv.org/pdf/math-ph/0306026)
was inspected directly: it concerns essential spectrum on specified
Sobolev/vorticity spaces, with arbitrarily long trajectories. It does
not give a positive physical energy or a nonzero tagged optical output.
Its transport-plus-compact discussion likewise preserves that distinction.
No unproved spectral observation density is imported here.

Nineteen first-pass exact checks use the accepted euler_fourier full
coadjoint and material-Jacobi APIs. They retain all product harmonics,
the full axial velocity, the mean-preserving pressure projection and
the forced first shell. The calculation establishes (3) and explains
why this particular positive-weight route does not supply the requested
positive physical action. It is not a no-go for Euler or for other
Kelvin generators.

## The exact positive passive choice

Choose the actual axial displacement h=-s. Equations(1)--(2) reduce to

    xi=(-s,J grad s),
    w=(Ts,0,0),
    s_t=-Ts.                                            (4)

This is an EXACT invariant sector of the full Euler/Lin system on the
nonlinear periodic cell. Axial pressure is zero and the full transverse
Euler equation is satisfied: w is X-independent and u has no X
dependence, so the only velocity evolution is (Ts)_t=-T(Ts). No local
pressure replacement or approximate straight core is used.

The full induced vorticity curl w is nonzero for the nonconstant profiles
used in0216. Since w=P(xi cross omega), it satisfies
delta omega=curl(xi cross omega). Thus the actual vorticity and the
material flux tag are varied by the SAME volume-preserving displacement.
In particular, the first variation of omega.grad chi=0 stays zero for
the moved tag and moved vorticity. This is the vortex-transport property
that a w=0 label variation did not establish.

At every time the transverse displacement is precisely the one used
in0216, and the axial velocity does not contribute to axial mechanical
spin. The actual theta, axial S and G are therefore EXACTLY its equations
(2),(5), including G_t=S and initial G. The new axial displacement -s
is retained in all full-fluid moments and forms. Its mean vanishes for
the selected nonzero streamline harmonics. No claim that its quadratic
axial response vanishes follows from that zero mean.

## Complete positive phase and physical energy

The full canonical forms, including the axial component, are now

    Omega_12=-rho integral s1 Ts2,
    Hphysical=rho/2 integral (Ts)^2.                    (5)

The phase agrees with the homogeneous-label phase, but the physical
energy is exactly HALF the value of0216's homogeneous label family.
This follows directly from w=(Ts,0,0): its self-helicity is zero, so the
full Beltrami coadjoint energy is rho||w||^2/2. The material Jacobi
stiffness plus the actual rate Gram gives the identical result.
Seventeen first-pass canonical checks expose both forms and the full
Euler evolution. The positive route is not selected by changing the
sign convention of (3).

For s1=f(E)cos(l theta_o), s2=f(E)sin(l theta_o),

    Omega_12=-pi rho l integral f^2 dE,
    H1=H2=pi rho l^2/2 integral omega(E)f^2 dE,
    H12=0.                                             (6)

For a narrow smooth band at fixed E_* the actual observed columns tend
to (C cos(nu_*t),-C sin(nu_*t)), nu_*=l omega(E_*). The inherited mass
M=pi rho l integral f^2/(nu_* C^2) is positive, and (6) tends to the
corresponding mechanical oscillator energy M nu_*^2 C^2/2. Unlike the
homogeneous-label family, there is no factor-two energy excess here.

With l=4, the ACTUAL nonlinear-cell spin coefficient remains the positive
j_* constructed in0216. Thus this is a field-changing, positive-action,
fixed-positive-tag physical clock. It is a smooth-band finite-time
prepared construction, not a delta-function eigenstate or a fixed
finite-dimensional oscillator subspace.

## Remaining normalization and the next actual control

For a fixed tag and fixed observed angle, narrowing the band increases
M like inverse width, while j_* remains finite. Equality M=j_* has NOT
been supplied by (5). An inertia assigned after observing the clock
would miss this concrete mismatch. The physical same-field common-K,
centroid/shape and acoustic cross-action construction also remains active.

The failure-generated next control uses the actual opposite-sign Kelvin
sectors already derived here: h=-s has positive energy, whereas h=0 on
the mean-free high-mode sector has negative energy. On fixed wrapped
bands away from the observed tag, their real quadrature energy matrices
can be normalized using the full exact forms. Unequal actual transport
frequencies make their phase-to-energy ratios different. Pairing the
opposite energies may therefore supply signed phase with zero total
energy, using genuine Kelvin fields and positive ensemble weights.
Full pressure, finite cross constraints, remote physical errors and the
compatible spatial scale are the next derivation, not assumed here.
This route is registered separately as0221 and is executed immediately.
