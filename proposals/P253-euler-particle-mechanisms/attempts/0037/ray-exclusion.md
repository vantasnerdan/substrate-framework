# Exact one-pass geometry excludes a high-frequency exponential mechanism

The low-frequency block is not the only possible source of instability. This
section checks the opposite limit directly on the actual solitary base flow.

## 1. Every material ray crosses the excitation once

In the solitary frame the exact base velocity is

    u_c=V(r)e_theta-c e_z+delta u_mu(r,z),

where `0027` proves

    (u_c)_z=-c+partial_r f/r<0.                         (1)

For sufficiently small `mu`, `(u_c)_z<=-c0/2`. Therefore `z(t)` is strictly
decreasing on every trajectory and no material ray is trapped in, or returns
to, the localized excitation.

The weighted high-regularity convergence in `0027` gives

    integral_R sup_r |grad delta u_mu(r,z)| dz
      <=C mu L_mu.                                    (2)

The exterior term is included through its exact Bessel minimizer; high radial
regularity and the one-`X` weight turn the `Z^s` bound into (2). Combining
(1) and (2),

    sup_ray integral_R |grad delta u_mu(x(t))| dt
      <=C mu L_mu ->0.                                (3)

This is an actual one-pass strain bound, not compact support assigned to the
velocity.

## 2. Background and perturbed ray exponents

For the column part `V(r)e_theta-c e_z`, trajectories have constant radius,
uniform axial drift, and angular rate `Omega(r)=V(r)/r`. In a co-rotating
cylindrical frame its deformation gradient is a shear with entries at most
linear in time, produced by `Omega'(r)`. The cotangent and full-pressure
bicharacteristic-amplitude systems consequently have zero exponential
Lyapunov exponent; this is the ray version of the exact positive column metric
in `0030`.

Variation of constants through the single crossing, with (3), gives an
integrable coefficient perturbation of the two polynomial column
propagators. This proves zero exponential rate. A literal finite scattering
matrix between Jordan/shear background propagators would additionally require
a polynomially conjugated or weighted integrability estimate, which has not
yet been proved. For each fixed sufficiently small solitary member,

    limsup_|t|->infinity (1/|t|) log ||BAS_mu(t)||=0    (4)

on every ray. The polynomial prefactor may depend on `mu` and inherits the
column shear. Equation (4) claims neither convergence of a scattering matrix,
a uniform bound, nor nonlinear stability.

Thus the actual solitary wave has no positive high-frequency ray exponent.
Any exponential linear instability must instead arise from a finite/low
frequency global mode or from a nonlocal coupling not represented by the
principal ray cocycle. In particular, the compact Gavrilov repeated-orbit
mechanism established in `0032/0039` cannot transfer: its essential feature is
a returned material/covector orbit, while (1) forbids a return here.

## 3. Remaining spectral interval

The exact scalar threshold limit in `threshold-reduction.md` excludes a
right-half-plane eigenvalue in the limiting critical block through its
one-negative-direction constrained coercivity. Equation (4) excludes a
positive exponential rate at arbitrarily high frequency. The unclosed region
is therefore concrete:

1. uniform convergence of the critical Euler Evans/resolvent family to the
   KdV family near spectral scale `mu/L_mu`;
2. absence or classification of intermediate-frequency point spectrum; and
3. the `Pi/Q` derivative-form coupling in
   `scattering-block-reduction.md`.

Route C now targets this bounded spectral interval. A finite-dimensional
Galerkin matrix cannot decide it. An analytic Fredholm/Evans family with the
exact Bessel exterior can; if its determinant leaves a small sign or root
location irreducible, that scalar will receive a separately preregistered
small-ratio computation rather than an unbudgeted stability simulation.
