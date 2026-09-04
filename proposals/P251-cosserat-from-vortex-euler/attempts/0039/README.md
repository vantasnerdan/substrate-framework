# 0039 — material translation and screw-stationary continuation

Continuation of 0035/0037. Exact fixed subtheorems; no comparator or numerical
remainder. The original N3 permits a slowly varying affine ensemble. Equality
with microscopic Euler at every finite wave number is not a new acceptance
condition: the conditional coarse-graining premises remain those in the frozen
proposal. The full objective still needs a derived physical angle action,
finite-core coefficients, and a compatible stationary tube ensemble.

## Material translation and angular reaction

For an advected parcel D(t), M=rho|D|, X=M^-1 integral_D rho x, and V=Xdot,
Reynolds transport and Euler give the exact identities

    M Xddot = -integral_boundary(D) p n,
    Ldot = -integral_boundary(D) (x-X) cross (p n),
    L = integral_D rho (x-X) cross (u-V).

The moving-center contribution to Ldot vanishes because integral_D rho(u-V)=0.
Thus pressure couples parcels and provides torque; a parcel is not a freely
prescribed cage. In a periodic full cell the global translation is Galilean,
as checked directly against the nonconstant-pressure ABC solution in 0037.

Let Q=integral_D rho r r^T, I=tr(Q)Id-Q, r=x-X. When I is positive definite,
Omega=I^-1 L is the least-squares angular velocity. With w=u-V-Omega cross r,
integral rho w=0 and integral rho r cross w=0. The exact energy decomposition is

    E = M |V|²/2 + L.I^-1.L/2 + integral_D rho |w|²/2.

This derives positive material translation and locked inertia without assigning
mass to a circulation centroid. Omega is an instantaneous velocity projection;
an independent angle coordinate requires the physical structure map and shape
connection, not a relabeling of Omega. Shape motion and pressure reactions
remain in the action. The verifier independently derives the tensor identity
for a general linear velocity jet and anisotropic second moments, including
the strain-induced angular momentum missed by a naive rigid-spin split.

## A genuinely stationary 3-D route

The planar pressure-absorption argument in 0037 is not available generically
in 3-D. A different construction uses screw symmetry and a Galilean boost.
For fixed h>0, let a profile v satisfy

    v(R_theta x_perp,z+h theta) = R_theta v(x_perp,z).

Differentiation gives Jv-(Jx.grad)v=h partial_z v. Consequently a rotating and
axially translating profile

    u(t,x)=R_(a t) v(R_(-a t)x_perp,z+b t)

is exactly v(x_perp,z+(ha+b)t). If this is an Euler solution, then

    w(x)=v(x)+(ha+b)e_z,    p_w(x)=p_profile(x)

is an actual stationary Euler solution, with the same vorticity. This changes
the physical uniform axial background; it is not an unsupported removal of
Coriolis terms. Pressure, material transport and kinetic coefficients must be
computed in this specified field. The source screw convention fixes the sign.

Candidate source: Averkiou–Musso, arXiv:2511.12296v1, equations (1.4)–(1.5),
Theorems 1–2, <https://arxiv.org/pdf/2511.12296>. Its polygon construction
provides a relevant finite-core helical family. The proof's positive-part
nonlinearity and precise differentiability are audited in 0036; the word
"smooth" in its abstract is not an automatic C-infinity or spectral license.
No result here imports a planar limit or a finite-core optical spectrum from
that theorem. The stationary field is generally not a constant-lambda Beltrami
field, so EPS compatibility still needs its own bridge. This is an additional
Euler construction, not a re-proof or replacement of EPS existence.

Status: in progress. Next: joint finite-core collective action, using the
material mass/reaction identities and a physical non-gauge angle sector.

Completed child receipt: verify.py exits zero with 9/9 exact checks, with
first-run output retained here. Route verdict: established for the material
mass/angular-reaction/locked-inertia identities and the stated screw-profile
stationarization. The joint angle/material action continues separately.
