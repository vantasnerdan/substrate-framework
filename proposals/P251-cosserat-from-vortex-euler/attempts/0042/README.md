# 0042 — same-fluid body angle and microscopic relative-angle action

Parent P251 / issue #198, original conditional affine smooth-Euler objective.
Owner: Codex `/root/construction_review`; write surface: this directory only.
Attempt 0040 and its positive orbit action are frozen inputs, not reinterpreted.

Positive child deliverable: compute the common material-body angle response
of the same fluid cell as 0040, including its angular-momentum connection;
combine it with the established internal angle/shape Routh action, and obtain
a covariant absolute microrotation field map without postulated inertia.

Candidate A: factor the material map as g=X+R_B h on a transported finite
representative volume. Derive the exact kinetic identity, moving-boundary
pressure work, axial locked inertia, and first internal derivatives of angular
momentum from the actual field/generators. Candidate B, if that construction
does not license the common coordinate: add a mean-sector generator and
conjugate shape to the Euler orbit and compute the enlarged KKS/Hessian.
Selection is physical field identification, same-fluid action and pressure
reaction, parameter economy, and positive coupled inertia. No comparator or
target modulus is inspected or fitted.

Sourced inputs: 0039 material center/locked-inertia identity; 0040 positive
orbit sector of u0=(−b sin y,2b sin x,2b cos x+b cos y), dimensionless
periodic coordinates with physical scale ell. xi_q and xi_s retain their
cos z/sin z cell shapes. The finite material cell is the initial cube
[−pi ell,pi ell]^3; its boundary is transported, and its Euler pressure
tractions remain in the action exchange with adjoining cells. A fixed torus
does not acquire a continuous physical rotation merely by naming B.

Analytic design: start from physical kinetic energy, keep the complete
Omega·L_internal term, define the momentum ensemble before Routh reduction,
and evaluate I_body and the internal linear momentum/shape connection rather
than assuming them zero. Distinguish macro displacement evaluated at cell
centers from a continuum velocity field evaluated at every microscopic point,
so the affine rotation kinetic energy is counted exactly once. All present
checks are symbolic integrals and variations; there is no numerical remainder.

Pass licenses: the computed material-cell kinetic/action map with its explicit
affine and moving-boundary premises. The parent carries its derived stress,
wryness and EPS tube geometry through their remaining obligations.

## Exact result and route verdict

`verify_body_map.py` derives the kinetic decomposition and the complete
conditional change of fields. `stdout.txt` records its first successful
execution, 22/22 exact checks. At the stated initial cube,

    I_B/V=2 rho pi² ell²/3,  L0/V=3 rho b ell,
    I_q/V=6 rho ell²/25,
    B_body=beta+q,  Psi=beta+A q,  A=1+I_q/I_B,
    J_Psi=I_B²/(I_B+I_q),  J_beta=I_B I_q/(I_B+I_q),
    K_Psi=K_q/A².

The original q/s pair has zero axial-average first variations of angular
momentum and locked inertia, and zero averaged axial mechanical curvature.
These identities are real useful field-map results. They do not themselves
identify the sum of the locked material metric and the orbit-reduced metric
as one Euler action. The verifier exposes two concrete gaps in that glue:

1. The periodic cube is not an invariant material parcel: on x=±pi,
   u_x=-b sin y. Its transported off-diagonal moment immediately changes,
   `Qxy_dot/V=rho ell b`. A stationary Eulerian cube metric cannot be used
   as the constant metric of that advected parcel.
2. For a time-dependent physical rotation, `g=R_B g0` has velocity
   `R_B u0 + Bdot ez cross x`, with vorticity `R_B omega0+2 Bdot ez`.
   A rotation of the fixed isovortical leaf pushes vorticity to
   `R_B omega0`. Consequently simply adding the body's material locked
   inertia to the fixed-orbit Routh action mixes two momentum ensembles.
   The missing reaction is not a numerical error in the field map.

**Candidate A route verdict: refuted as a direct same-orbit construction.**
The exact material kinetic identity and conditional covariant map remain
established; the naive sum is not an Euler-derived common-angle action.
This verdict does not close the common-angle obligation.

## Executed boundary repair and representation continuation

The same prototype supplies an exact invariant material parcel:

    psi=2b cos x+b cos y,
    D={psi>2b} cross S1_z.

Because u.grad(psi)=0 and the flow is incompressible, D is preserved exactly.
Writing `Y(x)=acos(2-2 cos x)`, its transverse range is
`0<=x<=pi/3`, `|y|<Y(|x|)`. The axial locked inertia is constant and is the
positive same-fluid integral

    I_D/V_D = rho ell²
      integral_0^(pi/3)[x²Y+Y³/3] dx / integral_0^(pi/3)Y dx.

No reset of a transported parcel, periodic-face flux omission, or prescribed
inertia enters this repair. The stationarity of every spatial moment over D
follows from Reynolds' formula: u.n=0 on its sidewall and the axial faces
cancel. Positivity of the integral is immediate from Y>0 in the interior.

To keep generators inside this parcel, replace the global trigonometric pair
by compact curl-generated circular polarizations. Attempt 0045 carries this
repair farther: it proves, by an explicit finite-wave-number bound, a full
positive core/cage H and nondegenerate KKS pairing directly inside *any*
actual EPS material tube. This is an executed representation change, not a
claim that the old periodic matrices survive a cutoff.

For a z-independent compact cage bump phi the exact new generator cross is
`phi² ez+(phi/k)J grad_perp phi`. Its averaged axial mechanical curvature is
`2 rho integral(phi²)/I_D`, which is nonzero. Thus the old zero-curvature
calculation cannot be copied to this repaired parcel. One must retain this
same-fluid momentum connection in a full reconstruction, or construct the
common angle and its conjugate directly on the Euler orbit.

Candidate B is being executed independently in new attempt 0046: an
SO(2)-invariant material cylindrical Beltrami cell supports the actual
common-rotation symmetry direction, and a compact angular-momentum-changing
partner supplies its intrinsic Routh inertia. This changes representation
and computes the reaction on the same orbit instead of adding material mass.
The common-angle parent obligation stays active until that construction and
its coupling close. This attempt is preserved without editing 0040.
