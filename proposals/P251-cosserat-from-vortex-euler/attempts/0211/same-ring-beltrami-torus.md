# A literal constant-curl nonzero inner torus on the same global ring

The new seed and its full border are proved in `affine-inner-seed.md`.
Fix its small but positive delta FIRST. The large-radius construction
below uses the full global Green field of 0186/0195, not a toroidal
Dirichlet solution patched to an unrelated exterior.

## 1. The same global Green construction with the proved new seed

Use physical cylindrical coordinates (r,varphi,z), r=R+x, and phi=psi/R.
Set F(psi)=R G_sigma(psi/R) and take constant Bernoulli B. The full
steady Euler equation is

    -Delta*phi=f_delta(phi),
    u=(-psi_z/r,F(psi)/r,psi_r/r),
    curl u=G_sigma'(phi)u.                                (1)

In 0195's full stream kernel the effective source is
q_R=(R^2/r^2)f_delta(phi). Its complete bordered equation is

    phi=Ktilde_R q_R(phi)+a-U x-U x^2/(2R),
    integral q_R=Gamma_delta, integral x q_R=0.             (2)

On the fixed source disk R^2/r^2=1+O(1/R), and the full renormalized
kernel tends to K0 with O(log R/R) in the source-to-potential elliptic
norm used in 0186. The new derivative at the planar limit is exactly
the invertible map (12) of the companion proof. The small-parameter
IFT therefore constructs (2) for every sufficiently large finite R.
No unchanged Rankine inverse or continuation through an unproved
profile interval is invoked.

The extension and outer threshold are the actual full-space ones:
Logarithms use the fixed seed length a_ref as their unit, so log R
below means log(R/a_ref), with the constant part absorbed into mu_R.

    phi_global=K_R q_R-mu_R-U_R r^2/(2R),
    U_R=Gamma_delta log R/(4pi R)+O(1/R)>0,
    mu_R=3Gamma_delta log R/(8pi)+O(1)>0.                   (3)

The new source factor contributes only a source-first-moment term
to the leading logarithmic drift, fixed to zero by the center border.
Thus these coefficients are the same full-kernel coefficients as 0195,
not a supplied translation speed. Negativity on the outer source-disk
boundary, axis and infinity, followed by the exterior maximum principle,
excludes an extra positive island. The continued source in (2) is
therefore exactly f_delta of the global solution, with its displayed
R^2/r^2 factor, rather than a local source cut off after solving.

F and all its derivatives vanish on the axis and exterior, where phi
is negative. The same streamfunction gauge correction as 0186 gives
a regular global vector potential. The velocity and vorticity are
C-infinity, with compact toroidal vorticity and the ACTUAL uniform
far velocity -U_R e_z. The steady-frame action retains E-U_R I_z.
Adding that uniform velocity to move to the rest-at-infinity frame
would change a constant-curl identity; no such frame change is hidden.

## 2. An actual open constant-lambda region and nonzero closed core

Choose a small fixed inner radial disk s<s_b in the planar seed,
with phi>delta by a positive margin and J0(lambda s)>0. Smooth
convergence of the global construction preserves this margin on a
slightly larger disk about its unique elliptic meridional center.
There, EXACTLY, not just at leading order,

    G_sigma'=sigma lambda,
    curl u=sigma lambda u.                                (4)

In particular the actual toroidal center speed is nonzero: it converges
to W_sigma(0)=sigma lambda A. The center revolves into a physical
closed circular streamline and vortex line on this SAME global field.
The nearby level curves of phi are smooth nested closed curves; their
revolutions are actual invariant tori. Their entire positive-radius
inner neighborhood, including the core, lies in (4). This is not the
quiet cavity of 0200 or merely a generalized-curl factor on that region.

The global factor G_sigma'(phi) is generally variable in the outer
taper and zero outside the core. A globally constant-lambda compact-
vorticity field is not claimed. The positive result is a literal local
Beltrami torus naturally embedded in a global smooth Euler solution.

## 3. Full transit, flux action and a strict twist margin

It suffices to take sigma=+1; the meridional reflection sigma=-1 gives
the opposite swirl/curl sign and the same geometric conclusions with
the corresponding flow-oriented return. Orient the poloidal angle
along the flow. The straight inner limit has

    V(s)=lambda A J1(lambda s), W(s)=lambda A J0(lambda s),
    F0(s)=V(s)/(s W(s))=J1(lambda s)/(s J0(lambda s)),
    J0_flux(s)=integral_0^s W(t)t dt=A s J1(lambda s).        (5)

J0_flux is the ACTUAL toroidal-section flux action, not unweighted
meridional area. Expanding the Bessel equation gives

    F0=lambda/2+lambda^3 s^2/16+lambda^5 s^4/96+O(s^6),
    dF0/dJ0_flux=[lambda^3/8+O(lambda^5 s^2)]/W(s)>0        (6)

on a sufficiently small fixed positive-radius annulus. Its central
limit is lambda^3/(8W(0))=lambda^2/(8A)>0. First fix such an
annulus with a strict positive minimum; only then increase R.

For the actual ring, on a regular level C={phi=h}, the toroidal
advance per poloidal circuit and the flux action are exactly

    Delta_varphi_R(h)=G(h) integral_C dl/[r|grad phi|],
    rho_R(h)=2pi/Delta_varphi_R(h),
    J_R(h)=(1/(2pi)) integral_inside(C) u_varphi dr dz.      (7)

These include nonuniform transit, metric and toroidal velocity. They
satisfy the exact coarea identity `J_R'(h)rho_R(h)=-R` in this
stream-value normalization. Identifying nearby
level curves with the bounded normal field grad phi/|grad phi|^2
and differentiating the integrals proves convergence, including the
required first action derivative,

    rho_R/R -> F0, J_R -> J0_flux,
    (1/R)d rho_R/dJ_R -> dF0/dJ0_flux>0.                   (8)

All contour denominators are bounded away from zero on the fixed
annulus. This is the exact flux twist required for invariant-circle
persistence, not an unweighted frequency derivative selected numerically.

## 4. Arithmetic-free radius and boundary selection

At the actual meridional center (r_c,z_c), the positive poloidal
linear frequency is (R/r_c)sqrt(det Hess phi), while the toroidal
angular velocity is R G(phi_c)/r_c^2. Hence the actual core return
rotation is

    rho_core(R)=r_c sqrt(det Hess phi_c)/G(phi_c)
                         =lambda R/2+O(log R).             (9)

The fixed new-profile bordered inverse is uniformly bounded for large R.
Differentiating the exact kernel and source factor gives
partial_R(phi_R-phi0)=O(log R/R^2) in each required finite core norm,
with the center derivative obtained by its own nondegenerate Hessian.
Consequently

    d_R rho_core=lambda/2+O(log R/R) !=0.                   (10)

This uses the actual continuous radius, not an integer harmonic or an
unproved irrationality property of a Bessel zero. Choose a sufficiently
large finite R with a Diophantine core rotation, or simply separated
from half-integers. The core's transverse multipliers are then elliptic
and different from both +1 and -1. At this now fixed R, nonzero twist
(8) lets one choose an inner boundary flux action with a Diophantine
rotation. The same chosen tube lies wholly inside (4).

## 5. Exact EPS persistence hypothesis, on this same inner region

On the affine region the scalar equation (1) has analytic coefficients
and affine nonlinearity, so its solution and the velocity are analytic
there by interior elliptic regularity. The boundary return is analytically
conjugate to its rigid rotation by its flux action-angle coordinates.
The actual positive section measure is u_varphi dr dz. Its normal
torsion is nonzero by (8), and its boundary rotation was selected
Diophantine. These are precisely the hypotheses of EPS Theorem 7.6,
not those of its unrelated weak-lambda source construction.

The primary theorem and the actual nearby-measure repair were reread
in `../../sources/1210.6271.pdf.txt`, lines 2825 onward and 3220-3330:
[Enciso--Peralta-Salas, Theorem 7.6 and equations (7.28)-(7.32)](https://arxiv.org/pdf/1210.6271).
For nearby analytic divergence-free fields on a slightly larger tube,
the section measure changes. The local Moser identification of those
positive analytic two-forms is made BEFORE applying the same-measure
theorem. It supplies a nearby invariant boundary. The nonunit core
multipliers separately give its nearby periodic core by IFT, and its
strict ellipticity persists. A Diophantine core is not claimed to remain
Diophantine under every perturbation; the robust conclusions are the
nondegenerate elliptic core and the KAM boundary.

Thus the SAME smooth global Euler ring has a literal constant-curl,
nonzero closed core and a robust inner unknotted vortex tube. No
global replacement field is needed to earn this local EPS-compatible
geometry. A later global Beltrami approximation would be a different
background and would need its own dynamical/action transfer. This
child has not identified the new Bessel-core modes, material clocks or
coupled continuum merely by establishing their correct geometric home.
