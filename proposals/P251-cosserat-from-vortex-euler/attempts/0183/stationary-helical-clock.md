# An exact finite-core stationary helical material observation

This route changes the geometry responsible for the painted-tag dephasing
in0176. It establishes a stationary Euler field and material observation,
not an Euler eigenmode or an EPS continuum by itself.

## Exact Cartesian field and its admissible localization

Write s=x²+y², h=(-y,x,c), c>0 and u=f(s)h. Direct differentiation gives

    div u=0, (u.grad)u=-f(s)²(x,y,0),
    curl u=(2cy f_s,-2cx f_s,2f+2s f_s).

Therefore EVERY smooth radial f solves stationary incompressible Euler
with pressure per density p_s=f²/2. There is no radial material flow.
The force-free condition curl u=lambda(s)u requires

    f+(s+c²)f_s=0,
    f=C/(s+c²), lambda=2c/(s+c²),
    p=-C²/[2(s+c²)]=-|u|²/2.                         (1)

This is smooth at the axis and has an open finite-radius core. Its curl
factor is VARIABLE. It cannot be passed to a constant-lambda EPS theorem
or a constant-lambda canonical API without a new construction.

For a globally smooth radial velocity localization, choose a smooth
compactly supported f(s) agreeing with (1) on the chosen core. The actual
pressure is

    p(s)=-1/2 integral_s^infinity f(t)² dt.           (2)

The Euler residual still vanishes exactly, while the vorticity changes
by the displayed f_s terms. Force-free compatibility is lost in the
cutoff annulus: the first-order equation in(1) has no nonzero compactly
supported solution. This is an explicit velocity localization, not a
Biot–Savart truncation claimed to preserve the old vorticity.

The localized column is periodic in z with any chosen period, hence has
finite energy on its periodic cylinder and on a transverse periodic
array with disjoint supports. Smooth zero velocity outside each support
makes such an array an exact Euler solution; the pressures also match
their zero exterior values. It is not compact in all of R³, and its
helical trajectories on a periodic cell are not a constructed Euclidean
closed EPS torus. Averaging whole arrays under translations and O(3)
gives a stationary isotropic probability law, not interaction stiffness.

## One invariant smooth material clock on an open core

For each integer m>=1 define the globally smooth complex scalar

    F_m=(x+i y)^m exp(-i m z/c).

It obeys h.grad F_m=0 and hence u.grad F_m=0 for EVERY f above. The
angular shorthand theta-z/c describes this invariant away from the axis;
F_m is the actual regular observable used at the axis. On an axial
period 2pi c take

    chi=g(s)[1+b Re(F_m/a^m)],
    0<=g<=1/(1+|b|), support(g) subset {s<a²}, 0<|b|<1.

Then 0<=chi<=1 and u.grad chi=0. Its finite mass fraction and all its
reference moments are stationary. Angular or axial phase integration gives

    Q0=integral chi F_m dV
      = b/(2a^m) integral g(s) s^m dV !=0.           (3)

Choose b>0 for a positive real reference. Under a physical rotation about
the column axis by phi, Q becomes exp(i m phi)Q. Thus arg(Q/Q0)/m is a
unit-normalized small material angle. Axial translation changes the same
relative helical angle; it is an observation of the registered material
helix, not an absolute vorticity director or an independent axial rotor.

## Actual material perturbation and current

For a solenoidal Euler velocity perturbation v, let the tag perturbation
obey its actual passive-density equation

    delta_chi_t+u.grad delta_chi=-v.grad chi.

Use periodic axial boundary conditions and compact transverse tag support.
Integration by parts, div u=div v=0, and u.grad F_m=0 yield exactly

    delta_Q_t=integral chi v.grad F_m dV.            (4)

Equivalently, for a Lin displacement delta_chi=-xi.grad chi,
delta_Q=integral chi xi.grad F_m. The terms involving [u,xi] in its time
derivative cancel by the invariant-scalar identity. No choice of a
corotating time coordinate and no radial dephasing approximation is used.

If an ACTUAL Euler mode v=exp(-i omega t) vhat with omega!=0 exists,
(4) has a monochromatic particular solution. Initial tag displacement
must match that solution; an arbitrary initial tag also leaves a constant
offset. This implication is exact but does not supply the mode, its
Hamiltonian sign, its curvature, or its spectral persistence in EPS.
For a real mode take real and imaginary parts of the complexified result.

The current in(4) is a weighted helical moment. In particular
partial_z F_m=-i m F_m/c, so axial perturbation enters it. It is NOT
automatically the mechanical angular momentum integral rho chi r cross v
nor the full displaced-domain spin. A positive canonical action and a
nonzero physical spin overlap still require derivation from that same
Euler mode. This distinction preserves the actual0181 current achievement.

## Actual helical operator interface and next construction

For equivariant vector perturbations [h,v]=0, Cartesian basis rotation
gives h.grad v=Jv, Jv=(-v_y,v_x,0). Thus full linear Euler reduces to

    v_t=-P[2 f Jv+2 f_s(x v_x+y v_y)h].             (5)

The Leray pressure and the actual domain remain in P. Dropping the first
term because the scalar helical phase is invariant would erase a genuine
Coriolis contribution. The screw Killing momentum tau=u.h obeys material
transport under helically symmetric Euler evolution, so its perturbation
satisfies delta_tau_t+v.grad tau=0. For (1), tau=C and any nonzero-frequency
helical mode has delta_tau=0. For localized f, tau=(s+c²)f is variable and
this reduction changes. Attempt0185 derives the remaining quotient
vorticity, phase/action and spectral problem from these actual equations.

No magnetic evolution or stability theorem is imported. The algebraic
profile is also known as the Gold–Hoyle force-free profile in the magnetic
literature; source provenance is Solov'ev and Kirichek (2021), section4.2,
https://academic.oup.com/mnras/article/505/3/4406/6291197 . Only the profile
identification is taken from that source; equations(1)–(5) are derived here.

route_verdict: established as stated (stationary Euler/material-clock
construction and exact observation/operator identities).
evidence_scope: exact analytic geometry and conditional actual-mode
observation; not a spectral existence or parent-continuum theorem.
The next positive obligations are executed in0185 and0186, while0184
continues the same-cell coupled physical-current closure.
