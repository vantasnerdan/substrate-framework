# Actual joint physical positions, rate chart, mass and retained current

## 1. The apparent initial position bracket cancels with the true means

Consider one real transverse/helicity component, with complete real
Bloch partners retained when necessary. Let z=(D,V,z1,z2), and let the
optical initial generator and actual Euler-velocity means be m_a,n_a.
The actual joint symplectic matrix is

    Ω=[[rho J,C],[-C^T,beta J]],
    C=[[rho n1,rho n2],[-rho m1,-rho m2]].               (1)

The ACTUAL initial full-fluid position and velocity rows are

    X0=(1,0,m1,m2)z,       Xdot0=(0,1,n1,n2)z.         (2)

As usual {F,G}=−F Ω^(-1)G^T in the convention Ω(D,V)=rho.
The V column of Ω is rho X0^T and the D column is −rho Xdot0^T.
Thus, whenever the joint phase is nondegenerate,

    {X0,Xdot0}=1/rho,
    {X0,theta0}=theta_V(0)/rho=0.                      (3)

The last equality follows because the common-V column has ZERO initial
material displacement. The actual central material angle depends on
the displacement, not its cotangent momentum, at that instant.
Dropping m,n in (2) while keeping their cross form in (1) creates a
false position-bracket defect. Conversely dropping the cross form while
keeping the means also changes the actual physical bracket.

The Schur denominator is beta+rho det([n;m]); its nonvanishing,
not separate positivity of two diagonal entries, decides the initial
joint phase rank. The compact reaction construction in the companion
file sets m=n=0 exactly using actual exterior Euler generators. This
gives an actual block-orthogonal phase before the physical chart is
chosen. Its small changes to beta and to physical tag rows are measured
on the packet's own scale.

## 2. Exact initial physical rate chart and its positive mass

After that construction choose optical INITIAL coordinates (q,r) as
the actual tag angle and its actual initial rate on the optical plane.
This is the physical 0147 chart, not a newly named unobserved phase.
Let its positive scalar mass at the reference instant be j, including
the declared cell-volume density factor. Then Ω_opt=j J.

Write b for the ACTUAL initial tag-angle response to acoustic D and
d for its actual initial time derivative. The two full physical rows
and their rates are

    Q=(X,theta)=(D,q+bD),
    R=(Xdot,thetadot)=(V,r+bV+dD).                       (4)

The coefficient b in the common-V angle rate is the same measured b:
eta_V(0)=0 and eta_V,t(0)=V. There is no fitted angular response.
The chart [Q;R] is invertible for every b,d. Direct use of the joint
form gives

    {Q,Q}=0,
    N={Q,R}=[[1/rho,b/rho],[b/rho,1/j+b²/rho]],
    M=N^(-1)=[[rho+j b²,−j b],[-j b,j]],
    det M=rho j>0.                                     (5)

Thus the ACTUAL initial physical mass is positive without assuming
the tag ignores affine strain. The initial position-position block
of the pulled-back symplectic form is j d J, and its rate-rate block
is zero. The d connection is retained. On a coherently time-reversed
population d=0 at the symmetric initial time, but that does not set
later connections to zero.

The scalar density normalization can be made explicit for a whole-field
population of axial packet marks. Prepare each local angle with the
coherent input n·Phi and reconstruct the physical vector as
(E[n⊗n])^(-1)E[n theta_n]=3E[n theta_n]. The actual quadratic action,
tag spin and initial G are averaged BEFORE this reconstruction is used:
each carries the factor E[n⊗n]=I/3. Thus one identical packet per cell
has j_macro=j_packet/(3 Vcell), and its literal spin/G laws have that
same coefficient when their packet laws match. The reconstructed angle
is not the raw mean, and rho is still the complete ambient-fluid mass.
For nonidentical marks use their actual covariance/measure in place of
these displayed identical-packet factors.

For the actual macro displacement D exp(iK·x), the local first gradient
is iD⊗K. Applying 0148's directly differentiated central-quadrupole
formula and averaging the complete in-plane material mark gives the
rigid-rotation coefficient one: the initial axial angle is
n·curl D/2. Reconstructing the vector from E[n⊗n]=I/3 gives b=h/2
in a transverse helicity sector curl=h, up to the actual higher spatial
jet. This is an initial field/observation calculation, not a claim that
the acoustic tag follows an affine history for all t.

The new in-tag displacement-moment control gives the actual initial
hybrid map U=X−j h q/(2rho) at first gradient, with the acoustic
displacement/shape second-gradient terms still retained. Pulling back
(5), or equivalently the diagonal (X,q) kinetic form, gives leading
mixed kinetic coefficient [rho*j/(2rho)−j/2]h=0. This is the genuine
initial physical normalization previously unavailable from spin alone.
It does not supply the later optical spatial curvature or an autonomous
locking potential by matrix naming.

## 3. Full-time physical chart: exact conditions and complete action

Propagate ALL joint columns by the same actual linearized Euler/Lin
operator on the same stationary field. Average their actual phase
actions on common initial coefficients before eliminating anything.
The conserved resulting Ω is the computed full joint matrix. In these
initial coordinates the effective Hamiltonian is zero, as proved in
0154. Let Q(t) be the actual physical (X,theta) observation rows and
let R=Q_t; include all tag registration, frame and current derivatives.
Define

    T=[Q;R],      B_Q=−Q Ω^(-1)Q^T,
    N=−Q Ω^(-1)R^T.                                    (6)

The exact rate chart exists where det T is nonzero. It has

    Ω_y=T^(-T)Ω T^(-1),
    B_y=T_t T^(-1),
    H_y=−Ω_y B_y−(Ω_y)_t/2.                            (7)

The complete first-order action is −y^TΩ_y y_t/2−y^TH_y y/2.
No connection or complementary observation is removed. Equations
(6)-(7) apply also to the full realified Bloch/multipolar matrix.

There is a precise additional condition for the ordinary position/rate
mechanical action: the physical configuration rows must satisfy
B_Q(t)=0 throughout the window. Then N is symmetric, because
B_Q,t=N−N^T. If N is positive definite, M=N^(-1)>0 and the rate-rate
block of Ω_y vanishes. Writing its position-position block as A=-A^T,
the action can be expressed, with its actual time boundary, as

    L=Q_t^T M Q_t/2+Q_t^T A Q/2−Q^T K Q/2,
    M Q_tt+(M_t+A)Q_t+(K+A_t/2)Q=0.                   (8)

Here K is symmetric and determined by (7); it is not independently
inserted from the accepted conditional C009 pencil. If B_Q is nonzero,
the actual rate-rate symplectic term survives and a variation in the
rate contains its derivative. Retain (7); dropping that term would not
be a valid algebraic Legendre elimination in the registered positions.

The exact moment conditions at finite time can be written without
unknown canonical coordinates. In a block-orthogonal acoustic/optical
initial phase, let x_a,x_o be the actual X rows and t_a,t_o the actual
theta rows. For a single optical pair of density form beta J,

    {X,theta}=det([x_a;t_a])/rho
                          +det([x_o;t_o])/beta.             (9)

This is the genuine observed-response condition. Initial orthogonality,
separate scalar masses and matching tagged spin do not by themselves
prove it at later times. The actual whole-field mean optical response
starts at second spatial order under the axial/isotropic preparation
and the repaired zero means. The corresponding acoustic tag row still
requires its own actual local response; a small BULK L2 comparison
cannot be divided by its small quadrupole to replace that row.

At K=0, the acoustic columns are exact material translations and
Galilean translations: eta_D=D, eta_V=tV. Their registered central
tag-angle rows vanish; the repaired optical full mean remains zero by
momentum conservation. Thus T is block diagonal at K=0, with the
positive actual acoustic and optical scalar charts on their declared
finite windows. Continuity at a SELECTED finite cell gives a regular
full joint rate chart for sufficiently small K. It does not prove
B_Q=0 at nonzero K, and does not make an inverse-period estimate uniform
as the cell changes. The full action (7) is the established retained
state at that point; (9) identifies the next actual reduction condition.

## 4. Hybrid positions and the dilute optical error scale

Use the ACTUAL initial row and integrated current of the companion
file, and denote its accumulated angular row by
H(t)=G0+∫_0^t S(s)ds. Then, to the proved first gradient,

    U=X−h H/(2rho)+[symmetric-shape/transport rows],
    {U,theta}={X,theta}−h{H,theta}/(2rho)
                                +[shape/transport brackets]. (10)

If S=I_s theta_t+chi_s theta+e, the exact decomposition is

    H=I_s theta+[G0−I_s(0)theta0]
                       +∫(chi_s−I_s,t)theta+∫e.             (11)

The new initial tag control fixes the square bracket at the selected
normalization. The term proportional to theta commutes with theta;
the integral connection and real spin error do not automatically do
so. Their 0147 relative bounds remain part of (10). This prevents a
coordinate-only shift of X from being mistaken for a correction of a
nonzero physical-position bracket.

For every kinetic or phase comparison use the natural block scale
D=diag(rho,j). A sufficient positive-mass perturbation bound is

    ||D^(-1/2) delta M D^(-1/2)||<epsilon<1.            (12)

In particular mixed errors must be smaller than
epsilon sqrt(rho j), and the optical diagonal error smaller than
epsilon j. The equivalent inverse-mass condition uses
D^(1/2) delta N D^(1/2). For beta_packet and all physical spin/G
moments, use packet norms first and divide every action and observation
by the SAME Vcell. The relative margins therefore do not vanish merely
because a selected finite construction is dilute. Nothing here asserts
a fixed positive j limit while Vcell tends to infinity.

## 5. Strongest construction and next executable join

Established together: actual initial cross rows; exact initial canonical
mean subpair; constructive compact exterior reaction moments; an extra
in-tag initial displacement control preserving the earlier spin match;
positive complete initial physical mass; and the exact full-time joint
moving action and configuration-bracket conditions in the actual mean,
registered material angle and hybrid-current observations.

The remaining full-time two-position reduction is (9)-(11), with the
actual acoustic tag response and the retained currents. Its next routes
are a controlled same-cell affine/rotational acoustic tag history or
additional actual response/moment columns that close the physical row
without replacing it. The retained joint phase (7) is already a genuine
same-action construction. Neither this open reduction nor the unassumed
0155 optical long-wave license is a verdict against the positive packet,
accepted conditional C009 or the parent objective.
