# Transfer to the same global generalized-force-free ring

## 1. Order and actual operator

First fix the smooth source and a sufficiently small NONZERO axial
carrier k from the positive principal branch in
`radial-mode-and-action.md`. Its all-poloidal eigenvalue has a fixed
isolated contour, a nonzero KKS value and a positive phase energy.
Only afterward let the actual ring radius R tend to infinity and
choose integers n_R with n_R/R tending to k. These are the global
0195 force-free rings, not superposed columns.

Use their actual smooth flux coordinates (I,theta,varphi), including
the finite-axis regularity condition. The full compact-active-vorticity
linear Euler generator, including EVERY poloidal Fourier channel, is

    L_R=T_R+K_R,
    T_R=-Omega_R partial_theta-i n_R V_R+N_R,
    K_R=alpha_R partial_theta B_R+i n_R beta_R B_R
             -(alpha_R',beta_R') (B_R)_I.          (1)

Here N_R has the two real transport shears in its first column, N_R²=0;
B_R is the COMPLETE free-space inverse curl in the toroidal harmonic.
There is no reflecting wall or suppressed exterior pressure. The
straight limit has the corresponding W, Omega and both vorticity
components. The active-support domain is the isovortical compact-source
space used in 0195; arbitrary exterior passive vorticity is not silently
included as an isolated spectral sector.

## 2. Why the new derivative term still has a compact Fredholm limit

The 0186 toroidal-kernel estimate gives B_R to B_infinity convergence
in L² operator norm on the fixed active cross-section, and a uniform
one-derivative gain locally. Pullbacks and the volume weight converge
in each fixed smooth norm by the actual bordered ring construction.
The new alpha partial_theta B term cannot simply be called compact:
before a transport resolvent it has order zero.

Instead write the spectral equation as

    [1-(z-T_R)^(-1)K_R] w=0.                       (2)

On the chosen contour the full transport bands remain separated. In
poloidal channel m the explicit inverse is a scalar shifted denominator
plus N_R divided by its square. It follows, uniformly in m and R, that
both (z-T_R)^(-1) and partial_theta(z-T_R)^(-1) are bounded on L².
All the coefficients in (1) depend only on I; the transport inverse
commutes with partial_theta. A coefficient moved across the triangular
matrix inverse contributes only another bounded, I-dependent matrix.
Thus (2)'s new term can be written as a bounded order-zero transport
multiplier applied to the compact B_R. It is compact. This argument
uses the vorticity shear rather than omitting it.

The same explicit denominators prove norm convergence of these bounded
multipliers when Omega_R, n_R V_R and N_R converge: for large |m| the
denominators grow linearly in |m| and cancel the one numerator derivative;
the finitely many remaining channels use the fixed spectral separation.
Together with the actual B_R kernel convergence this proves uniform
operator-norm convergence of the analytic Fredholm families on the
contour. No separate unsupported norm estimate for partial_theta B_R
is needed. Uniform smoothing follows after this factorization because
B_R gains one derivative and the bounded transport multiplier has the
usual fixed Sobolev bounds on this noncritical contour.

The simple isolated straight pole therefore continues to a unique
actual ring pole. Real Hamiltonian symmetry and simplicity keep its
frequency real; its nonzero KKS and positive quadratic phase energy
persist. Eigenfunctions converge in every fixed local smooth norm by
the equation and the above derivative gain. This is a full-poloidal
pressure continuation at fixed k, not an m=0 truncation of the ring.

## 3. The measured section rows, and what the limit does not identify

The actual ring mode has nonzero toroidal harmonic n_R. Take a fixed
meridional section and a finite smooth transverse observation window.
The physical angle is the Euclidean covariance of an actual advected
nonnegative marker, initialized by the column marker in its local
Euclidean frame. Its full passive transport is used at finite R.
The spin is the complete-fluid relative Euclidean angular moment in
that window, with the advective and pressure torque flux through every
boundary, not the angular moment of the marker alone.

For any fixed optical interval T and observation buffer L, the smooth
background/mode convergence just proved gives convergence of their
linear angle and spin observation ROWS, and of the inherited phase
action. This is a statement about both initial phase columns uniformly
on [0,T], not division by an instantaneous angle/rate which can vanish.
Then let L grow: the column's true k>0 Bessel tail makes its section
rows converge exponentially, with the flux retained at each finite L.
Thus the actual global generalized-force-free ring supplies the positive
angle-action and positive section spin/rate overlap of the column with
arbitrarily small, explicitly ordered finite-window observation error.
At finite R the transported marker clock/current connection is the
actual one; equality (6) of the column companion is not asserted
exactly for all time on the curved ring.

These rings have the twisted material/vorticity tori and nonzero closed
core of 0195. Their force-free factor is variable G'(phi), not constant
lambda. This establishes an actual global geometric supplier with a
positive full Euler pole and measured section response; it does not
assert literal constant-lambda EPS membership.

Nor does fixed k and n_R~kR establish the n=1 or n=2 global Euclidean
centroid/covariance response. Full axisymmetric-cell tensor selection
retains only particular low toroidal harmonics in those complete-cell
moments. That next geometric route has the joint limit k=n/R and needs
a uniform low-carrier ring operator/current estimate. Three neighboring
large integers do not substitute for it or for a common continuous
laboratory Bloch wavevector. This is the next candidate, not a defect
in the fixed-sector positive pole proved here.
