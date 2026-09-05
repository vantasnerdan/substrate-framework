# Whole-field reflection and an exact threshold response bound

This executes the failure-generated parity repair registered in the README.
It is a new positive operator estimate, not a relabeling of the earlier
embedded branch or an assumed full-ring eigenvalue.

## 1. Actual reflected field and one-sided transport

Reflect the complete global 0195 field in a meridional plane. Euler is
invariant under this orthogonal transformation. Its meridional velocity,
translation, source f=G G' and torus geometry are unchanged, while the
toroidal velocity changes sign. Vorticity is transformed as an axial
vector, so the generalized force-free factor changes sign as required;
it is not held fixed while reflecting only selected terms.

In the straight core write W=-w, w>=0, and Omega>0. Then

    Z=2Omega+r Omega', Kappa²=2Omega Z,
    w w'=-r Kappa²/2, Z=-w w'/(r Omega).            (1)

For positive toroidal n the m=0 transport interval is now nonpositive.
The positive-frequency positive-energy bending branch lies above it.
The other poloidal transport intervals remain separated by Omega_min
at large R. Moving outside the band does not by itself give a uniform
threshold inverse; the following calculation supplies the required
response-specific estimate on a precisely stated forcing class.

## 2. Exact forced equations and a particular-displacement cancellation

Force the COMPLETE axisymmetric linear momentum equations by

    F=eta cross omega0,
    (r eta_r)'/r+i k eta_z=0.                      (2)

This is the coadjoint velocity-source class: the full Leray projection
is represented by the actual pressure. No radial/axial pressure term
is deleted. Put omega=k c, d=c+w and v_r=-i k d f. Direct elimination
gives

    v_theta=-Z f+i F_theta/(k d),
    v_z=P/d+w'f+i F_z/(k d),
    P'=(k²d²-Kappa²)f+F_r+2i Omega F_theta/(k d),
    f'+f/r=P/d²+i F_z/(k d²).                     (3)

For (2), F_r=Z eta_theta-w' eta_z,
F_theta=-Z eta_r, F_z=w' eta_r. Make the EXACT change

    f=Y/r-i eta_r/(k d),
    P=d²Y'/r-d eta_z.                              (4)

Using only the divergence condition in (2), every apparent w'/(kd)
singularity cancels, leaving

    [(d²Y'/r)]'+(Kappa²-k²d²)Y/r
                   =Z eta_theta+d(eta_z'-i k eta_r).              (5)

The reconstructed ACTUAL velocity is

    v_r=-i k d Y/r-eta_r,
    v_theta=-Z Y/r,
    v_z=(dY)'/r-eta_z.                             (6)

The particular terms in (4),(6) are retained. An estimate for Y alone
is not a license to erase them from a measured velocity or spin.

## 3. Uniform weighted variational bound

For real c>=0, the negative of (5)'s left side has quadratic form

    Q_c[Y]=integral[(c+w)²(Y'^2+k²Y²)/r-Kappa²Y²/r]
      =Q_0[Y]+integral(c²+2cw)(Y'^2+k²Y²)/r,
    Q_0[Y]=integral w²[(Y'-Y/r)²/r+Y²/r³+k²Y²/r].                 (7)

This follows from the actual force-free identity (1) and integration
by parts, not an asserted centrifugal stability criterion. At k>0
the full exterior tail is used. For c=0 use the completion of regular
compact test functions in the displayed weighted seminorm, identifying
its zero-velocity exterior-only representatives. On that physical active
quotient it is a norm. For c>0 the c² term restores coercivity on the
full global Hardy space. This also makes precise the c=0 positivity
shorthand in frozen 0201: its negative-c exclusion and c>W0 principal
pole use the strictly positive added terms and are unaffected.

Let L_eta(Y) be the pairing of the right side of (5) with Y. By (1),
the first term is bounded using the weight wY/r^(3/2):

    |integral Z eta_theta Y|
       <= ||sqrt(r) w' eta_theta/Omega||_2 sqrt(Q_c[Y]).

Hardy and (7) also give

    integral (c+w)²Y²/r³ <= 2 Q_c[Y].

Consequently

    |L_eta(Y)| <= C_eta sqrt(Q_c[Y]),
    C_eta=||sqrt(r)w' eta_theta/Omega||_2
               +sqrt(2)||r^(3/2)(eta_z'-ik eta_r)||_2.             (8)

These are actual smooth compact-generator norms, independent of the
distance c from the transport edge. Lax--Milgram in the weighted
completion gives sqrt(Q_c[Y])<=C_eta. The same argument applies in a
fixed right complex sector |Im c|<=theta Re c, theta<1: the real part
of (7) retains a fixed fraction of its positive terms. This is suitable
for an eventual positive-frequency contour, not merely a real test.

In particular the axisymmetric positive-c homogeneous problem is absent
for this sign of W. Its zero-c form is positive, so a positive bending
pole would not be crossing an unnoticed axisymmetric positive mode.

The estimate is for the weighted response and its weak force pairings.
It does NOT claim an unweighted L² resolvent bounded on arbitrary
forcing as c tends to zero.

## 4. Reconstructing physical velocity with the actual flat tail

The radial and azimuthal rows of (6) are bounded in L²(r dr) by (7),
the smooth profile and eta. The axial row additionally needs

    integral w'^2 Y²/r <= C_w Q_c[Y].              (9)

This is not true by comparing w'/w pointwise at a flat boundary.
For the actual 0195 tail w=sqrt(A)exp(-1/phi(r)), phi(a_c)=0,
phi'(a_c)<0, one has w w''/w'^2 tending to 1. Choose a fixed outer
annulus where w w''>=w'^2/2. Integration of

    (w w' Y²/r)'=
       [(w'^2+w w'')/r-w w'/r²]Y²+2w w'Y Y'/r

then bounds its w'^2 Y²/r integral by a fixed inner trace and
integral w²Y'^2/r. Both are controlled by (7); the trace lies where
w has a positive lower bound. Density extends the bound to the
weighted completion. This proves (9) with a finite profile-dependent
constant, and hence a uniform L² bound for ALL reconstructed velocity
rows in (6), including their particular terms.

The weighted form (7)-(8) needs only (1) and smoothness. The stronger
unweighted reconstruction statement explicitly uses this actual flat-tail
Hardy property. A different source profile can supply the same property
by its own tail proof; it is not inferred from the word smooth.

## 5. Exact interface to the remaining bending/current construction

An approximate isovortical velocity v=V xi has actual full Euler residual
V(xi_t-A xi), where V xi=P(xi cross omega0) and A is the Kelvin generator.
Thus the source class (2) is physically motivated, rather than chosen
from an arbitrary forced ODE. To use (8) for the ring, the COMPLETE
curved bending residual must be expressed in that class with its true
norm, after the fast poloidal and exact Euclidean reductions.

There is also a useful improvement to the coarse background estimate:
the log(R)/R target-dipole term of the exact Green kernel is canceled
by the ACTUAL -U_R x translating-frame term. Its source-dipole term
vanishes by the fixed center border. The remaining first-order inner
background distortion is O(1/R); its next correction may contain
log(R)/R². This cancellation should be used before estimating the
bending source, not discovered by a numerical tolerance adjustment.

If that same-source form bound is established at O(1/R), its Schur
feedback is a bounded finite contribution rather than another
logarithmic line tension. The positive n>=2 logarithmic form in the
first companion can then dominate it on a full frequency contour.
That is a concrete next proof step. It still requires the actual
curved residual and the literal Euclidean material current; neither
the positive restricted action nor (8) alone is declared that result.
