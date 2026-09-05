# Actual initial tagged displacement dipole and an additional material control

## 1. The initial position moment is not the matched spin

For one actual tag, write r for its centered background material
position and xi for its actual relative displacement. Define

    G_xi=∫tag rho r×xi,
    D_xi=∫tag rho(xi⊗r+r⊗xi).

These are literal displacement moments. The linear initial
centroid-plus-ambient displacement filter obeys

    rho(X0−U0)=i K×G_xi/2−i D_xi K/2+O(|K|²).        (1)

The centroid phases and the entire ambient complement are retained,
exactly as in 0117's velocity formula. The initial material displacement
filter is defined with the actual background label positions; evolving
U from the hybrid momentum retains the subsequent transport currents.

Differentiating the defining material integral gives

    delta S=G_xi,t+2∫tag rho xi×u0,                    (2)

with centered variables and their actual centroid/frame corrections.
Thus matching delta S to an optical canonical momentum does not imply
G_xi=j theta. In a circular core with xi_r=f, xi_theta=i f and time
factor exp(−i sigma t), one has G=i Lambda and
delta S=(sigma+2Omega)Lambda at the reference instant. The cancellation
near sigma=−2Omega can be order one in G even when spin is small.
Exterior returns cannot change G on the tag at the initial instant.

The symmetric D_xi density of a covariantly prepared axial optical
vector vanishes at zero macro K after the COMPLETE SO(3) average:
there is no nonzero isotropic axial-vector-to-symmetric-tensor map.
This does not remove individual shape moments or their higher jets.
The remaining initial first-gradient position condition is the actual
axial G row, not a presumed equality with canonical momentum.

## 2. Its exact radial row follows from the Euler/Lin field

In the 0142/0147 pressure convention, Z=2O+rO' and

    r(v_theta+Z xi_r)=m P/s,
    xi_theta=i v_theta/s−rO' v_r/s²,  xi_r=i v_r/s.

Eliminating v_theta and v_r gives

    r xi_theta=i[m P/s²−2rO xi_r/s].                    (3)

For a real radial displacement/pressure convention this is an imaginary
row, whereas the initial spin row mP/s is real. Integrate (3) against
the SAME nonnegative helical material fraction, spectral packet and
axial controls of 0147. In its parity coordinates the rows are

    theta0=(0,c0),   S0=(−eta beta_packet/c0,0),
    G0=(0,Lambda),
    Lambda=rho mu eps pi ∫∫[same packet/axial weights]
                    chi b [mP/s²−2rO xi_r/s] r dr ds.    (4)

The known metric, centroid and toroidal-frame corrections are included
in the actual finite-field integral. They have the controlled 0147
error, not an exact straight-frame identification.

For the eta=1 physical phase normalization, let
j_packet=−beta_packet/(gamma0 c0²)>0. The desired initial position
condition is precisely

    Lambda=j_packet c0=−beta_packet/(gamma0 c0).         (5)

If a standing-pair convention is used, replace this target by that
pair's ACTUAL action/spin coefficient after averaging, not by silently
changing eta. Equation (5) concerns one packet before division by cell
volume; dividing both sides by that same volume preserves the identity.

## 3. A new in-tag control actually exists

For the actual leading circular displacement in 0147,
xi is proportional to e_+ r exp(−x/2)L_n^1(x), x=g r².
Direct angular integration of r×xi against the tagged quadrupole
shows that G has radial row exp(−x/2)L_n^1(x). This differs from the
pressure-spin row exp(−x/2)P_n, P_n=L_n−2L_n'.

At the selected n=8,J=7, the existing sixteen exponential rows are

    x^j P_n,
    −x^j P_n/2+(3x/2)[(x^j P_n)'−x^j P_n/2], 0≤j≤7,  (6)

after removing their common exp(−x/2). Their exact polynomial rank
is 16 in the degree≤16 space. Adding L_8^1 raises that rank to 17.
The preserved first calculation derives these ranks from the actual
Laguerre polynomial, rather than counting proposed bump controls.
Together with the four independent reference rows 1,x,x²,x³,
the resulting 21 analytic radial functions are independent on every
open annulus. As in 0147, sufficiently narrow smooth bumps at points
with a nonzero evaluation minor give an invertible actual moment map.

Choose one additional radial direction for b0 annihilating all twenty
old radial functionals and with nonzero G functional. The second axial
control b2 is unchanged. At the reference carrier, its normalized
zeroth axial moment vanishes, so the leading G condition is a new b0
row. The enlarged matrix is block triangular: the old forty equations
and one nonzero new row. This supplies a 41st ACTUAL material-marker
control, not an exterior velocity chosen to change a tagged integral.

The angle numerator at t=0 is independent of the quadrupole marker
coefficient by angular selection, while its denominator is the fixed
Q_* and the reference angular-speed moment is already constrained.
Thus the target in (5) remains bounded in the normalized small-delta
system. It is of order delta in the natural G scale; canceling the
uncontrolled order-one G requires bounded order-one marker changes,
not a perturbation proportional to an exterior transfer error.

The exact forty-row map and its Jacobian were controlled in 0147.
The additional exact G row is a smooth perturbation of its new
independent leading row, including pressure, spectral and geometry
corrections. The finite-dimensional IFT therefore supplies the enlarged
41-row solution for sufficiently small delta and subsequent sufficiently
accurate transfer. The bounded new direction preserves Q_* and the
three reference dephasing moments. The existing full-time spin error
estimate still applies because it used those moments and bounded
absolute marker norms, not a unique forty-bump solution.

Choose the common small marker scaling to retain a strict nonnegative
fraction margin. No physical density is signed and no carrier-dependent
tag is introduced. This constructs (5) while preserving the existing
finite-time spin/action matching and its declared two carrier jets.
It does not claim that G(t)=j(t)theta(t) for every time or every carrier;
the separate current in (2) remains, and additional such conditions
would be distinct moment constraints rather than automatic corollaries.

## 4. The hybrid current with its actual initial row

Let A_S(t)=∫_0^t delta S(s)ds. With (1) as the actual initial datum,
0117 gives, retaining shape and pressure/stress remainders,

    U(t)=X(t)+(U0−X0)−i K×A_S(t)/(2rho)
                 +i[delta I(t)−delta I(0)]K/(2rho)+O_T(K²).

For a scalar helicity and a literal spin law S=I_s theta_t+chi_s theta+e,

    A_S=I_s theta−I_s(0)theta0
           +∫_0^t(chi_s−I_s,t)theta ds+∫_0^t e ds.       (7)

Only when the actual initial G target and the retained spin coefficients
match does this reduce to U=X−(j/2rho)curl theta plus its explicitly
known connections. The zero-frequency integration row is no longer
silently assigned. A coordinate-only shift proportional to theta does
not change {X,theta}; the integral/current row in (7) must also enter
the full joint Poisson and action calculation.
