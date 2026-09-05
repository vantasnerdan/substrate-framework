# Axisymmetric two-polarization hybrid response on the fixed ring

This is 0250's physical gain/history use of the root-authored axisymmetric
extension in `0260/vector-lift.md`.  The exact cohomological lift, its
mean-zero proof and its symbolic check belong to 0260.  Here they are inputs
to the literal material-momentum observation.  A tagged centroid is not used
as a substitute for the full point-to-hybrid coefficient.

## 1. Same-frequency centroid/current triangular pair

Use the circular thin-ring frame

    r=R+s cos(theta),       z=-s sin(theta),
    u=V(s)e_theta+W(s)e_varphi,       nu=m V/s,            (1)

and the axisymmetric harmonics `n=0`, `m=+1` or `m=-1`.  Put

    f=A(s)exp(iNI+im theta),       C=N/(lambda nu).         (2)

The imported exact preparations associated with

    b_u=u,                  b_p=m partial_varphi=m r e_varphi

have, uniformly at every fixed derivative order,

    xi_u^I=-f/lambda,       xi_u=C u f+O(1),
    xi_p^I=0,               xi_p=C m r e_varphi f+O(1).    (3)

The `O(1)` terms are finite-carrier amplitudes, while the displayed tangent
terms are `O(N)`.  The whole-space Leray field and every pressure correction
are already retained by the exact lift and the 0112 continuation.

For the stationary action tag, every scalar configuration moment depends
only on `xi^I`.  Hence the `p` column has exactly zero initial centroid and
shape variation.  The `u` column has the axial centroid coefficient

    delta X_z=(rho/(lambda M_tag))
                 integral J chi'(I) f z dI dtheta,

whose first poloidal coefficient is nonzero because
`integral z exp(im theta)dtheta=-i pi m s`.
Toroidal symmetry makes every transverse centroid and covariance-tilt row
zero for both columns.

Their leading axial angular-displacement currents are

    G_u,z=rho C integral chi f r W dx+O(c_N),
    G_p,z=rho C m integral chi f r^2 dx+O(c_N).            (4)

The angular coefficients

    integral r exp(im theta)dtheta=pi s,
    integral r^2 exp(im theta)dtheta=2pi R s              (5)

show that both gains are nonzero and

    gamma_N:=G_u,z/G_p,z
       =[W/(2mR)][1+O(1/N)+O(s/R)]                       (6)

in the thin-ring chart.  Use the *actual* smooth finite-carrier ratio in the
combination

    C_A=C_u-gamma_N C_p.                                  (7)

It has `G_z=0` exactly at each action frequency, retains the nonzero tagged
axial centroid of `C_u`, and has zero covariance angle.  This is the promised
triangular centroid/current row, not yet the hybrid conclusion.

## 2. Complete material second-momentum variation

Let `B_ijl=integral rho chi U_i r_j r_l` be 0241's *central material
momentum* tensor, with absolute material velocity.  At the centered
axisymmetric reference state, its first variation is

    delta B_ijl=integral rho chi [
        (D_t xi)_i r_j r_l
       +u_i((xi_j-delta X_j)r_l+r_j(xi_l-delta X_l))]dx.  (8)

This formula includes the moved material domain and the background velocity;
replacing `D_t xi` by the Eulerian perturbation `v` would omit the
load-bearing `O(N)` term.  It is the Lagrangian derivative of the defining
positive material integral, so pressure enters through the exact Lin
kinematics `D_t xi=Du xi+v`.

For the axial output, the reference contractions
`integral chi u_z x_j dx` vanish by axisymmetry and poloidal parity.  Thus
the centroid subtraction in (8), and therefore the moving-centroid phase in
0241's degree-two coefficient, vanish in this row.  In the circular leading
geometry

    u_z=-V cos(theta),
    (D_t xi_u)_z=C(V^2/s)sin(theta)f+O(1),                (9)

whereas the leading `p` displacement is purely azimuthal and

    (xi_p)_z=(D_t xi_p)_z=0.                              (10)

Substitution in (8), with the exact `phi` average, gives

    [delta B_zzz]_u
       =C V^2 s integral exp(im theta)
             sin(theta)[1-3cos(theta)^2]dtheta (2pi)+O(1)
       =i pi^2 m C V^2 s/2+O(1),                         (11)

and

    [delta B_zxx]_u=[delta B_zyy]_u
       =i pi^2 m C V^2[R^2/s+3s/4]+O(Rs)+O(1).           (12)

The common positive action-Jacobian, density, radial envelope and fixed-tag
dual coefficient have been divided out in (11)--(12).  Their signs are the
same.  The omitted `O(Rs)` terms are the smooth finite-thickness corrections
to the circular action chart, and the `O(1)` terms include the exact lower
carrier and pressure contributions.  In contrast,

    [delta B_zjj]_p=O(1),                                (13)

because both leading terms in (8) vanish.  Consequently the exact
`G`-cancelled column (7) retains the nonzero `O(N)` coefficients (11)--(12)
for one sufficiently large fixed `R` and all sufficiently large finite
carriers.  In particular, cancelling `G` does **not** cancel the full
material second-momentum row.

This calculation also exposes why the velocity-only expression previously
used as a diagnostic in `gain-matrix.md` was insufficient for this branch:
the actual material velocity variation is `D_t xi`, and its leading term is
one carrier order larger than `v`.

## 3. Literal point-to-hybrid coefficient and pressure ordering

For output along the ring axis and arbitrary macroscopic `K`, axisymmetry
reduces the degree-two term of 0241 to

    [Delta J_z]_(2)=-1/2{
       delta B_zxx |K_perp|^2+delta B_zzz(K.n)^2}.         (14)

The reference center is zero and the entire reference first-momentum row
`A_zj` vanishes, so (14) includes rather than drops the centroid-phase term.
The untagged point coefficient has no matching sparse tag harmonic and is
smaller than every preselected inverse carrier power by integration by
parts.  Choose the full-pressure order after the finite synthesis and tag
costs, exactly as in the fixed-tag diagonal already established in 0250.

On the selected transport band, `partial_t` multiplies the resonant material
coefficient by `-i omega`.  Equations (11)--(14) therefore give the actual
point-to-hybrid acceleration gain

    Gamma_A(n,K,omega)=
       h_perp(omega)|K-(K.n)n|^2+h_parallel(omega)(K.n)^2, (15)

where, after one common phase choice,

    h_perp h_parallel>0,
    |h_perp|+|h_parallel|>0                              (16)

uniformly on a smaller compact frequency band.  The ratio in the circular
limit is

    h_perp/h_parallel=2R^2/s^2+3/2+O(s/R)+O(1/N).         (17)

The initial check reported half of this ratio because it combined the
`cos(phi)^2=1/2` average in the transverse integrand with the full `phi`
integral used for the parallel row.  That script/output is preserved with the
`-initial` suffix.  The active check differentiates (8) in Cartesian
components and integrates both `theta` and `phi` with one normalization.  The
correction doubles only (12)/(17); their nonzero common sign and every
downstream rank statement are unchanged.

Thus the pressure, lower-carrier, first-moment and centroid corrections
cannot cancel (15) after `R` is fixed sufficiently large and then `N` is
chosen sufficiently large.  Neither parameter depends on the subsequent
finite-window target accuracy: higher WKB order and the sparse carrier
diagonal absorb that accuracy on the same ring.

The degree-one material dipole is odd in `K`.  Averaging the complete
reflected O(3) law eliminates it because no O(3)-invariant polar rank-three
tensor exists.  This is a whole-law parity statement, not a deletion of a
local tag moment.

## 4. Covariant acoustic vector reconstruction

Rotate the entire ring, tag, source, pressure and current as in 0257.  For a
ring axis `n`, define the exact normalized quadratic gain

    q_N(n,K,omega)=Gamma_A(n,K,omega)/|K|^2.              (18)

Equations (15)--(17) give `|q_N|>=c_*>0` uniformly in `n`, on the smaller
frequency band and for `K!=0`.  Given a desired polar acoustic spectral
amplitude `D_hat(omega)`, feed that realization the signed source amplitude

    a_n=3 (D_hat.n)/q_N(n,K,omega).                        (19)

The positive orientation and material laws are unchanged.  Only coherent
linear initial data are signed.  Since `E[n tensor n]=I/3`, the averaged
physical acceleration is exactly

    E[q_N a_n n]|K|^2=|K|^2 D_hat.                        (20)

The cost of (19) is bounded by `3/c_*` times the already finite reference
source cost, uniformly in orientation.  Applying the constructive 0254
finite-band formula first and the full-pressure/sparse-tag diagonal second
therefore gives both real parities and arbitrary prescribed finite-window
`C^r` acoustic correction histories on this one fixed ring.  No frequency
approach to zero and no radius growing with accuracy is used.

## 5. Reflected whole-law block and the surviving odd rows

Before whole-law averaging the acoustic column (7) is not optically silent
merely because `G=0`.  Its literal mechanical spin is

    S_A=2rho integral chi (xi_A cross u)dx,               (21)

because `G_A=0`.  The `u` leading tangent is parallel to `u`, but the
subtracted `p` tangent satisfies

    b_p cross u=-m r V e_s,                              (22)

so the local row (21) is generally nonzero.  Conversely, a local optical
column can have a quadratic polar hybrid moment.  These local coefficients
must be transformed with the complete reflected source map before deciding
the common gain.

Pair every orthogonal frame `R` with `-R`.  For an axial optical input `B`,

    B_loc(R)=det(R)R^T B,
    B_loc(-R)=B_loc(R),          K_loc(-R)=-K_loc(R).     (23)

A polar hybrid output is multiplied by `R`, hence by `-R` in the paired
realization.  Every even-`K` optical-to-hybrid coefficient therefore cancels
exactly.  For a polar acoustic input `A`, the local input changes sign under
`R -> -R`, while an axial output is multiplied by
`det(R)R`, which is unchanged.  Every even-`K` acoustic-to-angle/`G`/`S`
coefficient—including (21) at order zero or two—also cancels exactly.

This uses the full O(3) law.  The positive tag fractions `f_1,f_2` are scalar
and unchanged by reflection, so the optical determinant is not sign-averaged
away.  Likewise `q_N` in (18) is a scalar quadratic function of `K` and the
unoriented axis: `q_N(-n,K)=q_N(n,K)` and it is unchanged when the local `K`
reverses.  Formula (19) changes sign with the local polar input, exactly as
the rotated acoustic source must.  Its polar hybrid output then acquires the
second sign from `-R`, so (20) adds rather than cancels.

Odd cross rows survive and are not called zero.  With the 0241/0246
convention `C=i[K cross]`, their required leading physical observation is

    U_H=A-j C B/(2rho)+o(|K|^2),
    Phi=C A/2+B+o(|K|^2).                                (24)

The current representative also retains

    S_int=S_full-div Q,
    Q_ij=q(t)epsilon_ijk U_t,k,
    G_full-G_int=Delta G(0)+div integral_0^t Q(s)ds.      (25)

Thus the acoustic odd spin/current row and the optical `Q` contribution are
the literal curl and endpoint-memory terms in (25), with their actual signs
and time derivatives.  Reflection removes only forbidden even cross parity;
it does not delete (24)--(25).

## 6. One common frequency band and joint inverse

The axisymmetric acoustic frequency is

    omega_A(I)=m Omega_1(I),                              (26)

while the `m=1,n=-1` optical band uses

    omega_O(I)=Omega_1(I)-Omega_2(I).                    (27)

On the fixed thin ring, `Omega_1` is nonconstant on the selected action
interval and `sup|Omega_2|=O(1/R)`.  Choose `R` once so that this shift is
smaller than one quarter of the `Omega_1` range, then shrink away from the
endpoints.  The exact images of (26)--(27) contain a common nonempty open
frequency interval.  This is an overlap on one finite ring, not an `R` limit.

On that interval, at even spatial order the reflected four-output gain for

    (hybrid ; theta,G,S)

is block diagonal: (15) is a nonzero scalar acoustic block and equation (20)
of `two-polarization.md`, followed by 0257's whole-law reconstruction, is an
invertible three-row optical block.  The odd entries are `O(|K|)` and are
kept explicitly in (24)--(25).  Therefore the exact finite-`K`, full-pressure
gain remains invertible for all sufficiently small `|K|` by its finite
Schur/Neumann inverse.  Apply that inverse to the complete desired 0241/0246
observation vector, not to four independently idealized local columns.

The 0254 finite exponential construction now acts on this common band.  Its
finite coefficient and gain-inverse costs are selected before the WKB order,
sparse carrier and long-wave `K`; the existing common diagonal makes the
full-pressure, observation and odd-row remainders `o(|K|^2)`.  This closes the
joint acoustic/optical history range on the exact fixed 0211 ring.  The full
KKS/Jacobi normalizer is still a separate source calculation.

## Result

`route_verdict: actual fixed-ring joint acoustic-hybrid and full-vector
optical history inversion established with complete reflected parity and
surviving odd curl/current rows`

`evidence_scope: root-owned exact n=0 cohomological lift, exact Lagrangian
second-momentum variation, nonzero same-frequency G-cancelled hybrid gain,
uniform O(3) vector reconstruction, common acoustic/optical band and
finite-K reflected block inverse with the 0241/0246 odd rows retained`

The minimum physical-gain repair is `none`: even cross leakages vanish by the
permitted reflected law, the actual odd rows perturb an invertible common-band
block, and their desired coefficients are the explicit 0241/0246 map.  The
separate Jacobi/KKS normalization, final 0145/0147 periodic transfer and
geometry/density obligation remain open.
