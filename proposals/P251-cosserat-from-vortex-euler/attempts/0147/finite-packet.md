# Actual finite-action Lundquist packet, with physical spin normalization

This is the selected construction after the failure-derived corrections
in intermediate-candidate.md. The global/toroidal replacement is proved
in toroidal-transfer.md. No arbitrary infinite-fiber action measure remains.

## 1. Higher excitation retains the positive physical clock

Use0142's exact constant-lambda Lundquist field and actual Kelvin-prepared
linear Euler histories. Put k=-p, delta=sqrt(lambda/p),
ell=(2/(lambda p³))^(1/4), Omega=lambda U/2, m=2. The finite-order
quasimode construction works for each fixed radial excitation n. Write

    L_n(x)=L_n^1(x), f_n=r e^(-g r²/2)L_n(g r²), g=ell^-2,
    P_n=L_n-2L_n', cD=sqrt(2)Omega delta,
    sigma_n=-2Omega+(2n+2)cD+O(delta²).

For the same fixed physical radial marker chi proportional to r²,
the exact leading transported quadrupole integral is

    a_n(w)=(-1)^n (w-2)^(n-1)(w-n-2)/[(n+1)w^(n+3)],
    w=1+i cD t.                                               (1)

Derive (1) by differentiating the Laguerre Laplace integral
integral x e^-tx L_n^1(x)dx=(n+1)(t-1)^n/t^(n+2).
Its only zeros lie at real w=2,n+2, away from Re w=1. At w=1,

    a_n=1, (log a_n)'=-2n-2-1/(n+1),
    gamma_n=-2Omega-cD/(n+1)+O_T(delta²polylog),
    p² partial_p² gamma_n²=3sqrt(2)Omega²delta/(n+1)
                           +O_T(delta²polylog)>0.              (2)

Thus higher n supplies more independent physical moments without
losing positive physical curvature. Intrinsic sigma_n² curvature
remains negative; no frame winding changes that fact.

Choose J=7 and n=8 before checking the enlarged moment matrix. These
are finite construction parameters, not numerical basis cutoffs for
a claimed exact eigenvalue. The leading positive KKS follows from
integral x e^-x (L_n^1)² dx=n+1; all corrections are controlled as0142.

## 2. Exact finite-action preparation

Fix p_*>0 and a smooth spectral cutoff chi_band, one for
|q-p_*|<p_*/8 and zero for |q-p_*|>p_*/4. Hold it fixed when p varies.
Set

    a_L(q-p)=L exp[-L²(q-p)²/2]chi_band(q),
    Xi_p(x,s)=(2pi)^(-1/2)integral a_L(q-p)xi_N(q,x)e^(-iqs)dq. (3)

Use the exact isovortical initial velocity P(Xi_p cross omega0), and
the actual full Euler/Lin evolution after preparation. The fibers
commute with the full-space pressure projector. Thus (3) is smooth,
solenoidal, of finite full-space H^s action, and C² in p. It is not
a single frozen-frequency product. Fourier superposition retains
all longitudinal envelope/end-pressure effects exactly.

For the real and imaginary preparations, Plancherel gives the FULL
KKS

    beta_packet(p)=integral |a_L(q-p)|² beta_fiber(q)dq>0.       (4)

Here beta_fiber is the actual per-unit-length KKS in the complex
normalization of0142. For a constant profile the integral is exactly
sqrt(pi)L beta_fiber, derived from (3), not an assigned Lz. The
q-dependent integrand and both carrier derivatives are retained.
Actual Euler evolution preserves (4). Spectral-cutoff errors are
bounded by a polynomial times exp[-c(p_*L)²].

Fix the optical time |Omega t|<=T0. Choose a polynomial hierarchy,
for example p_*L=delta_*^-3; L is fixed during carrier differentiation.
Arbitrarily high finite-order0142 quasimodes allow their actual Euler
remainder to be made smaller than the physical error below, including
the finite polynomial p,L and marker-conditioning losses. The global
Lundquist coefficient bounds and Euler energy constant do not grow
with p or L. No exponential-in-L stability estimate is assumed.

## 3. Actual distributed helical material marker

Fix a small constant c>0 and ell_z=cL. Let chi_z be a nonnegative
Gaussian of width ell_z with a remote smooth cutoff. Use the actual
initial material fraction

    w0=mu chi(r)chi_z(s0)[1+eps b(r,s0)cos(2theta-p_*s0)],
    mu<=1/2, |eps b|<1.                                      (5)

The label s0 labels actual transported material sheets. On each sheet
take the physical central quadrupole in the transverse plane, register
its known initial phase p_*s0, and compare its angle with the same
sheet under the unperturbed flow. Average these physical angle
variations with the fixed nonnegative chi_z weights. This is an
actual registered collective marker angle. It is not the unregistered
global quadrupole and not an absolute core director.

The spin S is the actual total tagged angular momentum about the tag's
physical centroid, projected along the longitudinal axis. All position,
velocity, moving-tag and ambient exchange terms are retained. For
each fiber its axial row is exactly rho integral w0 mP/s, using
the full variable-O identity r[v_theta+(2O+rO')xi_r]=mP/s.
The helical factor in(5) makes this total spin coherent across the
packet, so it scales as L, just as(4). Neither physical density nor
the definition of total spin is signed or reweighted after evolution.

Construct b so that on every sheet

    Q(s0)=integral chi(r)b(r,s0)r³dr=Q_* !=0.                  (6)

The angle denominator is therefore an actual nonzero quadrupole.
The small matching value Q_* has normalized size delta, not delta/L.
Its background dephasing has size O(delta³)/Q_*=O(delta²), with
no growing axial-length amplification.

## 4. Physical finite-packet carrier curvature

For uncut Gaussian envelopes the longitudinal angle weight is
exp[-L²(q-p)²/2-ell_z²(q-p_*)²/2]. Completing the square gives

    q_bar=a p+(1-a)p_*, a=1/(1+c²),
    var_q=1/(L²+ell_z²),
    angle amplitude factor=(1+c²)^(-1/2)exp[-D(p-p_*)²/2],
    D=L²ell_z²/(L²+ell_z²).                                  (7)

Taylor's formula under this exact Gaussian integral, with two carrier
derivatives of the actual profiles, gives

    p_*² partial_p² gamma_packet²|_(p_*)
       =a² [3sqrt(2)/(n+1)]Omega²delta_*
        +O_T(Omega²[delta_*²polylog+delta_*/(p_*L)²]).           (8)

The positive a² is part of the measured carrier response and is not
silently replaced by1. The large real factor in(7) cancels from the
phase derivative but remains in the mass and parameter connection.
Remote Gaussian cutoffs and their differentiated tails can be chosen
below the natural margin after these finite scales are fixed.

## 5. Full-time matching: two axial controls and 20 radial controls

The simple radial-only tag fails: spin contains exp[-D(p-p_*)²/2]
and Pi contains its inverse, so their second carrier jets differ by
2D times the initial matched row. Forcing old radial controls to
cancel this makes their inverse grow as(pL)² and invalidates the
small-Q dephasing bound. Two bounded even axial controls repair this
actual mechanism.

Under the normalized weight chi_z(s)exp[-s²/(2L²)], the functions
{1,s²/L²} have zeroth/second moment determinant2v²>0, with
v=c²/(1+c²). Dualize them into phi_0,phi_2 whose normalized M0,M2
moments are (1,0),(0,1). On the compact marker support they and
their smooth extensions are bounded. Write

    b(r,s)=b_0(r)phi_0(s/L)+b_2(r)phi_2(s/L).                  (9)

Set Q(b_0),Q(b_2) equal to the expansion coefficients of the constant
function1 in this dual basis, times Q_*. This gives (6) pointwise,
not merely after axial integration.

Otherwise the two radial controls give different background quadrupole
phase rates on different sheets. Even an O(delta²) sheet-phase error
can acquire the factor(pL)² when the carrier changes the longitudinal
weighting. Impose for BOTH b_0,b_2 three additional exact moments

    integral chi b_j r³[O(r)-Omega]^h dr=0, h=1,2,3.           (9a)

Every sheet's reference quadrupole is then constant through three
time orders after removing its common Omega rotation. Its normalized
remainder is O(delta^12/Q_*)=O(delta^11), or O(delta^5) after two
scaled carrier derivatives in the chosen hierarchy. These are exact
background moments, not a frozen angular speed. Signed lobe orientation
makes the cancellations admissible without negative material density.

Match the ACTUAL total spin to eta times the canonical momentum
from the actual two-column action, with eta fixed in advance
(eta=1 or the coherent standing-pair eta=1/2). Impose the first
three carrier jets at p_* of each time Taylor coefficient through
order J AFTER removing the common intrinsic factor exp(i sigma t).
Both real phases are included, using parity to remove identically
zero rows. These are 3(J+1) physical moment equations, plus two
reference conditions. All coefficients are derived by differentiating
the actual pressure/material observations and the exact action.

The leading radial row family is

    1,x,x²,x³, e^-x/2 x^j P_n(x),
    D_r[e^-x/2 x^j P_n(x)], 0<=j<=J,
    D_r=-1/2+(3/2)x partial_x.                                (10)

Constant shifts of D_r from the j-dependent prefactors are invertible
row transformations. The normalized exact reference rows approach
1,x,x²,x³, since O-Omega has a nonzero quadratic leading term.
For n>J these2J+6 rows are independent. To prove
this, divide out the exponential. A relation is of the form
A(x)P_n+B(x)xP_n'=0, with deg B<=J<n. P_n has simple positive roots
and P_n(0)=(n+1)². Indeed, between consecutive simple positive roots
of L_n, the strictly decreasing logarithmic derivative L_n'/L_n
crosses1/2 exactly once, with one further crossing beyond its largest
root and none below its first. Thus gcd(P_n,xP_n')=1. P_n would have
to divide B, forcing B=0 and then A=0. A polynomial cannot equal a
nonzero exponential polynomial, separating all four reference rows.
This proves (10).

For n=8,J=7 choose20 radial bumps with a nonzero evaluation minor;
narrow smooth bumps preserve that determinant. The exact check
derives a nonzero Wronskian minor, rather than assuming generic rank.
At leading narrow-band order b_0 controls the zeroth carrier/time
rows, b_2 controls the second carrier rows divided by L², and the
first carrier rows scaled by p use D_r on b_0-b_2. This last term
comes from differentiating the profile-envelope correction as well
as the carrier phase; treating the packet as a product would miss it.
Add the J+1 choices that the D_r rows of b_2 vanish. The resulting
4J+12=40 equations, including(9a), have block-triangular leading
matrix with two copies of(10). The exact averaged sheet-angle row
is not globally linear in b: higher time derivatives contain the
inverse of each sheet's reference quadrupole. By(6),(9a), these
inverses are uniformly close to1/Q_*. The exact normalized moment
map AND its Jacobian differ from the leading affine map by
O(delta polylog+(pL)^-1). The ordinary finite-dimensional IFT about
the explicitly invertible leading solution supplies the exact root.
No nonlinear equation is replaced by its initial linearization.

Solve these EXACT finite equations by that IFT, using beta_packet and actual
initial Euler derivatives, not frozen mode frequencies. Set Q_* of
normalized order delta so the coefficients stay bounded independently
of L. Common small scaling enforces the nonnegative fraction(5).
Every radial/axial control remains fixed in physical space when p or
t varies. This is a derived physical marker, not a fitted inertia.

## 6. Natural-scale full-time estimate and positive action

After extracting the common intrinsic phase, all radial Doppler
differences on the marker are O(Omega delta polylog). Taylor's
integral remainder through order J therefore costs
O_T(delta^(J+1)polylog). The normalized inverse, pressure moments
and Gaussian weights have bounded absolute moments. Two carrier
derivatives may cost(pL)², and that cost is explicitly retained:

    max_(j<=2) |p^j partial_p^j(S-eta Pi)|/Pi_scale
      <= C_T delta^(J+1)(pL)² polylog + high-order Euler error.
                                                                  (11)

For J=7 and pL=delta^-3 this is O_T(delta²polylog), strictly below
the order-delta curvature signal in(8). The reference quadrupole
terms satisfy the same bound: (9a) removes the first three time
orders, leaving an O(delta^5) differentiated contribution in this
hierarchy, while the exact reference phase remains in the action.
This is the reason for the higher
excitation and J, rather than an arbitrary inflation of a tally.

The actual two Euler columns have constant positive beta_packet and
an actual angle row c(t,p),psi(t,p), with psi_t<0. Its exact action is

    L=M/2[(theta_dot-(c_t/c)theta)²-psi_t² theta²],
    M=-beta_packet/(psi_t c²)>0, Pi=M(theta_dot-(c_t/c)theta).   (12)

All time and parameter connections remain. Their full coefficient
bounds follow by differentiating the nonvanishing physical row, not
by replacing its phase with an autonomous oscillator. The second-
carrier curvature of the physical equation retains(8), including
the time-connection correction already derived in0142.

The arbitrary-order0142 Kelvin/Euler residual is now chosen below
(11), including the polynomial losses from pL=delta^-3, the finite
40-row inverse and the necessary observation derivatives. This is
possible at a fixed sufficiently high finite order before delta is
selected. Its energy constant is independent of p,L. Choose the
nonlinear disturbance amplitude last at fixed geometry and T.

The achieved object is a finite-action actual Euler preparation and
controlled material history, not an assigned fiber norm. Its total
tag spin is literal physical angular momentum; its angle is explicitly
the registered sheet-quadrupole collective observable. The remaining
global geometry/pressure transfer is the next file, not an assumption
inside the present construction.
