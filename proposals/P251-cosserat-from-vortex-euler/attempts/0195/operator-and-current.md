# Full swirling-ring pole and the actual stationary material current

## 1. Complete fixed-harmonic operator

Use the actual global family and volume coordinates of swirling-ring.md.
At each fixed integer toroidal harmonic n, the compact active-vorticity
space includes every poloidal channel. Vorticity sources are supported
in the continued invariant solid torus, but their velocity is the full
R³ inverse curl B_n=curl(-Delta)^(-1), not a velocity truncated there.

With u=Omega partial_theta+V partial_varphi and
omega=alpha partial_theta+beta partial_varphi, direct Lie differentiation
gives the exact linearized Euler generator

    L_n w = -Omega partial_theta w-i n V w
             +w^I(Omega' partial_theta+V' partial_varphi)
             +alpha partial_theta(B_n w)+i n beta B_n w
             -(B_n w)^I(alpha' partial_theta+beta' partial_varphi). (1)

The new poloidal-vorticity derivative is retained. Dividing by a scalar
pressure mode or silently reusing the no-swirl compact K would lose it.
Both brackets preserve divergence and the actual active support. The
coordinate-vector formula is interpreted in a smooth disk chart at the
center, with its correct volume norm, not an unweighted singular dI basis.

Its Lie-transport part T has stretching matrix

    Nw=w^I(Omega' partial_theta+V' partial_varphi), N²=0.

In poloidal harmonic m its resolvent is exactly

    (z-T_m)^(-1)=d_m^(-1)Id+d_m^(-2)N,
    d_m=z+i(m Omega+n V).                            (2)

For fixed n and |Omega| bounded below, all full transport bands are
explicit; away from them the resolvent and partial_theta times that
resolvent are bounded. Compactness of the full inverse curl restricted
to a bounded source domain is the 0186 local H1/Rellich result. The
term alpha partial_theta B_n is order zero, but is relatively compact
against the transport resolvent. Commute partial_theta through B_n:
B_n partial_theta(z-T)^(-1) is compact, while the commutator kernel
gains the smooth coefficient difference and remains order minus one.
Thus the full pressure/stretching response yields the same meromorphic
Fredholm factorization off the shifted bands, not a discarded exterior.

## 2. Continuation at fixed finite R and n

Fix all three 0186 integer harmonics n-1,n,n+1 and their isolated positive
contours before choosing epsilon. Flatten the continued source boundary
by a smooth volume-coordinate map. At this fixed R, the streamfunction
and meridional velocity change by O(epsilon²), while V, alpha and the
new vorticity derivatives are O(epsilon). The physical inverse curl
remains the same global operator; its pulled-back coefficients converge.

The transport difference has its only unbounded derivative in theta,
controlled by the no-swirl Euler resolvent's graph norm. The other new
derivative alpha partial_theta B_n is bounded L2 to L2, because B_n
gains one derivative locally. Consequently on each fixed contour C_n

    sup_C_n ||(L_epsilon,n-L_0,n)(z-L_0,n)^(-1)||
                                          <= C_R,n |epsilon|. (3)

The constant includes the actual finite toroidal n and all smooth core
derivatives. Choose epsilon so the right side is less than one for ALL
three contours. The Neumann resolvent identity preserves their simple
Riesz projections and excludes a shift into an unexamined ambient band.
Their continuous positive KKS/phase form prevents Hamiltonian escape
from the imaginary axis. The resulting modes have actual global Euler
velocity and pressure, and positive translating-frame action with the
physical E-U_epsilon I_z subtraction.

This perturbation does not claim uniform control at arbitrarily large n.
At fixed R, the positive 0186 discrete squared-frequency curvature has a
finite nonzero margin. The three continued eigenvalues and the chosen
normalizations depend smoothly on epsilon; taking epsilon still smaller
preserves that actual discrete sign. The R² in the discrete second
difference belongs in its finite constant. It is not hidden by an O(epsilon)
localization estimate or asserted to define a continuous laboratory K.

On the active core the off-band Kelvin reconstruction follows exactly
as in 0186: solve (lambda+ad_u)xi=v and use [u,omega]=0 to obtain
w=[omega,xi]. The full Euler KKS is evaluated on this actual tangent,
including the new interior poloidal vorticity. Its sign is continued,
not replaced by an absolute value. Exterior material relabeling freedom
does not truncate the global pressure or alter the compact vorticity
phase form.

## 3. General advected tags and the repaired stationary annuli

For clarity, an arbitrary positive material fraction is transported by
the ACTUAL base flow. In straightened coordinates, if initially
chi_0(I,varphi), then

    chi(t,I,theta,varphi)=chi_0(I,varphi-V(I)t).       (4)

Unless V is constant on its support, the angular patch shears. A full
torus average chi(I) would be stationary but generally annihilates the
high-n optical moment: Euclidean centroid/spin/covariance weights have
only finitely many toroidal harmonics. It is not a replacement for the
finite-arc optical tag.

For ANY actually advected fraction, let X be its base mass centroid,
r=x-X, eta=xi-delta X the centered material displacement and u_c=u-X_t.
The actual displacement and spin variations are

    G=rho integral chi r cross eta,
    S=rho integral chi [eta cross u_c+r cross D_t eta].

Reynolds transport, rather than freezing chi, yields the exact identity

    S=G_t+2rho integral chi eta cross u_c.            (5)

This includes the moving tag, moving position and centroid terms. It
reduces to the stationary formula only when the base fraction is indeed
invariant. For generic (4), both the covariance-angle chart and all
moments are time dependent. No scalar Euler energy/clock identification
follows merely by changing to an angularly moving tag frame.

The inner-supported g of swirling-ring.md executes a different repair.
On the OPEN tag neighborhoods F=g=g'=0. Therefore their actual
chi(I,varphi_phys) is stationary, not constrained against the Euler flow.
Their base centroid is constant: integral chi u_i=integral div(chi x_i u)=0
for a compact stationary fraction. The covariance has the same nonzero
gap and the same unit rigid-rotation response as its continued 0186 tag.

## 4. One actual current/control map for all three modes

For each isolated mode choose its phase by one fixed scalar Euclidean
covariance-tilt component theta. Its complete actual material columns are
xi=xi_q theta+xi_p theta_t. Other measured components, including conjugate
tilts, are not deleted. Time reversal takes the actual stationary field
to -u, and registers the COMMON INITIAL scalar angle/rate pair with the
conjugate fluid phase reversed, just as in 0186.

For the fixed physical component define from the actual full fields

    G_q=rho integral chi r cross xi_q,
    C_p=rho integral chi xi_p cross u.

Applying (5) to both ACTUAL preparations gives

    G_average=G_q theta,
    S_average=(G_q+2C_p)theta_t.                     (6)

The individual G_p and C_q terms cancel only after the full prepared
TR average. The condition to recover literal current is the measured
C_p=0, not an assumption about a supplied rotor mass.

At epsilon=0, 0186 has one positive three-parameter annular tag with
C_p(n-1)=C_p(n)=C_p(n+1)=0 and invertible tag Jacobian. The mode, active
material reconstruction, centroid/covariance chart, fraction and actual
velocity all vary smoothly in epsilon. The same finite-dimensional IFT
therefore supplies a nearby positive tag fixing these three rows exactly.
Choose its parameter neighborhood inside the predeclared open region
where g is identically zero. Thus the control adjustment does not break
the stationarity that licensed (5), and the derivative is taken on the
same actual field, not an independently transplanted pressure profile.

The overlap G_q remains positive. Hence, on each of the three actual
linear modes and its actual registered TR law, for every time,

    G_average=J_tag theta,
    S_average=J_tag theta_t,
    G_average(0)+integral_0^t S_average=J_tag theta,
    J_tag=G_q>0.                                    (7)

The initial Euler phase forms and E-U_epsilon I_z energies of the two
prepared partners pull back to the SAME physical initial angle/rate
form. They are averaged with their positive weights BEFORE eliminating
the conjugate coordinate. This yields the positive mode action with its
actual mass M; J_tag/M is retained, not forced to one. Because the tag
is actually stationary there is no unaccounted moving-frame angular
momentum subtraction or tag-Doppler replacement in this action.

The result is a positive scalar optical supplier on a genuine swirling
closed ring with a nonzero periodic core and shared inner streamline/
vortex-line tori. The full ambient mean/current assembly and a continuous
common-K optical response are not inferred from three discrete modes.
