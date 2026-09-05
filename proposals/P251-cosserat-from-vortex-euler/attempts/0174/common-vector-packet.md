# Common-vector optical packet: preparation, pressure, and physical observation

This construction uses the finite-action packet of0147, including its
distributed material-sheet marker, complete spin, and all clock connections.
Its angle remains that registered material observation, not an absolute
vorticity director. The common macroscopic wavevector is K=k*kappa in a
fixed laboratory frame. The same K is used before averaging whole fields.

## 1. A physical preparation, not a change of the measured clock

Fix the0147 packet, its reference p, its finite time interval, and its
positive axial curvature margin

    b_* = a_packet^2 * [3 sqrt(2)/(n+1)] Omega^2 delta / p^2.

All geometry, observation denominators, and time/parameter connections are
included in this margin. Let n0 be the tangent at the packet center xc.
Write its compact, solenoidal initial displacement as curl A_p; all
transverse and axial cutoff errors are already below the chosen margin.
For a periodic cell of side P containing this support in its interior set

    p(K)=p-n0.K,
    A_K(x)=exp[-i K.(x-xc)] A_{p(K)}(x),
    xi_K=curl_K A_K.                                      (1)

The phase is defined only on the compact support and extends smoothly by
zero. Thus A_K and xi_K are legitimate periodic Bloch amplitudes. Exactly,

    exp[i K.x] xi_K=exp[i K.xc] curl A_{p(K)}.               (2)

Equations(1),(2) specify the positive internal-helicity component. A real
common-K preparation has both signed carrier bands. If E_p denotes the
two real phase columns and J is their real complex structure, use

    Pi_+=(I-iJ)/2, Pi_-=(I+iJ)/2,
    E_K=exp[-iK.(x-xc)]
        [E_(p-n0.K) Pi_+ + E_(p+n0.K) Pi_-].              (2a)

Then E_(-K)=conjugate(E_K) exactly. The potential version of(2a) is used
before applying curl_K. A single real column shifted by p-n0.K would
not satisfy this reality condition. The two signed reference carrier
bands are disjoint; the full action uses their actual pairing, not the
sum of independently eliminated scalar oscillators. Subsequent small
geometric/periodic cross forms are retained in the full matrix.

The measured disturbance in this cell is the actual packet with its axial
carrier preparation changed by -n0.K. In a translated cell it has the
additional physical phase exp[i K.Pm]. The construction describes coherent
cell-center modulation; it does not claim that a macroscopically uniform
phase gradient is present at every point inside the material tag.

Choose the actual Kelvin velocity P_K(xi_K cross omega), then evolve the
full Euler/Lin equations. On full space, multiplication by exp[i K.x]
conjugates the Bloch equations to the ordinary equations exactly. Equation
(2) consequently persists as an identity between the actual full-space
histories. In particular no moving compensator is frozen into an Euler
ansatz: the physical solution is evolved, and its Bloch representation is
formed afterwards. In a curved tube the phase canceled in(1) is the full
Cartesian K.(x-xc), not K projected onto a changing tangent. Hence its
transport includes u.grad[K.(x-xc)] exactly. The difference between n0
and the local core frames is the already retained finite-arc geometry
error of0147. There is no omitted K times curvature commutator.

For this full-space preparation the complete phase action, material sheet
angle, tagged spin, and ambient currents are exactly the two0147 signed
families composed as(2a), with their actual cross pairing. The pure signed
reference bands have zero cross form by axial/angular orthogonality. The
second derivative of either signed parameter is n0 tensor n0. Odd spatial
terms and their physical clock connection remain until the declared whole
time-reversal pairing is applied. Neither a logarithm branch nor a physical-
frame winding has been chosen to make an energy positive.

## 2. The small-wave pressure issue and its actual repair

One cannot differentiate an absolutely convergent image series twice and
then assume it is still absolutely convergent. The Euler pressure return
has algebraic tails. We instead retain its Fourier symbol, including the
zero harmonic, and use the additional derivative in the Euler stress.

First the initial Kelvin force needs zero integral. On the exact0147
axisymmetric torus the packet contains only angular harmonics |m|>>1.
The background has harmonic0 in cylindrical components; conversion of a
vector to Cartesian components shifts this by at most1. Consequently

    integral xi_p cross omega = 0                         (3)

for every p in the preparation interval and both real phases. This is
exact for the sampled high-harmonic packet, not exponential smallness.
The Kelvin velocity P(xi cross omega) then decays as |x|^-4, whereas a
generic nonzero force integral would give |x|^-3. Its first two carrier
derivatives have the same property. Smooth cutoff and final approximation
errors are handled by the correction below, rather than called zero.

Fix the finite toroidal geometry first. Take a compact nonnegative cutoff
chi covering an open set where its actual velocity has all three
independent translation derivatives, and define

    G_ij=integral chi (partial_i u).(partial_j u),
    A^(j)=chi partial_j u,
    xi^(j)=curl A^(j).                                    (4)

The analytic genuinely three-dimensional CK torus has G>0. Indeed a zero
vector c in this Gram would give (c.grad)u=0 on an open set, hence globally
by analyticity. A nonzero translation invariance is incompatible with the
isolated closed circular core and its nondegenerate transverse return.
One can equivalently choose finitely many small supports witnessing this
linear independence. Their finite inverse is fixed before approximation.

Integration by parts gives

    integral (curl A) cross omega
       =lambda integral (Du)^T A                         (5)

on these constant-lambda backgrounds. Thus subtracting

    sum_j xi^(j) [lambda G]^-1_(j i)
                     integral (xi_p cross omega)_i         (6)

restores(3) exactly, separately for both phase columns and smoothly in p.
All controls and their actual Euler histories enter the full action and
observables. The cutoff/approximation force defect can be made arbitrarily
small before solving(6); its fixed inverse, Sobolev norms, spatial moments,
and actual finite-time Euler growth are included in the error budget.
No low-frequency return is propagated through an R->infinity estimate:
R and the controls are fixed first, approximation accuracy second.

## 3. The localized pressure-jet estimate

Here is the estimate used for the selected final periodic field. Its global
C^s bounds are fixed uniformly in approximation accuracy as in0153. Norms
of localized perturbations in this section are unnormalized full-volume
norms, not cell-volume averages. Pick s sufficiently large for the physical
observations, and a weight exponent 3/2<sigma<5/2.

For smooth compact Kelvin data satisfying(3), actual full-space linear
Euler histories on a fixed bounded-derivative background have finite

    ||v||_H^s + ||<x>^sigma v||_H^(s-2),                  (7)

on every fixed finite time interval, uniformly over the selected backgrounds.
The Lin displacement has the analogous local and weighted bounds. For
clarity, the cancellation needed in this assertion is dynamic, not an
assumption that all force moments remain zero: write the pressure as

    Delta pi=-partial_i partial_j (u_i v_j+v_i u_j).

Its velocity contribution has the order-one symbol

    B(q)T=-i P(q) q_j T_(.j),                            (8)

up to the retained transport convention. At spatial infinity its kernel
is O(|x|^-4). Split that kernel outside a unit ball; its tail applied to
the L1 stress has weighted L2 norm for sigma<5/2. The near-diagonal part
is the usual order-one local operator, controlled together with transport
by the finite Sobolev energy estimate. A dyadic partition and integration
by parts in divergence-free transport give the same estimate for the
weighted derivatives. The condition sigma>3/2 gives L1 control of v by
Cauchy--Schwarz. This closes the stress estimate using bounded u, without
an inverse-period Sobolev embedding. It also explains exactly why(3) is
needed for the initial order-zero Kelvin projection.

The apparent danger in the second K jet is explicit in Fourier space:

    |partial_K^j B(q)| <= C_j |q|^(1-j), j=0,1,2.         (9)

For j=2 this is only |q|^-1, which is square integrable at the origin in
three dimensions when multiplied by the bounded Fourier transform of an
L1 stress. The initial projection has the same bound: (3) gives
Fhat(q)=O(|q|), so two derivatives of P(q)Fhat(q) are O(|q|^-1).
The constants use the fixed first moments of the compact force.

More explicitly, for an infrared radius h and the lattice q=2pi m/P,

    integral_(0<|q|<h) |q|^-2 dq <= C h,
    P^-3 sum_(0<|q|<h) |q|^-2 <= C(h+P^-1).              (10)

The second inequality follows by counting lattice points in integer
shells: O(j^2) points times j^-2 per shell. This is the precise low-q
estimate, not period^2 times a small energy. The zero harmonic is separate:
P(k*kappa)=P_kappa there. Its finite-dimensional mean evolution is retained
exactly. For the compensated initial data its zero value vanishes by(3),
and its first/second directional jets have the cell-volume normalization
P^-3/2 in full-volume L2. The stress-generated harmonic likewise carries
this normalization. No direction-independent projector value at K=0 is
invented.

Only the undifferentiated stress needs the L1 bound in this argument.
In the twice-differentiated equation B'' acts on T0, B' acts on T1,
and B acts on T2 within the ordinary Euler propagator. B' is order zero,
so T1 requires its L2/Sobolev bound, not an incorrect L1 estimate for
x*v. Correspondingly the physical phase derivatives have the weighted
hierarchy sigma_j<5/2-j. Choose sigma_0>3/2, sigma_1>1/2, and any
0<sigma_2<1/2. The same dyadic estimate gives uniform L2 tails for
all three columns, with the slowest bound O(D^-sigma_2). A claim that
the twice-differentiated velocity remains L1 would be false; its generic
tail is |x|^-2. No such claim is used in(11).

Away from |q|<h the kernel and two derivatives are smooth. On a fixed
window of diameter D their periodic/full-space Green differences converge
with an estimate C_(D,h)/P. Equivalently, rescale the zero-mode-subtracted
cell Green function: its regular order-zero Leray kernel has local jth
phase derivative bounded by C P^(j-3), for j<=2. The free-space singular
part is canceled by the physical phase in(1). Smooth high-frequency
cutoffs and the fixed Sobolev bounds justify the kernel differentiation.

Combine this local estimate, the weighted tails from(7), and(10), and
then apply the actual Euler/Lin Duhamel energy estimate to the two phase
columns and their two directional K derivatives. The resulting local
observation/action estimate has the explicit ordered form

    error_(j<=2,T)
      <= C_T [tail(D)+sqrt(h)+C_(D,h)/sqrt(P)],            (11)

with tail(D)->0. Constants include the fixed data derivatives and marker
condition numbers. The harmless weakened power P^-1/2 includes the
infrared remainder in(10). This proof estimates the Fourier multiplier
and the localized stress, not an absolutely summed differentiated image
series. The order-zero propagation is bounded by the fixed actual Euler
energy constants; no minimum nonzero wave number is used there.

The cell's higher Taylor remainder and analytic neighborhood can depend
on P. Equations(9)--(11) assert uniform convergence of the first two
observed directional jets only. The exact harmonic rows may be direction
dependent before averaging. Their actual rows are retained, not assigned
zero or treated as an analytic full-vector projector at the origin.

## 4. Natural-scale transfer and whole-field averaging

The ordering is now explicit. Fix the0147 physical packet and its strict
curvature/action margins. Choose its cutoff and the correction(6) so
their complete finite-time C2 errors are below a prescribed fraction of
b_*. Fix the toroidal geometry and all resulting finite constants. Choose
the local approximation error below the same margin, then D, h, and the
periodic quadrature cell in the order specified by(11). Period can grow
without multiplying a bare Poincare inverse into the error. Both acoustic
0170 and optical requirements can be satisfied by the same sufficiently
fine quadrature after all these finite constants have been fixed.

For each realization, compute the exact full-volume KKS of the corrected
columns, evolve actual Euler/Lin, and use the literal material angle/spin
rows. Their first two common-K derivatives differ by less than the chosen
margin from the real signed axial family(2a). This includes mass,
time and parameter connections, the full tag centroid, and exterior return.
The axial positivity is therefore not inferred from a leading optical gap.
The permitted error is proportional to Omega^2 delta/p^2 itself.

Rotate the entire field, packet, controls and tag, keeping the same lab K.
Average these actual two-column phase actions before elimination, exactly
as specified for the common acoustic phases. Haar averaging the leading
axial tensor uses

    E[n_i n_j n_a n_b]
       =(delta_ij delta_ab+delta_ia delta_jb
                                      +delta_ib delta_ja)/15. (12)

The bounded directional-jet error in(11) averages with the same bound.
Isotropy and the actual parity/time-reversal pairing give the ordinary
isotropic second-order tensor after averaging; raw per-cell harmonic
projectors are not falsely declared polynomial. In particular, if the
zeroth kinetic mark normalization is E[n tensor n]=I/3, the leading
longitudinal/transverse directional weights are3/5 and1/5 respectively.
These weights now follow from one common-K physical preparation. They
are not obtained by independently rotating K in each fiber.

The weights apply componentwise to the transferred RAW action tensors.
They are not permission to average frequency squared instead of action.
For example even in a scalar autonomous comparison with H(p)=M(p)g(p),
the retained gradient combination is

    H''(p)-g(p) M''(p)=M(p)g''(p)+2M'(p)g'(p),           (13)

not just M*g''. Time-dependent connections and noncommuting phase blocks
add their actual terms. Thus this construction supplies the requested
common-vector license and bounds transverse errors at the positive0147
clock-curvature scale, but does not by itself prove that the final averaged
momentum-Schur gradient coefficient is positive. That coefficient is
computed from these full rows by0172. In particular neither the fourth
Haar moment nor the negligible transverse error erases(13).

The final object is a controlled fixed-time common-K optical action and
actual tagged response on the same selected stationary EPS-compatible
field. It does not supply an autonomous band, an absolute core director,
an averaged momentum-Schur sign by scalar frequency averaging, or an
acoustic-time approximation. The joint physical acoustic/optical
phase and current block remains the0172 consumer, not a consequence of
the separate positive diagonal blocks.
