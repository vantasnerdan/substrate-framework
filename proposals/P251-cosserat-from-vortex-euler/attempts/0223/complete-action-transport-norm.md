# Exact axial action removes the unnecessary spatial-Sobolev cost

This is the constructive norm repair of the last section of
`compact-complement-preparation.md`. It uses the ACTUAL scalar Kelvin
representation, full periodic pressure and coadjoint forms. It does not
prescribe an observed clock or assert the remaining physical-row rank.

## 1. Complete phase, with its real odd and even spatial terms

Use a left column at -k and a right column at +k. Write their scalar
data as (h_L,s_L),(h_R,s_R), and f_i=H^-1 h_i. The actual normal
generators are J grad s_L-ik grad f_L and J grad s_R+ik grad f_R.
The canonical KKS form is rho int omega dot (xi_L cross xi_R).
Integrating the normal Jacobian against the actual psi gives EXACTLY

    Omega(k)=Omega0+ik rho I1+k² rho I2,

    Omega0=rho int [s_L T s_R+h_L T s_R+s_L T h_R],
    I1=int [psi(s_L h_R+s_R h_L)
             -(s_L+h_L)grad psi.grad f_R
             -(s_R+h_R)grad psi.grad f_L],
    I2=int f_L T f_R.                               (1)

For h_i=-s_i the zeroth coefficient is -rho int s_L T s_R,
agreeing with the actual positive source0218. The imaginary first
coefficient is symmetric in the two inputs; the other two are
antisymmetric. Thus (1) has the required real Bloch phase parity.
Neither that first term nor its acoustic/optical counterpart is erased
by considering only an axial self block.

Every term in (1) is bounded by the scalar L² and T-graph norms.
The potentially large radial gradients in xi have canceled by exact
integration by parts, not by discarding physical kinetic energy.

## 2. Complete physical Jacobi energy

Let w_i=(w_Xi,v_i) be the FULL Kelvin velocities in (4) of the operator
source. Their normal curls are exactly

    c_R=curl_perp v_R=-T(s_R+h_R)
                           +ik div[psi grad(H^-1 h_R)],
    c_L=-T(s_L+h_L)-ik div[psi grad(H^-1 h_L)].       (2)

The complete Beltrami Jacobi matrix, at lambda=-1, is

    H_LR=rho int [w_L dot w_R+w_XL c_R+c_L w_XR
                                      +ik v_L dot J v_R]. (3)

This follows by integrating the two normal derivatives in
int w_L dot curl_k w_R once. It is the FULL energy, including helicity,
not just rho||w||² or the positive transport quotient. The physical
quadratic energy is one half the diagonal of this matrix.

All velocity components in (3) are bounded by
||h||+||s||+||Th||+||Ts|| through the full pressure formulas.
The same holds for (2), because div(psi grad H^-1) is order zero.
Thus (3), including the genuine odd transport term, has an actual
T-graph norm bound. At k=0 and zero normal harmonic it reduces to

    H_LR/rho=int [-Ts_L Ts_R-Ts_L Th_R-Th_L Ts_R
                  +T(s_L+h_L) H^-1 T(s_R+h_R)].     (4)

The passive h=-s and the h=0 competitors consequently have exactly
the opposite energy signs already derived in0218. This is an exposing
check on (3), not a newly assigned energy convention.

## 3. The actual finite-wave response has the same norm bound

On the fixed axial reflection sector the exact scalar generator is

    L(k)=-T I+A(k),

where A(k) is order zero, analytic for fixed small k, including its
explicit mean-pressure row. All fixed k derivatives are bounded order
zero. The commutator [T,A(k)] is order zero, as are any fixed number
of its iterated commutators: the background is smooth, and differentiating
its symbols in the transport direction does not raise their order.
The zero normal mode was treated separately in the exact formulas;
no singular 1/k factor is part of A(k).

Therefore the norm sum_(j<=r)||T^j q|| is propagated for each fixed
finite r with C_(r,T) bounds uniform for small k. Differentiating the
actual equation in k gives the same bound for any fixed jet through
its inhomogeneous Duhamel equations. In particular the cubic remainder
of the actual scalar response is controlled by these norms of the
selected initial coefficients, not a full isotropic H^4 norm.
The physical angle, G, spin, centroid and shape are bounded test rows
in these norms, by the exact formulas already derived. Their required
time derivatives use a finite higher T-graph norm.

The phase and physical energy remainders follow from (1)--(3) and the
same analytic pressure multipliers. For the axial optical action it is
therefore sufficient that the spectral construction have polynomial
costs in these scalar transport norms. This is precisely the cost
established for its compact forcing. It is not necessary to demand a
polynomial bound on a full spatial norm of grad z merely to control
this linear action or its physical observations.

## 4. Smooth fields, true amplitude scope and the remaining output row

Smooth functions form a simultaneous core for any fixed finite set of
these transport graph norms (periodic Friedrichs smoothing and the
transport commutator give the graph-core approximation). Approximate
each constructed scalar preparation in that norm to its allocated
error, retaining means/parity by the bounded actual symmetry projections.
The norm in the estimates above stays within the chosen finite budget,
even when the approximant has large uncontrolled transverse derivatives.
Its Kelvin generator is now a genuine smooth solenoidal field; the
actual pressure, material history, phase, energy and observations
converge by the proved bounded formulas. No distributional generator
is silently used as a finite microscopic realization.

A finite-amplitude nonlinear history requires choosing its amplitude
small enough for that selected smooth field's actual stronger norm.
This proof establishes linear prepared response and second variation;
it does not make that nonlinear amplitude uniform along the sequence.

This norm repair advances the registered compact-pressure continuation.
It does not remove its physical-output task: the corrected actual
theta/G/S constant, secular and second-order coefficients must satisfy
the common-clock/current constraints. Form normalization by actual
off-tag Kelvin controls is a separate input and cannot force those
tag coefficients to agree. The generic-K acoustic angle is also not
supplied by the axial slice. Both remain active parent obligations.
