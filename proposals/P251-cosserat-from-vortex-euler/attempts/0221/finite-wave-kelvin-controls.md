# Full Bloch pressure orders and the physical-error window

This continues the actual zero-wave Kelvin normalization, retaining its
large transverse displacement. It does not transfer0210's axial-only
projection estimate to a different generator.

## Exact solenoidal Bloch generator and mean conditions

Let D_perp=grad_perp+iK_perp, H_perp,K=-D_perp^2. On a scalar s with
zero normal Fourier mean, use its inverse on the nonzero normal modes;
for |K|<1/4 these have a fixed positive denominator. For c=-1 or c=0 set

    xi_K=(c s, J D_perp s+iK_X c D_perp H_perp,K^-1 s).    (1)

Its complete Bloch divergence is exactly zero: the divergence of the
last term is -iK_X c s and cancels the axial divergence. At K=0 it is
precisely the positive or negative generator already constructed.
Its actual Kelvin velocity is

    v_K=P_K(xi_K cross omega).                           (2)

All pressure and vorticity changes in (2) are retained. In particular
an axial velocity alone is not claimed solenoidal for K_X nonzero.
Equation(1) is a preparation, not an assertion that the finite-K
two-scalar sector remains invariant under Euler evolution.

Impose P1s=0 on BOTH sectors' quadratures, as a finite set of linear
constraints in the zero-wave construction. All multipliers in(1)
preserve Fourier indices; hence xi_K also has zero first shell and
zero mean. Since the background vorticity has only first-shell modes,
xi_K cross omega has zero mean for every K. Thus the initial coadjoint
forms below encounter no singular mean projector. The subsequent full
Euler evolution retains the actual ray-wise P_kappa convention for its
physical mean; that mean is not declared zero at finite K.

## Actual scalar pulled-back forms and their pressure orders

For two generator columns, the exact forms are

    Omega_K,ij=rho <omega.(conj(xi_i,K) cross xi_j,K)>,
    H_K,ij=rho <conj(v_i,K).(v_j,K+curl_K v_j,K)>.         (3)

These are the full Beltrami coadjoint forms for lambda=-1 and agree with
the constrained material Jacobi forms. A conjugate Bloch pair restores
the real physical field. Pressure is in the FULL P_K of (2).

Here is the order reduction that matters. Regard(1)--(3) as scalar
pseudodifferential forms in s on the normal torus, with K a bounded
covariant-frequency parameter. In the cone of large normal covectors p,
the leading generator symbol is i(0,Jp). Its leading vector polarization
is one real vector times i. Consequently its wedge with its conjugate
vanishes. The apparent order-two principal symbol of Omega_K is ZERO,
leaving an operator of order at most one.

Similarly, the leading velocity symbol is the REAL Leray projection of
that same single vector polarization times i. The apparent order-three
helicity symbol is proportional to v_0.(p_full cross v_0), hence zero
for the defining vector polarization, including K in p_full. This is
an identity of the entire principal symbol, not cancellation only after
setting K=0. The full scalar H_K therefore has order at most two.
The surviving subprincipal terms retain background derivatives, the
axial scalar c s, and the exact divergence completion in(1).

All K dependence enters through covariant differentiation, the full
Leray symbol and the inverse in(1). Each K derivative lowers their
joint parameter-symbol order by one. Thus, on the fixed smooth torus,

    partial_K^j Omega_K has order at most 1-j,
    partial_K^j H_K has order at most 2-j,               (4)

for the finite j used here. Products with the smooth background and
their commutators preserve these orders. The finitely many low Fourier
rows are smoothing remainders; their coefficients on the fixed regular
band are superalgebraically small in N by angular integration by parts.
The zero mode was already dealt with explicitly above, not by extending
K K^T/|K|^2 smoothly through the origin.

The whitened and phase-normalized streams have L2 norm
O(sqrt(|b|/N)). Their phase gradient is nonstationary on the fixed band,
so each negative Sobolev order contributes its actual N power. Equations
(3)--(4) give for the paired control

    Omega_K=bJ+O(|b|K^2/N^2),
    H_K=O(|b|K^2/N),                                    (5)

after positive whole-field inversion averaging. Before that averaging,
the first-K forms are retained and may be nonzero. For the axial optical
input, the whole-field pair R and -R has the same axial input direction
det(R)R e_X and opposite body K. It cancels the odd form coefficients,
not the physical polar current or its allowed first-gradient response.
The action is averaged without multiplying by the reconstruction factor
three; raw b contributes b/3 to the isotropic vector phase as in0214.

The two signs in(5) do not erase the individual energies of size O(N).
Their K=0 cancellation is the exact full-form normalization. Their small
second coefficients follow from the actual scalar principal-symbol
cancellations, including full pressure. In particular merely counting
||xi_perp||=O(sqrt(N)) without those cancellations would not prove(5).

## Full Euler/Lin remote observations, including order zero

The initial negative-sector velocity already has nonlocal pressure tails.
Its kernel away from the chosen source bands is smooth. The full linear
Euler/Lin propagator is transport along the actual background flow with
matrix/pressure pseudodifferential terms of order zero. On a preselected
finite time interval its kernel is smooth off that flow graph. Since the
observed tag and source bands have a fixed invariant psi separation,
the propagator restricted from the latter to the former is smoothing.

Repeated source-angle integration by parts therefore bounds every fixed
finite list of literal tag angle, G, S, centroid, symmetric shape, time
derivatives and first two K derivatives by C_q sqrt(|b|)N^-q for any
chosen q. The derivatives in the generator and the normalization cost
are included by enlarging q. Full-fluid mean/current observations use
their propagated smooth adjoint tests and obey the same bound, uniformly
in the ray direction of the mean projector. A canonical momentum has
not replaced any of these physical observations.

Unlike the axial-only0210 control, this includes a generally NONZERO
but superalgebraically small K=0 tagged error. The preparation norms
and cubic physical remainder constants grow only polynomially in N
for fixed sources, b and time interval. Enlarge a finite exponent L to
cover that complete list. Choosing

    K_N=N^(-L-1), q>2L+2

makes the zero-order tagged error divided by K_N^2 tend to zero, as
well as the first-/second-order observation errors and the cubic form
remainder divided by K_N^2. This is the stronger window required by the
new negative Kelvin sector;0210's exact order-zero invisibility is not
an input.

## Narrow-band optical source and a compatible joint diagonal

For a source band of width h, replace a positive averaging kernel by
a FIXED smooth signed moment-flat kernel of any required finite order r.
Its integral is one and its first r centered moments vanish. Signed
initial field amplitudes do not change the positive material tag. Taylor's
theorem for the actual frequency-coordinate observable kernels gives
O(h^(r+1)) physical output error through the fixed time/K derivative list.
The source Sobolev exponent in h is set by that derivative list, not by r;
r changes fixed constants of the chosen kernel. Its quadratic phase and
energy still have their actual finite inverse-width costs.

Let h^-P bound the required finite phase/energy corrections. First choose
a fixed M>P and N=h^-M. All preparation and third-remainder costs are then
bounded by h^-D for a finite D, enlarged to include each physical norm.
Choose r>2D+2 and an integration-by-parts order q with Mq-P>2D+2.
Finally set K=h^(D+1). Then the moment-flat source error and remote
control error are o(K^2), the actual cubic remainder divided by K^2 is
O(h), and (5)'s energy second coefficient is O(h^(M-P)). All vanish
while the physical tag and its measured j stay fixed. Constants may
depend on r,q and the already fixed source profiles; no uniform bound
over their design orders is asserted.

This joint ordering, not 'make K small last', is what licenses the
prepared diagonal. Actual optical clock/current second coefficients,
including any acoustic cross-observation, still come from their supplier.
The control changes full action forms with vanishing physical-output
error; it does not supply an absent optical or coupled Euler equation.
