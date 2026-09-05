# Same-field physical optical window: action controls and their actual costs

This continues `band-center-current.md`. It concerns its actual homogeneous
tag supplier plus actual off-tag Euler controls on the SAME fixed C016
field. It does not claim that the tagged vortex geometry changes, nor
that its acoustic/optical cross forms have already become the parent's
Cosserat pencil. No source from root0218 is used to upgrade that scope.

## Initial forms, real law and normalization

For every fixed smooth band width h>0 and fixed small body K, all source
configurations and velocities in the companion proof are actual smooth
Euler/Lin data. Initial phase and full energy are the integrals(10), not
independent oscillator inputs. Use the real two-sideband preparation and
equal whole-field inversion pairing before eliminating the initial phase.
The even two-phase skew form is a scalar beta_h(K)J. Its full conserved
energy is a real symmetric two-by-two form H_h(K). Both have well-defined
second jets computed from their defining full-field integrals.

The observed angle and the corrected physical spin/G fix a target
J_eff(K)>0 near K=0. Normalize the restricted even phase to the actual
observed angle/rate Wronskian times J_eff. Normalize the separately
computed complete energy to the corresponding oscillator quadratic form.
This is implemented by actual0210 quadrature phase returns and0205
phase-null energy returns, not by changing a coefficient in(10).

The scalar phase return has either sign, exact zero leading energy, and
fixed positive actual band norm. For an arbitrary finite symmetric energy
difference, its spectral decomposition and the positive/negative actual
energy returns span that difference. To avoid a singular square-root
parameter when a target eigenvalue crosses zero, add and subtract a
common sufficiently large positive scalar matrix before factorization.
All resulting actual amplitudes are smooth in K near0.

Collect initial cross-phase and cross-energy rows through order2 against
all fixed source and acoustic columns before selecting the off-tag band
profiles. The homogeneous finite-kernel construction in0205/0210 provides
nonzero profiles on disjoint fixed regular bands satisfying these rows,
with a positive norm lower bound. Coefficients are selected once for the
entire jet, not separately as an uncontrolled K-dependent kernel vector.
This makes the added controls leave the ORIGINAL source/acoustic cross
forms unchanged. It does not falsely assert those original forms zero.

Whole-field inversion removes the control's genuine first-K transport
energy term; choose its phase-normalization amplitude from the EVEN
target jet. An odd target amplitude multiplied by that odd transport
term would create an additional second-order energy row. Such rows are
included in the source even energy before normalization, not omitted.
The even second-order control phase/energy errors are respectively
O(|beta|/N²), O(|beta|/N), as proved in0210. The actual0210/0205 initial
data and complete nonlocal Euler propagator, rather than a local pressure
replacement, are used for all finite-K controls.

## Explicit compatible own-scale diagonal

Here the tag, nonzero Q, annular center, literal j0>0, microcell and finite
time window are fixed BEFORE h,N,K. Their constants may be large, but
they are independent of this subsequent diagonal.

The moment-flat main and dipole sources have, for0<=r<=3,

    ||partial_K^r xi(0)||_L2 <= C h^(-r-3/2).

Their material derivative has the same bound: T differentiates their
angle but not their energy profile. Smooth coefficients of Du and the
pressure Hessian are fixed. Consequently the initial phase/energy forms
are O(h^-3), and their third-K derivative is bounded by C h^-6. Their
exact homogeneous Lin histories obey analogous fixed-time bounds by the
actual flow formula, including its axial shear. Thus no unspecified
source-dependent macro remainder is being hidden.

For normalized quadrature returns of0210, leading target size O(h^-3)
requires amplitudes O(h^-3/2) times its normalized O(sqrt(N)) velocity
and O(N^-1/2) configuration. Taking three fixed K derivatives adds at
most the corresponding h^-r amplitude costs. Linear Euler's L2/H1
finite-time estimates and Lin transport give the conservative full-form
cubic bound

    |R_form| <= C_T h^-6 N³ |K|³.                   (1)

Indeed the return velocity H1 norm is at most C h^-3/2 N^3/2; the
full Jacobi form is quadratic in a velocity L2/configuration H1 pair.
Radial K derivatives of P_K on each nonzero fixed-cell Fourier mode
are bounded operators; its transverse harmonic bundle is retained.
This estimate overbounds the lower-order pseudodifferential gains used
in0210. Finite energy returns obey the same upper bound. The finite
cross-kernel coefficients have unit norm on a fixed finite profile
basis, so their selection does not add unbounded K derivatives.

Use the explicit sequence

    N_h=ceil(h^-4), |K_h|=h^22,
    moment order m=45, off-flow integration order q=6. (2)

The signed moment-flat profile is fixed at that finite order before
h shrinks. The resulting bounds at the physical second-jet scale are

    band observation error / |K_h|² = O(h²),
    cubic full-form remainder / |K_h|² = O(h⁴),
    even control energy-second-jet error = O(h),
    first-K tagged return error / |K_h|² = O(h^(1/2)). (3)

The last row uses the actual off-flow bound
C h^-3/2 N^-q after0210's sqrt(N) normalization. Its zeroth tag error
is EXACTLY zero on separated invariant supports; it is not divided by
K_h² as a nonzero persistent error. Second-K return observation errors
are smaller. The return phase-second-jet error is O(h^5).

Moreover |K_h| times the largest prepared velocity H4 norm is bounded
by C h^22 h^-3/2 N_h^(9/2)=O(h^(5/2)), tending to zero. The source
H4 costs are smaller. Thus the physical phase and spin remain finite
and positive, while the actual growth of the finite preparations is
compatible with the displayed continuum scale. This is a fixed-time
limit; no acoustic-time estimate or uniform full-field norm is claimed.
In particular K_h times a norm tending to zero does not make an O(1)
optical angle a uniformly small nonlinear Euler perturbation. This is
the actual linear prepared-response/second-variation continuum. Every
selected finite-h datum is smooth and has a finite norm; a nonlinear
perturbation amplitude must be chosen against that norm separately.

The true moving physical chart keeps all finite-h connections. Since
the angle and its time derivatives converge at the rates in(3) and its
limiting Wronskian is nonzero, those connections converge to the stated
autonomous second-jet chart. No time-dependent action is silently called
autonomous at a fixed nonzero h.

## Generic laboratory K and the isotropic optical interface

The source's normal support is compact in a contractible elliptic island.
Replicate it in normal lattice cells with phase exp(iK_perp.R_cell).
Within its support the corresponding periodic Bloch coefficient contains
exp[-iK_perp.(x_perp-R_cell)]. This is a genuine smooth local normal
preparation and produces no pressure interaction for the w=0 source.
It leaves the wrapped axial factor and the exact shear in the companion
proof untouched. Off-tag Euler controls instead use full P_K and the
full-pressure estimates already included above.

The observed tag angle/spin are the Fourier amplitudes of actual cell tag
moments. They are not the pointwise Eulerian Fourier velocity. Complete-
fluid means retain the cell factor exp(-iK.x), with its actual first and
second moments; these remain in the cross/current interface. Normal
tag-density centroid and covariance rows follow from the actual transported
delta_chi, while the added tangent control changes none of them.

For the whole-field law take a single common lab K and common vector
optical input, with n the actual axial tag direction. In body coordinates
K_X=n.K and |K_perp|²=|K|²-(n.K)². The raw squared-clock correction is
2nu0[b_parallel(n.K)²+b_perp(|K|²-(n.K)²)]. The physical vector angle is
Phi=3E[n theta], while spin/current and action have no inserted factor3.
Therefore j_hom=j0/3 and the derived transverse/longitudinal squared-clock
curvatures are

    B_T=2nu0(b_parallel+4b_perp)/5>0,
    B_L=2nu0(3b_parallel+2b_perp)/5>0.               (4)

These follow from the actual common-input phase/action and the Haar
fourth moment. Each prepared physical first clock derivative is zero
by the implicit center, so no first-column frequency covariance has
been discarded. The physical gradient mass from the corrected S/G
rows is retained in both mass and stiffness; subtracting nu0² times
that mass leaves exactly(4)'s positive curvature numerator. This is an
optical prepared-continuum interface, not a frequency-average postulate.

The original acoustic microscopic cross forms, their physical mean/curl
map, and the full coupled phase-chart positivity still need their actual
join. Most importantly, the source inside the tag remains homogeneous
Lin with w=0: field changes occur in its real off-tag action controls.
The field-changing vortex supplier and full EPS mechanism are not proved
by this material optical result. They continue on their original scope.
