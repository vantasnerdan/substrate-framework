# Fixed-ring response: exact band criterion and the remaining Euler lift

This is the substantive calculation frozen by `README.md`.  It compares the
two registered routes without assuming a verdict.  The compact geometry and
stationary assembly owned by 0248 are inputs only; nothing below constructs or
reviews that assembly.

## 1. What a fixed transport band would actually prove

Let `J` be a nonempty open frequency interval and let
`G:J -> C^(d x d)` be smooth.  After restricting to a smaller open interval,
assume `G` has a bounded pointwise inverse.  For a compact time interval
`[-T,T]`, define

    (K a)(t)=integral_J exp(-i omega t) G(omega)a(omega)domega,
    a in C_c^infinity(J;C^d).                            (1)

For every finite `r`, the range of `K` is dense in
`C^r([-T,T];C^d)`.  This statement permits arbitrarily large, but finite,
source norms.

The proof is an exact annihilator argument.  If a continuous functional
`Lambda` on `C^r` kills the range, regard it as a compactly supported
distribution of order at most `r`.  Varying `a` in (1) gives

    G(omega)^* Fourier[Lambda](omega)=0 on J.             (2)

The Fourier--Laplace transform of a compactly supported distribution is
entire.  Invertibility of `G` and the identity theorem make it zero
everywhere, hence `Lambda=0`.  Hahn--Banach proves density.  Complex
conjugate preparations give real histories; sums and differences give the
cosine and sine parity sectors separately.  The same proof works for the
even or odd subspace after using its corresponding parity extension.

This result is stronger than persistence of every finite Taylor determinant:
it controls the whole finite window.  It is also deliberately weaker than a
uniform observability estimate.  For every requested target and tolerance a
smooth compactly supported `a` has finite norms, but those norms can diverge
arbitrarily fast down an approximation sequence.  That is allowed by 0250's
quantifiers so long as the microscopic field generated from `a` is an exact
solution on the one fixed background.

For several physical observations, `G` is the actual pointwise gain from
independent microscopic source polarizations to the independent target
channels.  Thus the useful fixed-ring theorem needs a full-row-rank gain for
the acoustic hybrid row and the optical angle/full-current rows, not merely a
nonzero scalar tag coefficient.  Phase and energy are quadratic rather than
linear output rows; after (1), their complete Gram correction must still be
solved in the observation kernel with the required positivity.

## 2. The exact positive band already present in the straight core

The source body 0222 supplies an actual field-changing, full-pressure Kelvin
lift in the straight literal-constant-curl region.  For a compact scalar `S`
independent of the axial coordinate,

    xi=(J grad_perp S,lambda S),
    w=P(xi cross omega)=-lambda (T_p S)e_z,
    S_t=-T_p S.                                          (3)

The discarded part of `xi cross u` is an actual compact gradient; `w` is
solenoidal, changes Eulerian vorticity, and solves Euler/Lin rather than a
passive-label equation.  For the `m=-1` pair, 0226's measured straight-core
history has the form

    h_F(t)=L[F(s) exp(-i Omega(s)t)].                     (4)

On any positive-radius subannulus on which the Bessel clock `Omega(s)` is
strictly monotone and the tag weight is nonzero, change variables
`omega=Omega(s)`.  Equation (4) is exactly (1) with a scalar nonzero gain.
The compact inverse-return row does not restrict this observed amplitude:
0222 places its correcting return outside the tag support, where it changes
the admissible Kelvin field but not `L`.
Section 1 therefore upgrades the previous finite moment statement at this
linear-observation level:

> On the fixed straight literal-curl supplier, compact smooth
> field-changing Kelvin preparations have dense measured angle histories in
> every finite-window `C^r` topology, with both real parities and finite
> preparation norm for each approximation.

This is not yet the 0250 object.  The straight column is not the compact
Euclidean ring, and the result does not make the full spin/current gain or
the acoustic hybrid gain pointwise invertible.  The separate exact quadratic
root in 0226 proves `H=nu beta>0` for its finite moment construction, but the
source does not prove a simultaneous arbitrary-history, complete mixed-Gram
normalizer.  Those distinctions prevent the new density lemma from silently
promoting more than its linear angle scope.

The new 0251 source receipt supplies a second exact positive input.  Baldi's
Theorem 1.1 gives an analytic action-angle chart on a toroidal shell of the
localized Gavrilov flow, with two particle frequencies and strictly varying
frequency ratio; the chart is independent of the localization cutoff.  At
least one integer frequency combination is nonconstant on an open action
interval, and hence supplies an exact compact-geometry passive transport band
to which section 1 applies.  This is a real fixed-geometry advance.  The
theorem concerns particle advection, however; it supplies neither (3) on the
curved field nor an intertwining map into the projected linear Euler
generator, nor the physical acoustic/optical gain matrix.

## 3. The finite-radius failure of the available lift

The exact straight identity (3) uses the constant axial translation field.
When 0222 maps it to the actual closed ring, the corresponding uncompleted
toroidal field has the exact finite-radius defects

    eta cross u=-grad(FS/r^2)+(TS/r)e_varphi
                              -(2FS/r^3)e_r,
    div eta=-lambda S_varphi/r^2.                        (5)

The Piola step repairs divergence of the initial displacement, and the
actual initial velocity and subsequent history use the whole-space Leray
projector.  It does not turn the transported comparison field into an exact
band solution.  The full-pressure estimate instead bounds its normalized
finite-window error by a radius-decaying quantity, with constants selected
after the finite preparation inventory.  The later 0226 recursion cancels
any specified finite collection of polynomial pressure tails and computes a
finite-order residual, but explicitly leaves the curved observation rows and
their remainder constants as order-dependent inputs.

Equation (5) is the load-bearing construction failure for the presently
available Route A: for generic nonzero `S` its radial term and divergence are
nonzero at every finite radius.  The existing proof earns

    for each fixed preparation inventory and accuracy, choose R large,

not the frozen order

    choose one R_* and then satisfy every preparation accuracy.           (6)

This is not a no-go for exact fixed-ring Euler response.  The actual Euler
evolution from each Kelvin initial datum exists, and its exact observed
family could still be complete.  What is absent is a proof that its full
pressure-coupled kernel retains the transport band's dense range.

The acoustic side has the same boundary more sharply.  C-CST-017's accepted
periodic construction uses both quadratures of an exact wrapped-streamline
passive sector to correct the complete point-to-hybrid acceleration.  Neither
0222 nor 0226 constructs the corresponding actual ambient-compensated
acoustic preparation on the compact ring; 0226 records that row as open.
Baldi's particle frequencies do not by themselves change this conclusion.

## 4. Why small curvature and finite ranks cannot repair the inference

The synthesis map in (1), viewed from a fixed `L^2(J)` source norm into
`C^r([-T,T])`, is compact: its kernel and one additional time derivative are
uniformly bounded on compact sets, so bounded source sets are relatively
compact by Arzela--Ascoli.  It cannot have a bounded right inverse onto the
infinite-dimensional `C^r` target, since the identity would then be compact.
Thus Route B cannot use an exact bounded Banach-space right inverse for this
natural smooth band synthesis.  It needs approximate inverses with an
explicit, accuracy-dependent source cost or a different noncompact source
topology.

Dense range is also not stable under a small operator-norm perturbation.  A
concrete model is

    T:e_n -> sigma_n e_n on l^2,  sigma_n>0, sigma_n ->0.                 (7)

`T` has dense range.  For any `delta>0`, choose `N` with `sigma_N<delta` and
set `E e_N=-sigma_N e_N`, zero on all other basis vectors.  Then
`||E||<delta`, while the closure of the range of `T+E` misses the `N`th
coordinate.  Consequently an `O(1/R_*)` full-pressure comparison, however
small, cannot on its own transfer density at a fixed radius.  Finite Taylor
matrices remain invertible under such a perturbation only up to the order
whose inverse norm it controls; the example shows exactly why those facts do
not control the infinite approximation ladder.

This refutes the registered bounded-right-inverse implementation of Route B
in the natural `L^2`-source/`C^r`-history spaces.  It does not refute dense
fixed-ring range, unbounded regularized inversion, or an exact spectral
representation of the full Euler operator.

## 5. Failure-derived exact route and minimum repair

Fix the candidate ring and let `A_R` be its actual divergence-free linear
Euler generator, including the whole-space pressure projection.  Let `S_R`
map action amplitudes on a regular Baldi/inner-core action interval into
actual smooth Kelvin initial data, and let `B_R` contain the complete physical
acoustic hybrid, optical angle, `G`/spin/current and required derivative
observations.  The exact source-to-history map is

    C_R a(t)=B_R exp(t A_R) S_R a.                       (8)

For a fixed window its range is dense precisely when the following adjoint
uniqueness statement holds: every compactly supported vector distribution
`Lambda` satisfying

    S_R^* integral exp(t A_R^*) B_R^* dLambda(t)=0       (9)

is zero.  Equation (9) is an exact fixed-radius criterion; pressure is inside
`A_R^*`, and no large-radius comparison remainder is left to dominate an
unbounded approximate inverse.

There are two sufficient ways to earn (9):

1. construct an exact intertwiner
   `A_R S_R=S_R[-i omega(I)]` with a full-row-rank physical gain
   `B_R S_R=G_R(I)` on an open action interval, reducing (9) to the entire-
   transform proof of section 1; or
2. prove (9) directly as an observability/unique-continuation theorem for the
   full pressure-coupled Euler generator, for example by establishing an
   absolutely continuous spectral band of adequate multiplicity and nonzero
   physical gains.

The minimum scientific repair is therefore not another finite moment row.
It is one of these exact full-Euler statements, first for the acoustic and
optical linear outputs and then with a finite kernel correction that solves
the complete phase/energy/current Gram conditions and preserves positive tag
fractions.  The acoustic block must include both parities and the literal
point-to-hybrid acceleration; the optical block must distinguish its positive
material tags from the ambient and have rank for angle and full current.

A high-frequency geometric-optics lift is a plausible new method, but it is
not yet that repair.  Its useful eikonal is provided by the action-angle
frequencies, while the Euler amplitude must solve the pressure-projected
transport/stretch system.  At principal order, writing `k=grad phi` and
`A=grad u`, the actual linearized Euler equation gives

    D_t phi=0,
    D_t k=-A^T k,
    D_t b=-A b+2k [(A b).k]/|k|^2,
    b.k=0.                                                 (10)

The final term is the pressure multiplier forced by the divergence
constraint; omitting it gives a passive-vector amplitude, not Euler.  On a
regular invariant torus (10) is a periodic cocycle, so the admissible clock
frequencies and polarizations depend on its actual Floquet monodromy.  If a
`q`th-order parametrix has residual

    C_(q,T) N^(-q) ||a||_(X_q),                           (11)

then closure requires an explicit diagonal making (11) smaller than the
requested error after the band-synthesis source `a` has been selected.  No
current source bounds the superoscillatory inverse cost against
`C_(q,T)N^(-q)`, and smooth fixed tag gains can themselves decay at high
harmonic order.  Without that joint estimate, geometric optics repeats the
same inverse-versus-remainder gap as Route B.

More precisely, the failure-generated WKB/Floquet route must localize on one
regular invariant torus and solve all three of the following on the actual
finite background:

1. the periodic pressure-projected amplitude equation and its Floquet
   solvability conditions, rather than only the passive eikonal equation;
2. a nonzero, full-rank measured gain for the fixed positive tag and the
   continuous ambient hybrid observation; and
3. a joint norm/error ordering of the form (11), after the finite-window
   inverse cost and every physical observation conditioning constant are
   known.

If the fixed smooth tag annihilates or suppresses the necessary fast radial
packets beyond that ordering, a family such as
`chi_N=chi_0[1+epsilon_N cos(N psi)]`, with
`0<epsilon_N<1`, is a plausible observation repair.  It preserves positive
material fractions, but it changes the tag with preparation scale.  It is
therefore explicitly outside 0250's fixed-tag route unless a later
append-only contract authorizes that scope departure.  It is recorded here
as failure-generated parent progress, not as a license or theorem.

There is a stronger fixed-tag diagonal if the first two items are earned.
Attempts 0112 and 0116 already give the required formal ingredient on an
exact fixed periodic Beltrami field: for every finite order `m`, compact
Kelvin data have a full pressure recursion and exact Euler evolution whose
finite-window error is `C_m N^(-m-1)`.  They also show how a strict elliptic
Euler-amplitude Floquet margin can survive periodic inverse localization.
Their existing shrinking-parcel observation is not the fixed physical tag
required here, and 0116 expressly leaves the resolved optical residue open.

To separate response accuracy from the clock band, choose a fixed angular
combination `ell` with nonconstant
`nu(I)=ell.Omega(I)` and use the exact eikonal

    phi_j=N_j sigma(I)+ell.theta-nu(I)t,                 (12)

where `sigma(I)` is invariant under the base flow.  The temporal band in
(12) is independent of the high WKB scale `N_j`.  After an ideal band
approximant has been selected, its finite source cost therefore does not
change when `j` is increased.

One fixed positive material tag can, in principle, contain all the matching
dual phases.  For bounded smooth dual profiles `q_j`, set

    chi=chi_0[1+epsilon sum_j a_j Re(q_j exp(-iN_j sigma(I)))]. (13)

Choose the WKB orders `m_j` increasing faster than `4j`.  After the finite
coefficient, time-window and derivative constants through level `j` are
known, choose `N_j` recursively so large that

    C_j N_j^(-m_j-1) <= N_j^(-2j),

and take `a_j=N_j^(-j)`, reducing it further to absorb the finitely many
`q_j` derivative norms.  Then (13) is one `C^infinity` tag because
`a_j N_j^k ->0` for every fixed `k`; choose `epsilon` so its absolutely
summed modulation is below one half, making the material fraction strictly
positive.  Its normalized WKB error at the `j`th dual harmonic is at most

    C_j N_j^(-m_j-1)/a_j <= N_j^(-j) ->0.              (14)

Enumerating integer time windows and derivative orders makes the same fixed
tag work on every finite requested window.  Since the band-synthesis cost is
chosen before `j`, (14) can also beat that arbitrary finite cost.  This is an
explicit repair of the inverse-versus-remainder ordering and does not require
a frequency band approaching zero.

Equations (12)--(14) are conditional on the still-load-bearing physical
facts: a stable continuous Floquet branch on an open action interval, smooth
dual profiles with a nonzero full acoustic/optical gain matrix, and a Kelvin
packet whose measured leading coefficient is actually the selected `a_j`
term.  They do not manufacture those gains.  If those conditions fail, the
preparation-dependent `chi_N` family above remains only the separately
declared scope departure.

The exact periodic Beltrami construction in 0145, strengthened by 0147, is
the strongest materially different candidate background.  It gives an exact
stationary periodic constant-curl Euler field with a contractible compact
Euclidean invariant torus, positive-density whole-field law and controlled
global derivatives.  Attempt 0147 additionally supplies a finite-action
actual optical Kelvin packet and one fixed nonnegative registered material
marker on that same kind of exact periodic background, retaining full
pressure, KKS and physical spin.  Its construction still chooses the torus
accuracy and radius after a prescribed finite packet/error inventory; it does
not prove arbitrary response accuracy on one already fixed field, an acoustic
hybrid response, or the gain/adjoint theorem (9).

The independently established 0252 center theorem excludes the 0248B
subclass that places a nonzero elliptic swirl core inside an axisymmetric
pressure-localizable region: at a nondegenerate poloidal center,
`u dot grad p=0` forces the swirl to vanish.  It does not constrain the
nonlocalizable constant-curl 0145/0147 core.  Consequently Route C should
first freeze one exact 0145/0147 periodic field and attack (8)--(14) directly
there.  This can remove the need for a new stationary shape inverse if the
response transfers and the separate geometry owner accepts that field's
compact-Euclidean scope; 0250 does not decide the 0248/0253 transaction.

## 6. Route verdicts and strongest supported scope

- **Route A, exact fixed-ring transport-frequency bands:** blocked at the
  exact projected-Euler intertwiner or equivalent adjoint uniqueness (9), and
  at the full physical gain/Gram calculation.  The straight field-changing
  band and Baldi's fixed-geometry passive band are established positive
  inputs; neither is silently relabeled as the compact-ring full response.
- **Route B, analytic continuation with a bounded right inverse:** refuted as
  stated in the natural smooth band spaces by compactness.  Finite-rank
  continuation and a small curvature operator error cannot substitute for
  dense-range stability, as (7) shows.  Unbounded regularized inversion or an
  exact spectral theorem remains viable.
- **Failure-derived Route C, exact adjoint uniqueness/full Euler spectral
  band on one frozen 0145/0147 periodic constant-curl background:** open with
  (8)--(9) as its precise construction target.  A WKB implementation
  additionally owes its periodic Floquet amplitude and the cost diagonal
  (11).

The strongest new supported statement is the exact arbitrary-finite-window
linear angle controllability of the fixed straight field-changing
constant-curl supplier, together with a necessary-and-sufficient exact dual
criterion for transferring that achievement to one compact ring.  The fixed
compact ring already has actual field-changing finite-window responses and
finite-order full-pressure corrections, but no reviewed source establishes
their arbitrary-accuracy range at one radius.  The compact Euclidean parent
and the distinct 0248 stationary geometry/density obligation therefore stay
open.
