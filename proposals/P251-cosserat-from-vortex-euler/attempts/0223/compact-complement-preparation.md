# Constructive compact-pressure preparation and its precise remaining cost

This is a continuation of the actual remainder (9) in
`observed-second-jet.md`, not an assumption of optical closure. It gives
a quantitative finite-energy spectral construction at a specified norm,
and identifies which extra estimate a smooth second-jet diagonal needs.

## 1. The actual input kernel is regular after the exposed derivative

Fix a closed regular streamline annulus inside the same nonlinear cell,
away from the center and separatrix. Its input energy interval is compact,
omega'(E) is nonzero there, and l=4. All spatial cutoffs below are fixed
smooth annulus cutoffs. No tag or tag density is narrowed in this step.
For s=f(E)exp(il theta_o), split the actual R of (7) as

    R s=R_a f'+R_b f.                                (1)

Both R_a and R_b are order-minus-two maps on the one-dimensional input;
this statement includes the explicit elliptic factors, not just a word
"compact". For instance

    R_a f'=H^-1[-(grad psi.grad E) exp(il theta_o) f'],

and the other terms contain either this same H^-1 with a smooth fixed
coefficient or H^-1 times a differential operator of order at most two
times a second H^-1. The other actual forcing H^-1 T(psi s) has the
same order. Mean removal is a smooth finite-rank term.

Here is a direct regularity argument adequate for the spectral estimate.
Take smooth cutoff input modes exp(in E). Since |grad E| is bounded
below on the annulus,

    H[e^(inE) a]=e^(inE)[n²|grad E|² a
                      -in(2grad E.grad a+Delta E a)+Ha].

Successively canceling powers of n constructs a fixed-annulus elliptic
parametrix for H^-1 with leading amplitude a/(n²|grad E|²).
The residual can be made O(n^-m) for any prescribed finite m. The
global mean-free inverse is bounded, so each actual input coefficient
has output L² norm O(n^-2). Applying T any fixed number of times does
not change this radial order: T E=0, and the differentiated parametrix
has the same estimate. Differential factors in R_b can be treated
between its two parametrices. This also proves the estimate for their
global tails, rather than replacing the periodic pressure locally.

Consequently their Hilbert-valued input kernels are H^s_E for every
s<3/2, in particular H^1_E, both in L² and in any prescribed finite
T-graph norm. The smooth diffeomorphism nu=l omega(E) transfers H^1_E
to H^1_nu. The norm of f' in (1) is an actual input cost, and remains
in every estimate below.

The order-zero term psi s lacks this kernel regularity and was already
separated into the explicit secular response in (8). Applying the
following inverse argument to that term would be a false shortcut.

## 2. A quantitative inverse needing no spectral-density hypothesis

On the mean-free normal torus put B0=1-H^-1. Its only zero Fourier
eigenspace is the complete first shell. The actual transformed planar
generator

    G=-sqrt(B0) T sqrt(B0)

is skew-adjoint on its transport domain: it is -T plus a bounded
skew-adjoint order-minus-one perturbation (including its finite-rank
first-shell terms). This is a property of the actual stationary Euler
operator. Its positive quotient norm is NOT its physical Jacobi energy.

For any self-adjoint A and a Hilbert-valued input kernel X(nu) in H^1
on a fixed interval, extend it smoothly with compact support and write
X=sum_m x_m exp(im nu/c). Its coefficient norms are summable by
Cauchy--Schwarz. The smooth spectral multiplier

    y_epsilon(nu)=(A-nu)[(A-nu)²+epsilon²]^-1 X(nu)

has inverse norm at most 1/(2epsilon). Its residual is
r_epsilon(A-nu)X, where r_epsilon(d)=epsilon²/(d²+epsilon²).
For each fixed x_m, the spectral theorem and Tonelli give

    int ||r_epsilon(A-nu)x_m||² dnu
                  <=(pi epsilon/2)||x_m||².

Minkowski then yields the explicit estimates

    ||residual||_(L²_nu;Hilbert)
          <=sqrt(pi epsilon/2) sum_m ||x_m||
          <=C sqrt(epsilon)||X||_(H^1_nu;Hilbert),
    ||y_epsilon||_(L²_nu;Hilbert)
          <=||X||/(2epsilon),
    ||partial_nu y_epsilon||
          <=C[epsilon^-2||X||+epsilon^-1||X'||].    (2)

These inequalities hold with arbitrary point, singular and continuous
output spectrum. The non-atomic variable is the actual smooth input
band nu(E); an assumption of absolutely continuous Euler output
spectrum has not been added.

## 3. Reconstructing the true scalar Kelvin state

For the scalar complement L_z=-B0 T, consider

    (-i nu+B0 T)z=g,  |nu|>=nu_min>0.               (3)

Solve the transformed equation approximately,

    (-i nu-G)q=sqrt(B0)Tg+e,
    z=(sqrt(B0)q-g)/(i nu).                         (4)

Direct substitution gives residual sqrt(B0)e/(i nu) in (3).
This exact reconstruction does not invert T, delete its streamline
kernel, or discard the forced first Fourier shell. The component -g/
(i nu) explicitly retains those rows. Its T-graph norm follows from
(4), the transformed equation, and the bounded first-shell projection
of T; higher fixed T-graph norms follow by the same identities.

Apply (2) to the actual kernels in section 1, with A the self-adjoint
multiple of G and the matching factor i. This constructs approximate
initial complement kernels with explicit polynomial inverse/residual
costs in these norms. In a source f_h of width h normalized to a
nonzero observed amplitude, the elementary costs ||f_h||=O(h^-1/2)
and ||f_h'||=O(h^-3/2) remain. For example a single residual estimate
is C_T sqrt(epsilon) h^-3/2 on a fixed time window; it is not declared
small by comparing it only with the much larger leading field norm.

The t-linear forcing in (9) is handled by a triangular pair of these
Sylvester equations. If w=Z0 s+t Z1 s and S(Z)=L_z Z+ZT, they are

    S(Z1)=-i H^-1 T psi,
    S(Z0)=Z1-i R.                                  (5)

When R is represented using f', its input evolution includes the real
commutator (partial_E f)_t=-i nu partial_E f-i nu' f. The additional
nu' term belongs to the same finite triangular system; it is not
silently commuted through the inverse. Successively chosen epsilon
values give finite polynomial graph-norm costs and arbitrarily small
residuals for each fixed smooth source. This is an actual prepared
Euler complement, not an independently imposed oscillator.

## 4. Physical output and the norm question resolved by the next source

Neither (2) nor (5) says that the corrected theta/G/S rows have the
required common harmonic. Their constant and secular coefficients
must still be substituted into (3)--(4) of the observed-jet source.
Band-center and amplitude derivatives can change certain physical
coefficients; equality of all three rows requires their actual rank,
not a generic controllability assertion. The extra second-order flux
continues to be present.

There is also a precise norm distinction. The estimates above control
L² and transport graph norms of the scalar complement. They do not
by themselves control a full spatial Sobolev norm of grad z, the
actual material generator. The global transport has hyperbolic
streamlines, and its spectral inverse has not been shown to admit a
polynomial spatial regularization cost. Smooth graph-core
approximations give smooth finite-cost initial data for every fixed
source and error. Their possibly growing full norms remain explicit
unknown constants, not a uniform finite-wave diagonal theorem.

The next source `complete-action-transport-norm.md` resolves this cost
question for the AXIAL linear supplier by changing to its actual
pressure/phase/energy representation: full spatial norms are not needed
for those estimates. It proves cancellation of the large gradient
terms and bounds the actual forms and response by the polynomial
transport norms above. This is a repair of the conservative sufficient
norm requirement, not a newly proved uniform bound in spatial H^s.

The remaining constructive choices are therefore concrete: prove the
needed physical-output rank with the full current row, or construct
the finite observed rows directly from actual smooth initial states
without a full graph inverse. The latter remains within the registered
output route and does not require stationarity of the entire complement.
