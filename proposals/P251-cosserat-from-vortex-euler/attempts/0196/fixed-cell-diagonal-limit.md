# Actual fixed-cell diagonal physical-action limit

This is the registered continuation after the finite-window construction.
The background Euler cell, A/B=1/100 and positive a are FIXED. Only the
actual prepared microscopic V velocity/configuration changes with the
chosen accuracy. Norm growth is retained before the macro k is selected.

## 1. Fixed-cell Bloch derivatives have actual finite constants

On the body-frame torus of period2pi, retain the physical transverse
mean plane at Fourier mode n=0. Its Leray multiplier is P_kappa, independent
of radial k. For n!=0,

    P_k(n)=Id-(n+k kappa)(n+k kappa)^T/|n+k kappa|².

For |k|<=1/4 its derivatives obey

    ||partial_k^j P_k(n)||<=10 j! 4^j |n|^(-j).      (1)

One proof uses its two complex denominator roots, each of modulus |n|.
On the complex disk of radius |n|/4 about such a real k, |z|<=|n|/2,
the denominator has magnitude at least |n|²/4 and the numerator norm is
at most9|n|²/4. Cauchy's derivative estimate gives(1). The n=0 derivative
is zero. Thus no changing inverse-cell-period constant is present.

To differentiate on a fixed divergence bundle, use at each nonzero mode
the minimal orthogonal rotation Q_k sending n/|n| to
(n+k kappa)/|n+k kappa|. Its explicit Rodrigues expression has denominator
1+their scalar product, bounded away from0 for |k|<=1/4. On the mean
transverse plane Q_k is identity. It is an H^s isometry for real k, with
the same type of derivative bounds as(1).

The exact generator on this fixed bundle is
Q_k^(-1)L_K Q_k. The full Euler H^s energy estimate gives

    ||exp(tL_K)||_(H^s to H^s)<=exp(Gamma_s |t|),

where Gamma_s is a fixed Sobolev product/commutator constant times a
finite W^(s+1,infinity) norm of the FIXED trigonometric background u.
Its transport term, including ik(kappa.u), is skew in the L² estimate;
the strain and Sobolev commutators supply Gamma_s. The fixed-bundle
generator's first three k derivatives have finite H^(s+1)->H^s norms
B_(s,j), obtained from(1), Q_k and this same fixed u.

Three differentiations of Duhamel's formula consequently bound the
third k derivative of the actual Euler velocity history by a finite
polynomial in T,B_(s,j),exp(Gamma_(s+3)T), times its initial H^(s+3)
norms and their first three k derivatives. The actual Lin displacement
equation is another first-order transport/strain equation with this
velocity as its source and satisfies the corresponding bound. These
are full-pressure Euler/Lin estimates, not a finite Fourier truncation.

For the requested finite time-derivative order r, take an integer s>=8+r
(the additional r derivatives absorb the actual transport losses) and define

    N_n=1+||g_n||_(H^(s+4))+||h_n||_(H^(s+4)).

The exact finite initial preparations are the0188 curl lift plus its
D correction and, for V, V+ik P_K(f g_n e_X) with initial material
configuration ik P_K(f h_n e_X). All mean conditions are kept by Lin's
formula. Their third k derivatives have bounds linear in N_n. Phase and
energy are quadratic forms of these actual fields, so one may take

    C_n=C(s,T,u,rho,a) N_n²,                         (2)

with C a finite constant formed from the preceding derivative and
Sobolev-product bounds, enlarged to cover the physical mean and required
time derivatives. The normalized third-order remainders of all these
linear responses and quadratic forms are bounded by C_n |k|³. The
constant can be enlarged to include a specified finite collection of
actual smooth tag observations with their fixed nonzero reference
normalizations; no optical tag or inertia is supplied by doing so.

For clarity, the constants in(2) can be obtained without an unspecified
uniformity assertion. Let Gamma be the maximum Euler/Lin energy constant
through the chosen three-derivative Sobolev hierarchy, B_j the corresponding
operator derivative bounds, and d_j N_n the explicit initial-data derivative
bounds from P_k,Q_k. Set

    p0=d0,
    p1=d1+B1 d0 T,
    p2=d2+(2B1 d1+B2 d0)T+B1²d0 T²,
    p3=d3+(3B1 d2+3B2 d1+B3 d0)T
             +(3B1²d1+3B1 B2 d0)T²+B1³d0 T³.

Duhamel induction bounds the jth state derivative by
exp(Gamma T) p_j N_n. The constants B_j and d_j are obtained from the
explicit multiplier/rotation derivatives and the finite four-mode
Fourier support of u, by the usual weighted convolution inequality.
If F_j bounds the jth derivative of one of the actual Jacobi quadratic
forms on that Sobolev space, its third derivative is bounded by

    exp(2Gamma T) N_n²
      sum_(j+l+m=3) [3!/(j!l!m!)] F_j p_l p_m.

Taking the sum/max of these displayed quantities for the finite list of
physical mean, time-derivative, phase and energy rows gives C_n. Every
F_j is read from the explicit Euler/Lin action; no field-dependent
favorable sign or empirical comparator enters this bound. This also
shows exactly where growth of g_n,h_n enters the actual remainder.

The norms N_n may grow extremely quickly as streamline bands approach
the saddle. Equation(2) does not claim a useful uniform growth rate; its
scientific role is that each fully selected finite control has an ACTUAL
finite bound before the macro scale is chosen.

## 2. Choose the macro scales after those constants

Choose delta_n=2^(-n), construct the passive/configuration controls from
the actual Euler target at this accuracy, and record their N_n,C_n.
Select successively

    0<k_n<min(k_(n-1)/2,1/4,
              2^(-n)/(1+N_n+C_n)).                  (3)

Every choice is positive. Then k_n N_n and k_n C_n tend to zero.
The microscopic displacement/velocity preparations k_n h_n,k_n g_n
tend to zero in the selected Sobolev norm, while their derived finite
second-order response remains nontrivial. The complete normalized
higher spatial remainder C_n k_n³/k_n² also tends to zero.

On this diagonal, the ACTUAL whole-law initial phase and energy satisfy

    Omega0=rho J+o(k_n²),
    H0=rho[|V|²+a k_n²|D|²]/2
          +o(k_n²)(|D|+|V|)².                       (4)

The physical mean columns and their required time derivatives obey

    F_n=1-a k_n²t²/2+o_(C^r)(k_n²),
    G_n=t-a k_n²t³/6+o_(C^r)(k_n²).                 (5)

Thus the exact moving physical phase chart gives

    M_n=rho+o_(C^(r-1))(k_n²),
    (M_n)_t=o_(C^(r-2))(k_n²),
    K_n=rho a k_n²+o_(C^(r-2))(k_n²).               (6)

These are a positive ordinary acoustic action's fixed-time second-order
asymptotics on ACTUAL prepared Euler/Lin histories, including its
conserved physical energy. A nonzero finite-n time connection was kept
until it tended to zero; it was never declared autonomous by fiat.
The exact raw packet momentum current(10) of the configuration proof
remains in the inherited action before its full restriction.

## 3. Precise continuum meaning and remaining physical joining

One may define a full radial family by using control n on
k_(n+1)<|k|<=k_n. Equations(4)–(6) still give a second-order asymptotic
expansion. That family need not be C² in k at its preparation-switch
points. This is therefore NOT a claim about differentiating a single
fixed-preparation Bloch family twice, nor a uniform microscopic C²
Bloch norm. Every finite preparation has its own genuine derivative
estimate(1)–(2), and only then is the diagonal taken.

If desired, select the nth control on the finite window |t|<=n and use
the corresponding C_n in(3). The same result is locally uniform on every
fixed compact time interval. It does not yield an acoustic-time window
t~1/k_n: that window can grow much faster than n, and no such estimate
has been assumed.

The first-cell response on a disjoint transported elliptic-core tag is
zero for every control. The condition k_n N_n->0 also controls the
normalized next local-flow remainder, but does not itself construct an
optical packet or its physical mixed phase/current rows. The whole-law
isotropic physical mean coefficient and its actual inherited phase are
now positive continuum inputs. Same-EPS optical insertion, hybrid
centroid/ambient observations, full physical spin/current coupling and
any acoustic-time requirement remain the exact parent joining work.
