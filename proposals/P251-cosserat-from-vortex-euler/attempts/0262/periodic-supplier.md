# Periodic supplier lemma: explicit Bloch identity and uniform estimates

This is the bounded response to 0264's three finite-K questions. The input
is the same actual smooth periodic steady field and constant-curl source
patch as in `local-normalizer.md`. All constants below may depend on that
fixed field, cell, finite profile inventory, and the already selected
baseline sources. They are uniform in large auxiliary carrier N and in
|K|<=K0, with K0 less than half the shortest nonzero reciprocal vector.
For baseline families having a ray-wise mean, use radial K derivatives on
a compact set of unit rays, retaining that mean matrix exactly.

## 1. The full Bloch product rule

Define div_K=div+iK dot, curl_K=curl+iK cross, and

    [u,xi]_K=(u dot grad+i u dot K)xi-(xi dot grad)u.

For div_K xi=0 and div omega=0, the ordinary product rule gives

    curl_K(xi cross omega)
      =(omega dot grad)xi-(xi dot grad)omega
        -omega div xi+i xi(K dot omega)-i omega(K dot xi)
      =(omega dot grad+i K dot omega)xi
                                      -(xi dot grad)omega.       (1)

The two terms proportional to omega(K dot xi) cancel, rather than being
omitted. On the support of xi, omega=lambda u; all fields and derivatives
of xi vanish outside its compact source patch. Therefore (1) is exactly

    curl_K(xi cross omega)=lambda[u,xi]_K.                       (2)

For xi=curl_K eta, div_K xi=0 identically. The Fourier Helmholtz projector
P_K has multiplier I-(q+K)(q+K)^T/|q+K|^2. Its removed part is parallel
to q+K, so curl_K P_K f=curl_K f for every mode, including q=0 at
nonzero K. Thus v_K=P_K(xi cross omega) obeys (2) with curl_K v_K
on the left. This proves the exact local Kelvin/Lin identity used in
both full forms. It does not replace the Bloch bracket by the ordinary
bracket at nonzero K.

## 2. Uniform homogeneous constraints need no right inverse

Fix M+1 disjoint smooth real bumps G_l in the selected patch. For each N
take the potential eta_N(G)=sigma G p_sigma exp(iN k.x)/(N|k|), using
its real and imaginary parts as the two real potential columns. Because

    f_K=(curl eta_N)cross omega+i(K cross eta_N)cross omega,       (3)

the mean of f_K is affine in K. Impose every real component of the mean
of both potential columns in (3), coefficient by coefficient in K, along
with the finite zero-wave cross-form constraints against prior columns.
These are M real **homogeneous** linear equations A_N c=0 for
G=sum_l c_l G_l. M may be enlarged before choosing the bumps; it is
finite and independent of N and K.

Choose any Euclidean unit vector in ker A_N. Such a vector exists since
there are M+1 columns. For every fixed s,

    ||G||_{C^s} <= (sum_l ||G_l||_{C^s}^2)^(1/2),
    integral w G^2 = sum_l c_l^2 integral w G_l^2.                (4)

Since w has a fixed nonzero sign, the absolute value of the second
expression has a positive lower bound independent of N. The same holds
for the energy weight |k|w^2. Hence both signature and ratio margins,
all envelope derivatives, and the exact zero-mean rows are uniform.

There is no K-dependent kernel choice: all K coefficients in (3) were
imposed at once. There is also no inhomogeneous constraint equation and
no use of A_N's smallest nonzero singular value. A uniformly bounded
right inverse for A_N is unnecessary and may fail as its oscillatory
entries tend to zero. The actual construction only needs (4). This
distinction resolves the review's constraint question without adding an
unsupported singular-value claim. It applies successively to the finite
number of physical control blocks. Smooth dependence on the discrete
carrier N is not used anywhere.

## 3. Full projector bounds, including the remote pressure

For q nonzero and |K|<=K0, |q+K|>=|q|/2. Direct differentiation of the
displayed multiplier gives

    |partial_K^alpha P(q+K)| <= C_alpha |q|^-|alpha|.             (5)

The q=0 force coefficient in (3) is exactly zero for every K; it does
not produce a discarded pressure or mean term. All remaining Fourier
modes are retained in (5). Each complex packet in f_K has the form
exp(iN k.x) times a compact envelope with uniformly bounded derivatives,
plus an O(N^-1) envelope. The K derivative of f_K is O(N^-1), and
its second K derivative is zero. Local compact support makes this a
smooth periodic function even if k itself is not reciprocal-lattice
valued. Integration by parts gives arbitrarily decaying Fourier tails
away from q=N k, uniformly for the profile family (4).

Split the reciprocal lattice into |q-Nk|<N|k|/2 and its complement.
On the first region (5) gives the appropriate inverse power of N; on
the second, arbitrarily many envelope derivatives dominate any fixed
power of q or N. Parseval's identity consequently gives, for every
fixed alpha,

    ||partial_K^alpha v_K||_L2 <= C_alpha N^-|alpha|,
    ||partial_K^alpha curl_K v_K||_L2
                                 <= C_alpha N^(1-|alpha|).      (6)

The real packets at -Nk obey the same estimates. These estimates are
for the full periodic projector, so its image pressure tails have not
been approximated by a local projector. Spatial Sobolev derivatives
multiply their bounds by the corresponding fixed power of N.

For two unnormalized auxiliary columns, (6) and the full H pairing give

    |partial_K^alpha H_c| <= C_alpha N^(1-|alpha|).              (7)

The local Omega pairing has bound C_alpha N^-|alpha| through order two
and is exactly quadratic in K. Actual whitening and the separated
phase/energy ratio in `local-normalizer.md` multiply energy columns by
O(N^-1/2) and leave the combined phase columns O(1). Their coefficient
matrices have uniform bounds. Thus (7) is uniform for every derivative
of order at least one after the exact zero-wave normalization. The
zeroth-order Gram matrices are exactly the fixed block matrices there.
It follows by the fundamental theorem of calculus that

    ||H_c(K)-H_c(0)||+||Omega_c(K)-Omega_c(0)|| <= C|K|,          (8)

with C independent of large N. This is the fixed neighborhood consumed
by the finite contraction, not a neighborhood shrinking with N.

## 4. Baseline cross forms and their derivatives

Freeze a finite baseline family smooth in space and in the retained K
parameters. Its norms may be large; they are recorded before choosing N.
The exact stationary Hessian with its second argument supported in the
constant-curl patch is

    H(B,C)=rho integral conjugate(v_B)dot
                                  (v_C-curl_K v_C/lambda).

The Bloch curl is self-adjoint for real K. Integration by parts and
self-adjointness of P_K move all derivatives and the projector onto
the smooth baseline test:

    H(B,C)=rho integral
       conjugate(v_B-curl_K v_B/lambda)dot f_C.                  (9)

The integral in (9) is now over the compact auxiliary force support.
The pressure in v_B is its actual full pressure, not a local surrogate.
Every retained K derivative differentiates a smooth test or (3).
Integration by parts against exp(iNk.x) proves an O(N^-Q) bound for
each prescribed Q, with the corresponding finite baseline derivatives
in the constant. The local KKS cross form has the same proof with
omega and B as the test. This supplies uniform bounds for all X_H,X_O
derivatives used in the affine Gram equation. It also explains why
deleting those cross terms by support separation would have been wrong.

The finite-K Kelvin product identity, homogeneous profile construction,
and full Fourier estimates (1)-(9) supply the three missing bounds named
in 0264. The remaining smoothing-to-observation and current consumer
retain their original hypotheses; these estimates do not construct the
compact stationary background.

`route_verdict: established as stated for the explicit periodic
auxiliary supplier and its carrier-uniform finite-K form bounds`

`evidence_scope: exact differential identities, homogeneous kernel
bounds and full periodic Fourier estimates; correction review pending`
