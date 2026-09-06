# Fixed-frame stationary filter and affine weighted propagation

## 1. Two independent questions

The `0072` fixed-frame symbol gives a smooth finite-energy Euler initial field
and an exact positive `d^(-1)` cross-energy atom. Two questions must remain
separate:

1. does its complete homogeneous `r^(-2)` angular field solve the stationary
   Euler sphere equations; and
2. does a smooth source-free `r^(-2)` asymptotic datum persist under the local
   Euler flow in a topology that remembers it?

A faster-decaying core can change a global stress moment, but cannot change a
nonzero homogeneous degree-minus-six stationary residual. A nonstationary tail
can nevertheless be a persistent boundary datum. Route C below closes the
second question through vorticity transport and a zero-mean Biot--Savart
reconstruction. It does not use a derivative-losing velocity/pressure
bootstrap.

## 2. Dyadic weighted spaces and the affine class

Fix integer `s>=4`, `0<alpha<1`, and `0<gamma<1`. For `m>0`, define

    ||h||_(Z_m^(k,alpha))
      =sup_(R>=1) [sum_(j=0)^k R^(m+j)||nabla^j h||_(L-infinity(A_R))
          +R^(m+k+alpha)[nabla^k h]_(C^alpha(A_R))],    (1)

where `A_R={R<|x|<2R}`, plus the ordinary `C^(k,alpha)` norm on `|x|<2`.
This records pointwise derivatives explicitly; bare `H^s` regularity is not
used as a weighted substitute.

Let

    U(n)=f(n)n+v(n),    div_S v=0,    integral_(S^2)f=0. (2)

The exact cutoff and annular divergence correction from `0072` gives a smooth
global divergence-free representative `U_*` satisfying

    U_*(rn)=r^(-2)U(n) for r>=R_*,
    U_* in H^k(R^3) for every finite k.                 (3)

Set `Omega_*=curl U_*`. The affine phase space is

    A_U^(s,alpha,gamma)
      ={u=U_*+R:
          div R=0,
          R in H^(s+1),
          q=curl R in Z_(3+gamma)^(s-1,alpha) intersection L^1,
          integral q=0}.                               (4)

The last row follows automatically for a smooth `R=O(r^(-2-gamma))`, but it
is retained because it is the low-frequency hypothesis in the reconstruction
lemma. On all of `R^3`, the `L^2` divergence-free field with curl `q` is unique,
so

    R=B q=curl(-Delta)^(-1)q.                           (5)

Changing `U_*` inside a compact set changes (4) by a compactly supported smooth
translation of the affine coordinate and leaves the asymptotic datum fixed.

## 3. Zero-mean weighted Biot--Savart reconstruction

Let `L(x)` be the matrix kernel of `B`; it is homogeneous of degree `-2`.
For

    q in H^s intersection Z_(3+gamma)^(s-1,alpha) intersection L^1,
    div q=0,
    integral q=0,                                      (6)

one has

    Bq in H^(s+1) intersection Z_(2+gamma)^(s,alpha),
    ||Bq||_(H^(s+1) intersection Z_(2+gamma)^(s,alpha))
       <=C[||q||_(H^s)+
            ||q||_(Z_(3+gamma)^(s-1,alpha))+||q||_(L^1)].   (7)

This is not an unrestricted inhomogeneous Sobolev gain. Use the cancellation

    Bq(x)=integral [L(x-y)-L(x)]q(y)dy.                 (8)

For `r=|x|`, split into `|y|<r/2`, `|x-y|<r/4`, the remaining comparable
region, and `|y|>2r`. On the first region,

    |nabla^j L(x-y)-nabla^j L(x)|
       <=C r^(-3-j)|y|,                                (9)

while

    integral_(|y|<r/2)|y||q(y)|dy<=C r^(1-gamma).      (10)

Their product is `O(r^(-2-gamma-j))`. The ball about `x` gains one derivative
by the standard odd-kernel cancellation after rescaling to a unit annulus;
the weighted `C^alpha` norm supplies the endpoint. The other two regions give
the same or faster order directly. Differentiating (8) and repeating the split
proves every row of (7), including the scaled Hölder seminorm.

The global Sobolev row includes the low frequency rather than hiding it in an
inhomogeneous order-minus-one slogan. From zero mean,

    |qhat(k)|<=integral min(2,|k||y|)|q(y)|dy
              <=C |k|^gamma ||q||_(Z_(3+gamma)^0 intersection L^1) (10a)

for `|k|<=1`; split at `|y|=|k|^(-1)`. Hence
`|k|^(-1)qhat(k)` is square integrable at the origin, while the `H^s` row
controls high frequency. This proves the `H^(s+1)` part of (7). Since
`div q=0`, Fourier inversion gives

    curl Bq=q,    div Bq=0.                            (10b)

Conversely, if `R in L^2` is divergence free and `curl R=q`, its Fourier
transform agrees almost everywhere with `i k cross qhat/|k|^2`; thus `R=Bq`.
This proves uniqueness in the declared global energy class.

Thus (4) is equivalently the affine velocity statement

    u(rn)=r^(-2)U(n)+O(r^(-2-gamma))                    (11)

with the derivatives through order `s` carried by (1).

## 4. Vorticity transport closes the weighted class

Let `omega=curl u=Omega_*+q`. The exact vorticity equation gives

    partial_t q+u dot nabla q-q dot nabla u
       =-u dot nabla Omega_*+Omega_* dot nabla u.       (12)

There is no pressure in (12). Under the bootstrap (4),

    u=O(r^(-2)),       nabla u=O(r^(-3)),
    Omega_*=O(r^(-3)), nabla Omega_*=O(r^(-4)),         (13)

so the right side of (12) is `O(r^(-6))`, strictly faster than the
`Z_(3+gamma)` weight. For a multi-index `|a|<=s-1`, commute `D^a` through
(12):

    (partial_t+u dot nabla)D^a q
      =-[D^a,u dot nabla]q+D^a(q dot nabla u)+D^a F,
    F=-u dot nabla Omega_*+Omega_* dot nabla u.        (13a)

On each transported annulus, `u=O(r^-2)` keeps radii comparable on a common
short interval. Multiplying (13a) by `r^(3+gamma+|a|)`, each Leibniz term puts
one derivative of `u` in `L-infinity` and the other factor in the matching
weighted row. For `|a|=s-1`, the `C^alpha` difference quotient obeys the same
estimate: the top commutator contains `nabla u` times the top quotient plus
lower weighted products. The reconstruction (7) controls all velocity
coefficients. Thus

    d/dt ||q||_(Z_(3+gamma)^(s-1,alpha))
      <=C[1+||u||_(H^(s+1))+Q]
             [Q+1+||u||_(H^(s+1))].                  (13b)

The `L^1` equation has the same product bound. Its zero row is exact because

    integral[-u dot nabla Omega_*+Omega_* dot nabla u]=0, (13c)

componentwise, by `div u=div Omega_*=0` and the declared decay; the transport
and stretching terms have zero integral by the same integration by parts.
Combining these estimates gives

    Q(t)<=C Q(0)+C integral_0^t
       [1+||u(tau)||_(H^(s+1))+Q(tau)]^2 d tau,         (14)

where

    Q(t)=||q(t)||_(Z_(3+gamma)^(s-1,alpha) intersection L^1). (15)

The direct curl identity provides a second check: once (7) is known,
`q=curl(u-U_*)` and the boundary integral of `u-U_*` is `O(r^(-gamma))`.

The approximation and common interval are explicit. Choose radial cutoffs
`chi_N=1` on `|x|<=N`, supported in `|x|<=2N`, and set

    q_N=curl(chi_N R),    R_N=Bq_N.                   (15a)

Smooth by convolution on a scale tending to zero if necessary. Every `q_N`
is a compact curl, hence divergence free and exactly zero mean. On the
transition annulus, the cutoff term `nabla chi_N cross R` and its `j`th
derivative are `O(N^(-3-gamma-j))`; therefore the norms in (6), including
`L^1`, are uniformly bounded. Since `R in H^(s+1)`,
`R_N=P_L(chi_N R)` converges to `R` in `H^(s+1)`; the mollified sequence has
the same limit and uniform bounds.

Standard `H^(s+1)` Euler theory applied to `U_*+R_N` now supplies one lifespan
depending only on the uniform `H^(s+1)` bound. Estimates (13b)--(14) are
independent of `N`. On every fixed compact set, Arzela--Ascoli for the
uniformly Hölder top derivatives, together with standard Sobolev stability,
identifies the limit with the unique `H^(s+1)` solution. Lower semicontinuity
on each dyadic annulus and then the supremum preserve the weighted bound; the
`L^1` and zero-mean rows pass using (13c) and the uniform tail estimate. This
closes the bootstrap from the initial datum itself. The transport argument has
no top-order pressure derivative and no infinite reserve of weighted
derivatives. This interval is common to the approximating sequence and is
determined by the chosen datum's norms; it is not uniform over the entire
affine class.

**Affine invariance theorem.** Every `u_0 in A_U^(s,alpha,gamma)` generates a
unique classical Euler solution on some `[0,T]` with

    u in C H^(s+1) intersection C^1 H^s,
    u(t) in A_U^(s,alpha,gamma)                         (16)

uniformly on that interval. By (7),

    lim_(r->infinity) r^2 u(t,rn)=U(n)                 (17)

uniformly with all angular derivatives licensed by (1). The complete radial
and toroidal degree-minus-two coefficient is therefore an exact local Euler
invariant on the constructed affine phase space. The theorem constructs the
weighted interval from initial data; it does not assume one.

## 5. Independent pressure multipole and time derivative

The pressure calculation supplies a second, independent check of the frozen
coefficient and the exact far-field force density. Put `T_ij=u_i u_j`. From
(11),

    T in L^1,       T=O(r^(-4)),
    M_ij(t)=integral T_ij(t,y)dy.                       (18)

For `G=1/(4 pi r)` and `K_ij=partial_i partial_j G`, the same four-region
kernel split as in (8)--(10) gives

    p(x)=K_ij(x)M_ij+O(r^(-4)log r),
    nabla p=O(r^(-4)).                                 (19)

The logarithm is the genuine borderline first moment of an `r^(-4)` stress.
No weighted top-derivative propagation is inferred from (19); that work was
done in vorticity variables. Since `u dot nabla u=O(r^(-5))`, Euler yields

    partial_t u=O(r^(-4)),                              (20)

consistent with (17). The reusable exact kernel is

    K_ij=[3x_i x_j-r^2 delta_ij]/(4 pi r^5).            (21)

## 6. Consequences and remaining stationary route

The `0072` fixed-frame Gaussian field belongs to the affine class associated
with its inverse-VSH coefficient. Its tail is therefore locally persistent
even though the minimal Gaussian core is not stationary. Compact or
faster-decaying data lie in the `U=0` sector and cannot create a nonzero tail
on the interval (16).

Euler preserves the total asymptotic coefficient, not a canonical split of
two overlapping tails. Similarity acts continuously by

    U ->(A/B^2)U.                                      (22)

No magnitude or action quantum is selected. The positive decreasing `C/d`
cross energy acquires a repulsive-sign effective-potential interpretation only
after `d` is constructed as an admissible relative coordinate in a finite or
renormalized KKS phase space, or material labels supply a joint action.

The fixed-frame stationary test remains active. Its inverse coefficients
carry `i^l`, so the inverse is a rotational spectral multiplier rather than a
classical pseudodifferential symbol. Split `l mod 4`, factor the elliptic
magnitude from the quarter-period phase, and evaluate the exact Shvydkoy
sphere residual. A rigorously tail-bounded nonzero projection refutes this tail
as steady and no core can repair it. Exact cancellation instead activates the
global isotropic-stress repair and full steady gluing problem.

**Route verdicts.** Route C is established by (7), (12)--(17). Routes A/B
remain active on the fixed-frame stationary residual. Nothing here licenses
all-time persistence, a stationary carrier, Euler force, scalar charge, action
selection, P5, electron, neutrino, or parent completion.
