# Degree-minus-two Euler tails: exact channel and evolution classification

## 1. Exterior divergence, flux, and vector harmonics

Write `x=r n`, `|n|=1`, and freeze the asymptotic class

    u(x)=r^(-2)[f(n)n+T(n)]+R(x),                       (1)

where `T dot n=0` and, for some integer `s>=4` and `gamma>0`, the remainder
and its derivatives through order `s` obey

    partial^alpha R=o(r^(-2-|alpha|)),
    <x>^(2+gamma+|alpha|) partial^alpha R in L^2

after subtracting a fixed smooth-core representative. These hypotheses are
stronger than finite kinetic energy; they are used to make the low-frequency
remainder and translated cross-energy remainder uniform.

The spherical divergence formula gives exactly

    div u=r^(-3) div_S T+div R.                         (2)

The radial derivative vanishes because `r^2 u_r=f(n)`. Hence the leading
field is divergence free precisely when

    div_S T=0.                                          (3)

On `S^2` there is no harmonic one-form, so

    f=sum_(l,m) alpha_lm Y_lm,
    T=sum_(l>=1,m) beta_lm n cross grad_S Y_lm.         (4)

The flux through every large sphere is `integral_S2 f`. A globally smooth
source-free divergence-free completion therefore has `alpha_00=0`. The
excluded field `q n/r^2=q x/r^3` obeys

    div(q x/r^3)=4 pi q delta_0,                        (5)

so it is exactly the only rotation-invariant polar tail and exactly the
forbidden source row. Every `l>=1` radial harmonic admits a smooth core:
for a radial cutoff `chi`, add in its transition annulus

    u_T=-(chi'(r)/r) grad_S Delta_S^(-1) f,             (6)

which cancels the cutoff divergence. A toroidal harmonic can be cut off
without correction. Thus (3) and zero mean are also sufficient for a smooth
divergence-free completion.

Under inversion, `Y_lm n` has parity `(-1)^(l+1)` and
`n cross grad_S Y_lm` has parity `(-1)^l`. The source-free tail space is
therefore a direct sum of nontrivial `O(3)` irreducible sectors; it contains
no invariant one-dimensional scalar charge line.

## 2. Exact low-frequency transform

For the convention `fhat(k)=integral exp(-i k dot x)f(x)dx`, homogeneous
distributions satisfy

    FT[r^(-d)Y_lm(n)]
      =c_(l,d)|k|^(d-3)Y_lm(khat),                     (7)

    c_(l,d)=i^(-l) 2^(3-d) pi^(3/2)
             Gamma((l+3-d)/2)/Gamma((l+d)/2).

For the radial tail, write `r^(-2)Y_lm n=x r^(-3)Y_lm` and differentiate the
degree-three transform. For `l>=1`,

    FT[r^(-2)Y_lm n]
      =A_l |k|^(-1) grad_S Y_lm(khat),
    A_l=i c_(l,3).                                     (8)

For the toroidal tail use
`r^(-2)n cross grad_S Y=x cross grad(r^(-2)Y)` to get

    FT[r^(-2)n cross grad_S Y_lm]
      =B_l |k|^(-1) khat cross grad_S Y_lm(khat),
    B_l=c_(l,2).                                       (9)

Both outputs are transverse to `k`, as required. These are total-angular-
momentum vector harmonics: for every `l>=1` the radial and toroidal copies are
two multiplicities of the same `SO(3)` irrep with opposite `O(3)` parity, not
two ambiguous scalar orbital labels. In particular `A_1=pi^2` and
`B_1=-4 pi i`. The radial coefficient `a_R` is polar under `O(3)` and the
toroidal coefficient `a_T` is axial. The two physical `l=1` representatives
and transforms are

    R_aR=(a_R dot n)n/r^2,    Rhat_aR=pi^2 P_T(khat)a_R/|k|,
    T_aT=(a_T cross n)/r^2,   That_aT=4 pi i(khat cross a_T)/|k|. (9a)

Their translated cross-energy multiplicity block is

    K_11(a_R,a_T;b_R,b_T;dhat)=rho_0/|d| times
      [[(pi^3/8) S(a_R,b_R), -2 pi C(a_R,b_T)],
       [-2 pi C(a_T,b_R),       2 pi S(a_T,b_T)]],     (9b)

    S=a dot b+(a dot dhat)(b dot dhat),
    C=dhat dot(a cross b).

The diagonal radial--radial and toroidal--toroidal rows have the same positive
transverse Oseen tensor with different exact weights. The mixed rows are
pseudoscalar, parity-sensitive chiral couplings. This `2 by 2` copy space is
the smallest honest `O(3)`-typed Schur block.

Let carrier one be `u_1(x)` and carrier two be `u_2(x-d)`. Parseval then has
the factor `exp(+i k dot d) F_1 dot conj(F_2)`. Let `F_a(khat)` denote the
angular vector in (8)--(9) and define

    H(khat)=F_1(khat) dot conj(F_2(khat))
           =sum_(L,M) h_LM Y_LM(khat).                 (10)

Since

    FT[r^(-1)Y_LM]=c_(L,1)|k|^(-2)Y_LM,

Parseval and translation give the complete leading cross kernel

    E_12(d)=rho_0/|d| sum_(L,M)
                h_LM Y_LM(dhat)/c_(L,1)+O(|d|^(-2)).  (11)

For the mixed `l=1` row, `c_(1,1)=-2 pi^2 i`, so
`FT^(-1)[i khat/|k|^2]=-dhat/(2 pi^2 |d|)`; this fixes the minus sign in
(9b). Reversing the definition of `d` flips precisely this chiral row.

The stated weighted remainder makes each cross term with the bounded
low-frequency remainder one order faster; smooth high-frequency pieces are
integrated by parts. Formula (11), rather than the power alone, is the charge
candidate: every surviving `L>0` coefficient records separation or internal
orientation.

## 3. The reviewed oriented tail and its sign

For

    phi=(1+r^2)^(-1/2),      u_a=-a cross grad phi
         =(a cross x)/(1+r^2)^(3/2),                   (12)

one has

    phihat=4 pi K_1(|k|)/|k|,
    uhat_a=4 pi i (k cross a)K_1(|k|)/|k|
            =4 pi i (khat cross a)/|k|+O(|k| log|k|). (13)

The inverse transverse Newton kernel is

    FT^(-1)[P_T(khat)/|k|^2]
      =[I+dhat tensor dhat]/(8 pi |d|).                (14)

Therefore

    E_ab(d)=2 pi rho_0/|d|
       [a dot b+(a dot dhat)(b dot dhat)]+O(|d|^(-2)).(15)

For `a=q_1 a_0`, `b=q_2 a_0`, the bracket is strictly positive for every
direction, so `q_1 q_2>0` gives a positive, decreasing `+C/d` kinetic cross
energy. Conditional on `d` being an admissible mechanical coordinate, this
has the repulsive-sign effective-potential convention. Its value changes by a
factor two between separation parallel and perpendicular to `a_0`, so it is a
Coulomb power rather than an orientation-independent scalar charge.

The `l=1` tail in (12) is a polar velocity parametrized by an axial toroidal
coefficient. Three basis copies in one Euler velocity
do not repair isotropy:

    u_(e1)+u_(e2)+u_(e3)=u_(e1+e2+e3).                 (16)

The trace over three independent vector channels would be isotropic, but it
adds independent fields or dynamically orthogonal sectors. A locked internal
frame, structured background, or autonomous fast orientation average is an
additional mechanism and scale, not a consequence of (12).

There is, however, a stronger one-field fixed-frame escape. For real unit
`n`, put

    F(n)=[P_n e_1+i n cross e_2]/sqrt(1+n_3^2).         (16a)

Then `n dot F=0`, `F(-n)=conj F(n)`, and `|F(n)|^2=1` exactly. Choose the UV
form factor `h_sigma(k)=exp(-sigma^2 |k|^2/2)` and define

    uhat_q(k)=q h_sigma(|k|)F(khat)/|k|.                (16b)

It is a real, smooth, divergence-free field in every finite `H^s`. Its inverse
degree-minus-two angular tail is explicit spectrally: expand

    F=sum_(l,m)[p_lm grad_S Y_lm+t_lm khat cross grad_S Y_lm]

and use (8)--(9) backwards,

    U=sum_(l,m)[(p_lm/A_l)Y_lm n
                   +(t_lm/B_l)n cross grad_S Y_lm].    (16c)

Smoothness of `F` gives rapid coefficient decay modulo the polynomial
multipliers. There is no `l=0` term, so (6) gives a source-free smooth-core
completion; the Gaussian Fourier construction supplies one directly. Its
translated cross energy is not merely asymptotic:

    E_12^sigma(d)=rho_0 q_1 q_2
       erf(|d|/(2 sigma))/(4 pi |d|).                  (16d)

This is finite at coincidence, positive for like signs, and Coulombic at large
distance. The construction is still not an SO(3)-scalar. It uses the frame
`(e_1,e_2,e_3)`, and for `R=R_z(pi)`

    F_R(n)=R F(R^(-1)n)=-F(n).                         (16e)

Thus a proper rotation flips the effective sign: differently oriented copies
have a relative-frame-dependent cross kernel. The exact positive atom is a
locked-frame or orientational-charge channel. It becomes an electric-like
scalar channel only after a same-substrate mechanism makes the common frame
physical and controls relative rotations.

This boundary has a general proof. Let
`D(R)F(n)=R F(R^(-1)n)` on the transverse angular Hilbert space. Suppose one
nonzero deterministic tail had the same positive scalar `d^-1` coefficient
for every separation direction and every independently chosen relative
carrier rotation `R`. Averaging the kernel over `dhat` isolates its `L=0`
coefficient, proportional to

    <F,D(R)F>=||F||^2 for every R.                      (16f)

Equality in Cauchy--Schwarz forces `D(R)F=F` for every rotation. The only
rotation-equivariant vector field on `S^2` is `c n`, while Fourier
transversality gives `n dot F=0`, hence `c=0`. Therefore no nonzero single
source-free Euler velocity tail gives a positive scalar law independent of
freely rotatable carrier orientations. This does not cover multiple
independent fields, mixed states, a dynamically locked common frame, or a
controlled autonomous orientation average.

Representation-theoretically, a deterministic scalar-amplitude family
`q -> q U_0` is rotation covariant without an internal orientation only if
`U_0` spans the trivial representation. Section 1 shows that its sole
degree-minus-two polar candidate is the forbidden flux source (5). Nontrivial
irreps can have invariant bilinear pairings, but a two-object law then depends
on their relative representation state unless a physical lock or average is
constructed.

**Routes A and B verdict.** The complete source-free tail space, its radial/
toroidal multiplicity block, and its `d^(-1)` kernel are established. Bare
one-field Euler supports both the oriented `l=1` channel and the exact
fixed-frame isotropic autocorrelation (16d). The latter remains orientation
and frame carrying by (16e), so a deterministic SO(3)-scalar channel is not
yet obtained. Independent-sector, background-lock and autonomous-average
escapes remain active constructions, not exclusions.

## 4. The leading tail is a local Euler asymptotic invariant

First apply an exact stationary filter. For a homogeneous leading field
`u=r^-2[f n+v]`, `div_S v=0`, a degree-minus-four pressure
`p=r^-4 P` can balance the leading acceleration only if

    v dot grad_S f-2 f^2-|v|^2-4P=0,
    nabla^S_v v-f v+grad_S P=0.                        (17a)

Equivalently `curl((u dot grad)u)=0` at degree minus six. In fact the entire
`l=1` radial/toroidal multiplicity has no nonzero stationary member. Write

    u=((a dot x)x)/r^4+(b cross x)/r^3.                 (17b)

If `b!=0`, rotate `b=B e_z` and use its residual axial rotation to put
`a=(A,0,C)`. On `y=0`, one exposing component is

    [curl((u dot grad)u)]_z
      =B(5x^2-2z^2)(Ax+Cz)/(x^2+z^2)^(9/2).           (17c)

Vanishing for every `x,z` forces `A=C=0`; another component then equals
`6B^2 xz/(x^2+z^2)^4`, forcing `B=0`. If `b=0`, rotate
`a=C e_z`; the corresponding component is
`-4C^2 xz/(x^2+z^2)^4`, so `C=0`. This covers both copies and their mixing.
A smooth core cannot repair the nonzero leading homogeneous residual.

The fixed-frame Gaussian realization (16b) also fails an exact global steady
test. Any sufficiently decaying stationary whole-space Euler field has
`M_ij=integral u_i u_j dx` finite. Fourier transforming
`div(u tensor u)+grad p=0`, projecting transversely at `k=epsilon n`, and
taking `epsilon->0` gives

    P_n M n=0 for every n, hence M=c I.                 (17d)

For (16a), Parseval separates the positive radial factor from
`Q=integral_S2 Re(F tensor conj F)`. Exact angular integration gives

    Q=diag(pi(8-pi)/2, pi(-8+3pi)/6, 4pi/3),           (17e)

and `Q_11-Q_22=pi(16-3pi)/3` is nonzero. Thus the minimal radial-form-factor
field (16b) is not stationary. A faster-decaying core correction can change
`M` and is not excluded by (17c), but it must supply the exact missing
isotropic stress. The stronger core-independent test (17a) for the infinite
VSH inverse (16c) remains a named analytic residual; it is not silently
inferred from constant Fourier norm.

Suppose a classical Euler solution on `[0,T]` remains in the weighted class
(1), uniformly with the corresponding derivative bounds. Put
`T_ij=u_i u_j`. Then `T` is integrable, `T=O(r^-4)`, and its first moment has
at worst a logarithmic divergence. With kinematic pressure

    pi=(-Delta)^(-1) partial_i partial_j T_ij,

the Newton-kernel split into `|y|<r/2`, `r/2<|y|<2r`, and `|y|>2r` gives

    pi(x)=partial_i partial_j G(x) integral T_ij(y)dy
           +O(r^-4 log r),                             (17)
    grad pi=O(r^-4),            div T=O(r^-5).          (18)

There is no pressure monopole or dipole of order `r^-1` or `r^-2`: the two
derivatives in (17) are load bearing. Euler consequently gives

    partial_t u=-div T-grad pi=O(r^-4).                 (19)

Multiplying (19) by `r^2` and taking the uniform angular limit proves

    partial_t[f(n)n+T(n)]=0.                            (20)

Thus the full degree-minus-two coefficient is conserved on every common
classical interval on which the frozen affine weighted asymptotic class holds.
Unweighted `H^s` local well-posedness alone does not supply that interval: the
displayed weight lies outside the automatic unrenormalized Riesz-transform
range. Propagation from initial data requires a bespoke estimate for
`R=u-r^-2U`, subtracting the time-dependent stress multipole in (17) and
retaining the logarithmic first-moment remainder. That weighted propagation
theorem remains open here. Conditional on the common interval, (20) is an
exact asymptotic label, not a formal coefficient, but it is vector or
higher-harmonic data rather than scalar charge.

This conservation acts as a boundary-sector label in both directions:
compact or zero-tail data cannot dynamically create a nonzero `V` on the
same weighted classical interval, and Euler conserves only the total
asymptotic coefficient unless a canonical material decomposition assigns
individual carrier labels. Under the exact Euler similarity

    u_(A,B)(x,t)=A u(Bx,ABt),

the tail and its cross coefficient scale continuously as

    V -> (A/B^2)V,              C_(1/d)->(A^2/B^4)C.   (20a)

Neither topology nor conservation selects a universal nonzero magnitude.
The tail is a continuous asymptotic sector unless an additional selector is
constructed.

For (12), direct calculus gives

    curl u_a=[(2-r^2)a+3x(a dot x)]/(1+r^2)^(5/2),     (21)
    div u_a=0,       u_a dot curl u_a=0,
    ||u_a||_2^2=pi^2 |a|^2/2.                          (22)

The energy is finite and helicity is exactly zero. The absolute velocity and
vorticity integrals diverge, ordinary angular momentum diverges quadratically,
and absolute vorticity impulse diverges. Symmetric-ball values of momentum or
impulse are regulator dependent and cannot replace absolute convergence.
The earlier finite-`j` KKS rotation-sphere/action construction therefore does
not transfer to this tail without a new renormalized asymptotic phase space.

## 5. Carrier and interaction verdict

Two translated copies give a legitimate smooth finite-energy Euler initial
datum and the exact initial cross-energy asymptotic (15). Their tails overlap
globally, so support separation does not label two invariant material
components. Neither (15) nor conservation (20) proves a two-carrier moduli
equation, stationarity, orbital persistence, or a force on persistent centers.
The initial pressure contains cross terms everywhere.

There is also an exact coadjoint obstruction. Translating only one noncompact
summand requires a displacement that approaches a nonzero constant at
infinity. The standard relabeling/orbit group used by the compact-carrier
construction is identity at infinity. Since Euler has one total vorticity
field and the two tails overlap, “translate component one while holding
component two” is not automatically a dynamically accessible tangent on one
coadjoint leaf. A mechanical separation coordinate therefore requires either
a finite or renormalized KKS asymptotic phase space with an admissible,
nondegenerate relative-translation tangent, or transported material labels
that canonically split the summands and derive their joint action.

Opposite leading tails cancel in the far field of a neutral pair and can
improve the composite localization, but the external `d^(-1)` label then
vanishes; its internal cross energy remains a finite separation-dependent
quantity. This neutral composite is a distinct faster-decaying global sector,
not two independently visible asymptotic charges. Turning the internal
quantity into a mechanical interaction requires one
joint Euler action evaluated on a controlled persistent translated family,
including deformation and radiation errors.

**Route C verdict.** On every common interval in the declared weighted class,
the degree-minus-two coefficient is an exact Euler asymptotic invariant and
the oriented cross energy is a positive, decreasing Coulomb-power term. It
has a repulsive-sign effective-potential interpretation only if `d` is first
constructed as an admissible mechanical coordinate. This is a stronger
bare-Euler candidate than a mere initial-data counterexample. It remains
blocked as a particle interaction on a persistent two-carrier family, finite
angular/KKS action, and a derived orientation-independent observable. It
neither refutes nor replaces the compact Cao and constrained U(1)
continuations.
