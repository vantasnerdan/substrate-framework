# Exact isovortical channel and its measured gain rows

This continuation checks the proposed divergence-free Lin generator on the
exact constant-curl action-angle region of the fixed 0211 ring. It repairs
the initial-data obstruction in `neutral-polarization.md`; it does not assume
that the resulting physical observation matrix has nonzero rank.

## 1. Exact divergence-free generator

Let `curl u=omega=lambda u`, with constant nonzero `lambda`, on a regular
invariant action-angle annulus. Write its volume coordinates as
`(I,theta_1,theta_2)`, with volume density `J(I)`, and

    u=Omega(I).partial_theta,      k=grad I,      k.u=0.       (1)

For a smooth compactly action-supported amplitude

    f=A(I) exp(i N I) g(I,theta),                              (2)

put

    a=-k cross u/|k|^2,
    xi_perp=-k f/(lambda |k|^2).                              (3)

The sign and normalization are exact:

    xi_perp cross omega
      =-[k f/(lambda |k|^2)] cross (lambda u)=a f.            (4)

Adding `beta u` does not change (4). Since `div u=0`, the condition
`div xi=0` for

    xi=xi_perp+beta u                                         (5)

is precisely the cohomological equation

    u.grad beta=-div xi_perp=:h.                              (6)

Its solvability condition holds for every zero-mean angular harmonic. The
normal flux is

    dI(xi_perp)=-f/lambda.                                    (7)

Periodic integration of the coordinate divergence therefore gives

    integral_T2 J h dtheta
      =-partial_I integral_T2 J dI(xi_perp)dtheta
      =(1/lambda)partial_I[J A exp(iNI) integral_T2 g dtheta]
      =0.                                                     (8)

Angle components of the Euclidean normal vector may be present in a general
action-angle gauge; they are periodic divergences and do not alter (8).
Thus, on a torus where the zero-mean cohomological inverse exists,

    beta_m=h_m/[i m.Omega(I)],       m!=0,                    (9)

defines a smooth solution after fixing its torus mean to zero. For a
Diophantine torus, `|m.Omega|>=gamma |m|^-tau` gives the explicit Sobolev
loss

    ||beta||_Hs <= gamma^-1 ||h||_H(s+tau).                  (10)

Because `h` contains one action derivative of (2), its fixed-order norms
grow polynomially; in particular `||beta||_C0=O(N)` and, for each fixed
spatial order `s`, `||beta||_Cs<=C_s N^(s+1)`. These finite costs must be
selected before the later long-wave scale.

There is a stronger open-annulus conclusion on the exact axisymmetric 0211
ring than a single Diophantine-torus statement. Choose a fixed nonzero
toroidal harmonic `n`. Axisymmetry preserves `n` in `h`, although its
poloidal Fourier index `m` need not remain single. With
`rho(I)=Omega_phi/Omega_theta`, all denominators have the form

    m Omega_theta+n Omega_phi=Omega_theta[m+n rho(I)].        (11)

On a compact inner annulus `rho` is bounded, so only finitely many integers
`m` can make (11) small. Strict twist makes each corresponding resonance
isolated. Reducing to a nonempty open subannulus away from this finite set
gives one uniform inverse for every poloidal Fourier coefficient. Hence
(5)--(9), not merely their value on one Cantor torus, are smooth on an open
action interval. Choose a pair `(m_0,n)`, with `n!=0`, for which

    nu(I)=m_0 Omega_theta(I)+n Omega_phi(I)                   (12)

is nonconstant; such a pair exists because strict rotation-number twist
precludes both frequency components being constant. On a still smaller
interval, `nu` is one-to-one and avoids zero. This supplies the nonzero open
transport band required by the 0254 bridge on one already fixed ring.

The exact initial Euler velocity is now

    v_N(0)=P(xi cross omega)=P(a f).                          (13)

It is field-changing, not a passive tracer. Its leading curl is
`iN k cross a f=iN u f`; consequently it is nonzero for every sufficiently
large carrier for a nonzero envelope. Equations (5)--(13) replace the
singular total-phase inversion (12)--(13) of `neutral-polarization.md`.
The same one-power source cost remains, but it is carried by the exact
tangential correction `beta u`, not by division by a small principal Kelvin
symbol.

## 2. Full-pressure finite-time response

The Leray field (13) has a standard nonstationary-phase expansion. Since
`k.a=0`, `div(a f)` has no order-`N` term, and solving the whole-space
pressure equation gives

    P(a f)=a f+N^-1 b_1 exp(iNI)+... .                       (14)

The arbitrary finite-order full-pressure recurrence in 0112 applies to the
transported phase `I` and to the transported angular factor
`g(theta-Omega(I)t)`. Its initial coefficients are fixed by the asymptotic
expansion of the actual projection in (13), rather than chosen freely. The
exact generator (5) replaces only 0112's regular `k.omega!=0` Kelvin
inversion. Evolving the exact linear Euler datum (13) and reconstructing
`xi` from

    xi_t=-u.grad xi+Du xi+v                                (15)

therefore gives, for every fixed `T`, derivative inventory and order `M`, a
full-pressure expansion with error `C_(M,T) N^-M`. The displacement
constants additionally contain the polynomial costs in (10). Increasing
`M` before choosing the sparse carrier makes these errors smaller than the
fixed-tag coefficient and than the finite synthesis cost in 0254. No
spectral stability theorem or unit-Floquet-multiplier inference is used.

This all-order conclusion is literal on the exact 0211 action-angle annulus.
It does not transfer automatically to the 0145/0147 periodic approximant:
there one still needs a finite-order invariant-action normal form whose
eikonal defect, pressure hierarchy and cohomological inverse obey the same
joint ordering.

## 3. Exact initial material-observation rows

Let `chi=chi(I)` be one stationary positive material tag. For any scalar
laboratory weight `F(x)`, its material moment has first variation

    delta M_F=rho integral chi xi.grad F dx
             =-rho integral F xi.grad chi dx
             =(rho/lambda) integral J chi'(I) f F dI dtheta. (16)

The second equality uses compact support and `div xi=0`; the last uses (7).
Thus the large tangential correction in (5) cancels exactly from every
scalar configuration moment. Formula (16) is an actual observation map,
not a passive-scalar analogy.

In particular, with `F=x_j` it is the tagged centroid numerator, and with
`F=x_i x_j` it is the complete symmetric shape/covariance row:

    delta X_j=(rho/(lambda M_tag))
                 integral J chi' f x_j dI dtheta,
    delta I_shape,ij=(rho/lambda)
                 integral J chi' f x_i x_j dI dtheta.       (17)

If `Q_0` is the nondegenerate reference covariance and `L_angle` is the
fixed differential of its Gram-angle chart, the exact initial optical row is

    Gamma_angle[N,A,g]
       =L_angle[(rho/lambda) integral J chi' f (x tensor x)dI dtheta].
                                                                    (18)

It can vanish: nonzero polarization does not imply a nonzero Fourier
coefficient of `x tensor x`, or nonzero image after `L_angle`.

The axisymmetric geometry nevertheless computes two useful nonzero rows.
In global cylindrical coordinates the complex centroid weight is
`x+i y=r exp(i varphi)`. A source with toroidal number `n=-1` and a radial
envelope supported where the resonant part of `chi'` has fixed sign gives a
nonzero transverse centroid coefficient. The 0217 physical tilt uses

    F_tilt=(x+i y)z,       theta_+=delta M_Ftilt/(i D_tag),
    D_tag=I_xx-I_zz>0.                                     (19)

In a thin ring `F_tilt=-R s sin(theta) exp(i varphi)+O(a^2)` in the 0217
orientation. Choosing `n=-1` and either matching first poloidal harmonic
makes (16) nonzero at leading order. Hence it remains nonzero on one
sufficiently large but fixed exact ring. This proves a transverse physical
covariance-angle gain, rather than merely naming its Fourier coefficient.

It also exposes an exact missing row. An action-only stationary tag is
axisymmetric, so `I_xx=I_yy`; rotation about the ring's symmetry axis has no
linear covariance Gram angle for that one tag. The construction above
therefore does not by itself supply the longitudinal component of the
unrestricted optical axial vector. A covariant whole-law assembly of rotated
transverse charts, or an actual additional stationary nonaxisymmetric
material observation, must be constructed before claiming that component.
The former belongs to the still-separate geometry/join scope and is not
inferred here.

One frozen positive smooth tag can supply all radial dual phases. For a
sparse sequence `N_j`, choose

    chi(I)=chi_0(I)[1+epsilon sum_j c_j cos(N_j I)],
    c_j=N_j^-j,                                           (20)

with `epsilon` small. This is positive and `C-infinity`. Pairing the
`exp(iN_j I)` source with the conjugate term in `chi'` makes the resonant
centroid/covariance gain scale as `N_j c_j` times the relevant angular
Fourier coefficient. Smooth nonresonant terms are superalgebraically small.
Localizing `A` near `I_0` makes the leading history arbitrarily close on a
fixed time window to `exp[-i nu(I_0)t]`; the localization norm is finite and
is chosen before the carrier. Varying `I_0` over the subannulus realizes
the open band (12) without changing the background or tag.

The antisymmetric displacement and mechanical-current rows do see the
cohomological correction. Their exact initial values are

    G_N=rho integral chi x cross (xi_perp+beta u) dx,       (21)
    S_N=rho integral chi {
          x cross[-u.grad xi+Du xi+v_N]+2 xi cross u} dx,   (22)

where `v_N=P(a f)` and `beta` is the solution (9). These are precisely the
repository's material definitions; in particular

    (G_N)_t-S_N=-2rho integral chi xi cross u dx.           (23)

For a pure `(m_0,n)` harmonic, the leading part of (6) is
`h=iN f/lambda+O(1)`, so

    beta=N f/(lambda nu)+O(1),
    G_N=(rho N/(lambda nu)) integral chi f x cross u dx+O(1),
    S_N=-i nu G_N+O(1).                                    (24)

The last relation uses (23) and
`xi cross u=xi_perp cross u=a f/lambda`, which is one order smaller than
the displayed tangential contribution. After pairing with the fixed tag,
the resonant leading gains in (24) scale as `N_j c_j`. The transverse
coefficient is explicit in cylindrical coordinates:

    (x cross u)_+
      =exp(i varphi)[-W z-i V(r cos(theta)+z sin(theta))].    (25)

It is not identically zero on the 0211 core (`V` and `W` are both retained),
so `n=-1` and a matching poloidal component select an individual nonzero
`G/S` column at sufficiently large fixed radius. What remains unproved is
the determinant separating that current column from the angle row at the
same frequency and, after any whole-law assembly, its full three-component
rank.

Equations (18), (21) and (22) are the literal angle, angular-displacement
and full-spin/current gain rows to place in the common source matrix. In
Fourier form each is the action integral of the conjugate `N_j` tag
coefficient times the `(m_0,n)` coefficient of its displayed physical
kernel. They provide a finite computable matrix, but neither symmetry nor
(4) proves that its angle/current determinant is nonzero. Likewise, the
centroid row (17) is the acoustic hybrid displacement contribution; the
point-minus-centroid acceleration requires differentiating the complete
Euler/pressure history and retaining the ambient point row before claiming
the 0243 hybrid coefficient. A nonzero centroid coefficient alone is not
that claim.

## 4. Action row and exact remaining construction

For the real and imaginary quadratures of one complex source (5), both
generators lie pointwise in `span{grad I,u}`. Their cross product is parallel
to `grad I cross u`, and hence

    Omega_KKS(xi_Re,xi_Im)
      =rho integral omega.(xi_Re cross xi_Im) dx=0.         (26)

This is an exact zero, not a small-carrier estimate. It shows that the two
transport quadratures cannot by themselves be advertised as a canonical
phase pair. It does not make the velocity response passive, and it is not a
parent no-go: their pairing with separate source/normalizer columns can be
nonzero, and 0228's full-form construction may be applied only after all
such columns and their cross entries are fixed. The Jacobi energy and every
cross form remain actual entries, not inferred from (23).

`route_verdict: exact fixed-ring isovortical Euler response band established;
the complete physical supplier remains open at measured joint gain and
normalizer closure`

`evidence_scope: exact divergence-free Lin preparation on an open twisted
action annulus of the fixed 0211 constant-curl ring; arbitrary-order
full-pressure finite-time Euler/Lin expansion with polynomial source cost;
exact centroid, covariance-angle, G and spin/current gain formulas; exact
zero mutual KKS row for the two transport quadratures`

The minimum remaining construction is to evaluate the finite determinant
separating (18), (21), (22) and the ambient point/pressure row for enough
source columns; supply the longitudinal covariance component by an actual
allowed common observation; and prove that the resulting acoustic-both-
parity/optical-angle/full-current map has the required range on a nonempty
subband. Then append separate normalizer columns and compute the complete
KKS/Jacobi/current Gram matrix before using the 0254 residual/gain diagonal.
For the 0145/0147 final periodic field, an additional finite-order invariant-
action/cohomology transfer is still required. No curvature radius is varied
in the result above.
