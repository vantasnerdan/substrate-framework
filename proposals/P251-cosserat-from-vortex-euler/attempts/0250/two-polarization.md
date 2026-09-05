# Actual gain/history use of the 0260 two-polarization lift

The exact vector cohomological construction and its symbolic check are
root-authored in `0260/vector-lift.md`. This 0250 artifact imports that bounded
source at its stated scope and computes the physical gain/history matrix. It
does not claim authorship of 0260 or use its nonzero KKS row as a positive
Jacobi normalization.

## 1. Imported exact commutator inverse

On the exact volume action-angle annulus write

    u=Omega(I).partial_theta,       D=Omega.partial_theta,
    f=A(I)exp(iNI+i ell.theta),       ell=(m,n), n!=0.     (1)

Each angle-coordinate field `b_j=partial_thetaj` commutes with `u`, because
the coefficients `Omega(I)` are angle independent. Thus

    D b=Du b,                                                (2)

and every constant linear combination `b` is an exact neutral leading
vorticity amplitude. With `k=grad I`, put

    a_b=-k cross b/|k|^2,       k cross a_b=b.              (3)

Let `w=curl(a_b f)`. For divergence-free `u,xi`,

    curl(xi cross omega)=lambda curl(xi cross u)
                        =lambda[u,xi].                     (4)

The 0260 source proves that the coordinate components of (4) are

    lambda D xi^I=w^I,
    lambda(D xi^theta-xi^I Omega')=w^theta.                (5)

Consequently the proposed solution is correct:

    xi^I=lambda^-1 D^-1 w^I,
    xi^theta=lambda^-1 D^-1 w^theta
                 +D^-1(xi^I Omega').                      (6)

The fixed nonzero toroidal number and the resonance removal already proved
in `isovortical-channel.md` give a uniform smooth `D^-1` on a smaller open
action annulus. Action-flat support of `A` is preserved by this inverse.

The divergence is not an extra assumption. As checked in 0260, since
`div w=0` and `div u=0`,

    0=div(lambda[u,xi])=lambda D(div xi).                  (7)

`div xi` has the same nonzero toroidal sector, on which `D` is invertible;
hence `div xi=0`. Equation (4) now says that
`xi cross omega-a_b f` is curl free. It is compactly supported in the
action annulus, so it is a gradient on `R^3`, and therefore

    P(xi cross omega)=P(a_b f)                             (8)

with the actual whole-space Leray projector. This is an exact finite-`N`
isovortical initial datum, not a symbol-only lift.

## 2. The angle-null tangential polarization

Let `J(I)` be the coordinate volume Jacobian. The tangential covector
components of (3) are exactly

    (a_b)_(cov,theta)=(J b_2,-J b_1).                     (9)

The curvilinear curl formula therefore gives

    w^I=J^-1[partial_theta1(a_b f)_theta2
                    -partial_theta2(a_b f)_theta1]
       =-i(ell.b)f.                                       (10)

For the two tangent directions

    b_u=Omega,             b_perp=(-n,m),                 (11)

one obtains

    ell.b_u=nu=ell.Omega,
    (xi_u)^I=-f/lambda,
    ell.b_perp=0,
    (xi_perp)^I=0.                                        (12)

Here `xi_perp` names the new tangential polarization, not the normal vector
with that name in `isovortical-channel.md`. The leading tangential curl and
inverse are

    w^theta=iN b f+O(1),
    xi^theta=N b f/(lambda nu)+O(1).                      (13)

Thus both source norms have the already declared one-power loss, with
`||xi||_Cs<=C_s N^(s+1)`. The velocity `P(a_b f)` is order one and the 0112
full-pressure recurrence supplies arbitrary fixed order. All inverse,
pressure and time-derivative constants are finite on the selected band.

For any stationary action tag `chi(I)` and scalar laboratory moment `F`,

    delta M_F=-rho integral F chi'(I) xi^I dx.             (14)

The `b_perp` column therefore has exactly zero initial centroid, covariance
and covariance-angle rows, whereas `b_u` has the nonzero angle row already
computed in 0250. In the finite-time WKB hierarchy the leading
`b_perp` action flux remains zero; exact-evolution leakage enters at lower
carrier order and is `O(1/N)` after its `G` normalization. It is retained in
the common error diagonal rather than called identically zero for all time.

## 3. Nonzero G and the full mechanical-spin row

Choose the optical harmonics `n=-1,m=+1` or `m=-1`. In the 0217 thin-ring
coordinates, `b_perp=partial_theta+m partial_varphi`, and

    (x cross b_perp)_+
      =exp(i varphi)[i s(R cos(theta)+s)-m r z].           (15)

Pairing (15) with `exp[i(m theta-varphi)]` gives a nonzero leading angular
integral `2pi i R s` for either sign of `m`. Hence the fixed radial tag ladder
gives

    G_perp=O(N c_N) !=0                                   (16)

on one sufficiently large fixed exact ring. This is a same-frequency
current column with zero leading configuration angle.

Its mechanical spin must keep the exact material correction. From
`G_t-S=-2rho integral chi xi cross u`, define

    R_b=2rho integral chi xi_b cross u,
    S_b=-i nu G_b+R_b+lower WKB sidebands.                 (17)

For `b_u`, the order-`N` tangential displacement is parallel to `u`, so
`R_u=O(c_N)` while `G_u=O(Nc_N)`. For `b_perp`,

    b_perp cross u=(sW-mrV)e_s,                           (18)

whose transverse harmonic is nonzero for `m=+1` when `RV>Ws`, and for
`m=-1` without that subtraction. Consequently

    R_perp=O(Nc_N) !=0.                                   (19)

Thus the proposed exact compatibility `S_perp=-i nu G_perp` is false in
general: its correction is leading, not a pressure remainder. This is not a
route failure; it supplies an independent physical current ratio.

Use the two fixed positive fractions `f_1!=f_2` from the permitted law. The
three same-frequency columns `C_u(f_1),C_u(f_2),C_perp(f_1)` have physical
rows

    [[theta_u,       theta_u,          0],
     [f_1 G_u,       f_2 G_u,    f_1 G_perp],
     [f_1 S_u,       f_2 S_u,    f_1 S_perp]].             (20)

Their determinant is exactly

    theta_u f_1(f_2-f_1)(G_u S_perp-S_u G_perp)
      =theta_u f_1(f_2-f_1)(G_u R_perp-R_u G_perp).        (21)

The first term on the last line is `O(N^3 c_N^3)` after including
`theta_u=O(Nc_N)`, while the second is one carrier order smaller. Every
displayed leading coefficient is nonzero by the thin-ring integrals above.
Therefore (21) is nonzero for every sufficiently large finite carrier,
uniformly on a smaller action-frequency band. This gives an actual
same-frequency inverse for covariance angle, `G` and full mechanical `S`,
with both real time quadratures and with the exact material-current identity
preserved.

## 4. Same-frequency finite-window history inverse

All entries in (20) are smooth functions of the exact action frequency on a
smaller compact band, and (21) is uniformly separated from zero for every
sufficiently large selected carrier. Let `Gamma_N(omega)` be this actual
three-row gain matrix after including its finite lower-order pressure and
normal-leakage terms. Its inverse exists with a finite declared norm; all
frequency derivatives needed by the finite observation inventory are finite.

Choose the finite exponential approximation of the desired
`(theta,G,S)` histories by the constructive 0254 theorem first. At each of its
finitely many frequencies apply `Gamma_N(omega)^-1`, using both real source
quadratures. Then select the full-pressure WKB order and sparse carrier so
that every normalized Euler/Lin and observation error is smaller than the
target divided by the finite synthesis and matrix-inverse cost. This gives
arbitrary `C^r([-T,T])` control of the three local physical histories on the
same fixed ring. The background and the two positive fractions do not change
with `T,r` or accuracy.

Finally apply 0257's per-realization covariance reconstruction and rotate the
complete source/current rows together. Its Gram `2I/3` has condition number
one and its costs are orientation independent. Thus the same construction
controls the full three-component whole-law optical angle, `G` and mechanical
spin histories. This statement still excludes the coupled hybrid-acoustic row
and every action normalization. `axisymmetric-hybrid.md` separately constructs
an acoustic hybrid band; it does not make this optical block's hybrid
coefficients vanish.

## 5. Cross KKS form

The previous null KKS result concerned the two phase quadratures of the same
`b=u` polarization. The cross pairing with `b_perp` is different. Pair real
columns with the same rapid spatial phase. In coordinate components, only
the normal component `(xi_u)^I=-f/lambda` contributes at leading order;
the determinant identity `ell.b_perp=0` cancels `nu` from (13). One obtains

    Omega_KKS(xi_u,xi_perp)
      =rho N/(2lambda) integral J(I)^2 A(I)^2 dI dtheta
         +O(1),                                           (22)

up to the fixed orientation sign and the nonzero angular normalization. It
is nonzero for sufficiently large finite `N`. All-tangent cross products are
normal and have zero dot product with `omega=lambda u`, so they do not alter
the leading coefficient.

Equation (22) supplies a genuine nonzero cross KKS row, but it is not the
complete desired action. After unit physical-output normalization its size
can grow like `1/(N c_N^2)`; that finite cost and the full Jacobi-energy/cross
matrix must still be passed to the separate 0228-type normalizer. No energy
sign or canonical mass is inferred from nonzero KKS alone.

## 6. Result and remaining construction

`route_verdict: exact same-frequency angle/G/full-spin rank and arbitrary
finite-window optical history control established on the fixed 0211 ring`

`evidence_scope: two exact commutator-lifted isovortical polarizations,
angle-null/nonzero-current separation, full material-spin determinant with
two positive fractions, nonzero cross KKS coefficient and explicit
carrier-dependent costs`

The construction closes 0250's angle/`G`/mechanical-spin separation and
finite-window optical history range on one fixed ring. The completed 0257
whole-law projector reconstructs all three axial-vector components with
uniform rotated costs. The axisymmetric continuation in
`axisymmetric-hybrid.md` supplies the full point-to-hybrid acoustic tensor and
its covariant vector inverse. Its reflected-law calculation also proves that
the even optical-to-hybrid and acoustic-to-axial leakages cancel, while the
surviving odd terms are the explicit 0241/0246 curl/current map.  The acoustic
and optical frequency images overlap on one sufficiently large fixed ring, so
these odd entries perturb an invertible common-band block and the joint
physical gain is closed. The separate normalizer must now match the complete
KKS/Jacobi matrix, including (22), before the 0254 finite-window diagonal is
applied. The final 0145/0147 periodic transfer and the independent
geometry/density obligation remain open.
