# Isotropic reconstruction from actual transverse covariance tilts

## 1. Per-ring physical observation

Let an axisymmetric tagged ring have unit symmetry axis `n` and reference
centered covariance

    Q_n=q_perp(I-n tensor n)+q_parallel n tensor n,
    D_tag=q_perp-q_parallel>0.                             (1)

The strict gap is the rotated version of the actual 0217 coefficient
`I_xx-I_zz>0`. It makes the symmetry axis a simple covariance eigendirection,
even though the two transverse eigenvalues agree. For an infinitesimal
physical rotation with axial vector `Phi`, covariance transport gives

    delta Q_n=[Phi cross]Q_n-Q_n[Phi cross],
    delta n=Phi cross n.                                  (2)

Equivalently, the off-axis covariance row is

    (I-n tensor n)delta Q_n n=-D_tag delta n,
    delta n=-(I-n tensor n)delta Q_n n/D_tag.              (3)

This is exactly the covariant form of 0217's
`Q=iD_tag(a_x+i a_y)` and its unit physical tilt chart, with the sign fixed
by `D_tag=q_perp-q_parallel`. The
component of `Phi` parallel to `n` is correctly absent from one ring's chart:
it rotates only the degenerate transverse eigenspace.

The new whole-law observation reconstructs after applying this literal
per-realization covariance chart:

    Phi_rec=(3/2) E_n[n cross delta n].                    (4)

The vector triple-product identity gives

    n cross (Phi cross n)=(I-n tensor n)Phi.               (5)

For normalized isotropic axes, `E[n tensor n]=I/3`, and hence

    Phi_rec=(3/2)(I-I/3)Phi=Phi.                           (6)

The gain Gram is exactly `E[I-n tensor n]=2I/3`, with condition number one.
The formula is well defined for an unoriented axis: replacing `n` by `-n`
also replaces `delta n` by `-delta n`, leaving `n cross delta n` unchanged.

This ordering is load bearing. Averaging the reference covariances first
would give a scalar matrix `E[Q_n]` and destroy the simple eigenaxis needed by
(3). The permitted common-law convention does not do that: C-CST-017 rotates
the whole state, source and physical observation together, and the inherited
0234 estimator likewise averages already measured per-realization angle rows
before applying its finite Gram inverse. Equation (4) is the vector-tilt
version of that same observable, not a new claim that one axisymmetric tag
measures axial rotation about itself.

## 2. O(3) covariance and the actual source map

For a proper rotation `R`, rotate the complete reference realization by

    u_R(x)=R u(R^T x),
    xi_R(x)=R xi(R^T x),
    chi_R(x)=chi(R^T x),
    p_R(x)=p(R^T x).                                      (7)

If `n=R e_z`, feed the reference two-column transverse source the local
amplitude `R^T P_n Phi`, where `P_n=I-n tensor n`. The `m=-1` complex pair is
equivariant under rotations about `e_z`, so this definition is independent
of the choice of transverse frame in `R`. Its observed axis change is

    delta n=(P_n Phi) cross n=Phi cross n,                 (8)

as required. For an improper `R`, `Phi` transforms axially as
`Phi_R=det(R)R Phi`; the cross-product rule then makes `delta n` transform as
a polar vector and `n cross delta n` as an axial vector. Thus (4) is an
actual O(3)-covariant axial observation. Inversion and time reversal retain
their separately declared preparation signs and do not change the projector
average (6).

Rigid rotations commute exactly with the whole-space Leray projector:

    P[R f(R^T x)]=R(Pf)(R^T x).                            (9)

They therefore conjugate the full linear Euler and Lin generators, including
pressure, rather than only their principal symbols. Every Euclidean `H^s`
norm and every finite time/parameter derivative norm of (7) equals the
reference norm. The material mass fraction stays positive and has the same
mass. Covariance gaps, centroid/shape moments and mechanical currents rotate
tensorially with unchanged magnitude.

For a source operator `S_n Phi=S_ref[R^T P_n Phi]` rotated by (7), let
`C_s` be any finite reference source/pressure/Lin norm. Then

    sup_n ||S_n Phi||_s <= C_s |Phi|,
    E_n ||S_n Phi||_s^2 <= C_s^2 E|P_n Phi|^2
                         =(2/3)C_s^2|Phi|^2.               (10)

Thus the continuum orientation law introduces no divergent source cost.
The same invariance holds for the KKS form and stationary Euler Jacobi form:

    Omega(Rxi,Reta)=Omega(xi,eta),
    H(Rxi,Reta)=H(xi,eta).                                 (11)

Their common-law averages and all cross entries with other covariantly
rotated finite sources are finite. Equation (11) does not normalize those
forms or change the exact null mutual KKS entry of 0250's two transport
quadratures; the separate full-form normalizer remains necessary.

## 3. Angle and current after averaging

If one local transverse source branch has physical angle gain `c_theta` and
mechanical displacement/current gains `c_G,c_S` on the plane `n_perp`,
axisymmetry permits, in each fixed real phase convention,

    delta n=c_theta (a cross n),
    G_n=c_G P_n a+c_G^J(n cross a),
    S_n=c_S P_n a+c_S^J(n cross a).                        (12)

The possible `J` terms record a chiral ninety-degree plane response. Under
the unoriented isotropic law, `E[n cross]=0`, while `E[P_n]=2I/3`. Therefore

    Phi_rec=c_theta a,
    E[G_n]=(2c_G/3)a,
    E[S_n]=(2c_S/3)a.                                     (13)

After the unit angle normalization `c_theta=1`, (6) supplies all three
components with no extra conditioning. The physical `G/S` rows remain their
actual rotated averages. In particular, whole-law reconstruction does not
separate angle from `G`: one local column with a fixed pair
`(c_theta,c_G)` remains one averaged column with pair
`(c_theta,2c_G/3)`. A second physical source direction or the vector-adjoint
uniqueness mechanism in 0250 is still required.

If a local branch obeys the exact material identity

    S_n=G_n,t+2rho integral chi xi cross u,                 (14)

averaging preserves every term of (14). Only a source combination that makes
the last averaged row vanish earns exact `S=-i nu G` for a monochromatic
phase. The large-carrier limit found in 0250 preserves the asymptotic ratio,
but (4)--(13) do not manufacture exact compatibility.

## 4. Supplier boundary

`route_verdict: established at full-vector covariance-observation and
uniform rotated-cost scope`

`evidence_scope: literal simple-eigenaxis covariance tilt, exact isotropic
projector inverse, O(3)-covariant whole-state source construction, and exact
uniform Euler/Lin/pressure/material/action norm preservation`

The precise supplier license is:

> Under the already permitted normalized isotropic whole-state law, the two
> actual transverse covariance tilts of rotated copies of one fixed-radius
> axisymmetric ring reconstruct a faithful three-component axial vector by
> `(3/2)E[n cross delta n]`, with gain Gram `2I/3` and no orientation-dependent
> growth of any declared finite source, pressure, material or action cost.

This closes the longitudinal covariance-observation component identified in
0250. It does not prove angle/`G` separation, the ambient point-to-hybrid
coefficient, vector-valued time density, the full KKS/Jacobi normalizer,
compact-ring density assembly or the issue200 parent.
