# Finite-core shear from the same actual EPS Euler action

## 1. Physical affine jets and the complete Hessian

Use exactly the smooth stationary EPS Beltrami field in 0045/0048, with
`curl u0=lambda u0`, `omega0=lambda u0`, and the invariant material tube D
of 0052. Fix its mass centroid X_D and write y=x-X_D. The finite coherence
region contains D and the compact returns prescribed below. No periodic
comparison field or preassigned shear modulus enters this construction.

For a symmetric tracefree matrix E define

    A_E=-y cross (E y)/3,
    R_E=curl(chi_core A_E).

The fixed smooth cutoff is one near the physical core points and supported
compactly inside D. Since `curl[y cross (E y)]=-3 E y` when tr E=0,
R_E equals the physical affine displacement E y near those points, with
gradient E. In the transition region its exact extra term is
`grad chi_core cross A_E`; it is retained. Its flow preserves volume
exactly, not just to first order. At the tracked cores its flow derivative
is exp(tE) for sufficiently small t, with determinant one. The finite
core is deformed, not replaced by a line singularity or a director label.

For every compact generator xi the exact transported vorticity satisfies

    delta omega=curl(xi cross omega0),
    delta² omega=curl[xi cross curl(xi cross omega0)].

The full relative Euler energy Hessian, including the second term, is

    H(xi,zeta)=rho integral_R3 [v_xi.v_zeta
                                   -v_xi.curl v_zeta/lambda],
    v_xi=P_R3(xi cross omega0),
    Omega(xi,zeta)=rho integral_R3 omega0.(xi cross zeta).

This is the same functional proved in 0045/0050. For a bounded region U
containing the supports of the transported vorticity changes, it is exactly
the displayed integrand integrated over U PLUS
`rho integral_R3\U v_xi.v_zeta`. Thus the entire exterior velocity Gram
and its finite decaying tail are included. No local inverse, omitted
pressure reaction, or norm-only first variation is used.

## 2. Executed moment projection, without a rank assumption

Apply the fixed-response construction of 0052 to the four physical
functionals consisting of cell mean momentum and common angular moment.
To match units, choose the already specified coherence length L_ref and use

    F_i(xi)=integral_D v_xi,i=integral xi.f_i, i=1,2,3,
    f_i=omega0 cross P_R3(1_D e_i),
    F_0(xi)=Omega(K,xi)/(rho L_ref), f_0=-vK/L_ref.

All f_a are analytic inside D: the Leray responses are harmonic and
curl-free there, and the actual EPS field and its rotation derivative are
analytic. Pick one nonnegative smooth cutoff chi_resp, positive on an open
ball away from all physical core jets, and put

    eta_a=curl(chi_resp curl f_a),
    G_ab=integral chi_resp curl f_a.curl f_b,
    Pi xi=xi-sum_a eta_a (G^+ F(xi))_a.

The selected admissible generators are curls of smooth vector potentials
compactly supported INSIDE D, as are every core, cage and response field
in this construction. Analytic continuation proves that ker G is exactly
the linear combinations of the four functionals vanishing on this space.
Indeed, a zero Gram norm makes curl(sum c_a f_a) vanish on an open set,
hence throughout D; pairing with xi=curl A then vanishes by compact curl
integration by parts. The converse follows by testing compact potentials.
Therefore F(xi) lies in range G and
`F(Pi xi)=0` exactly, even if G is singular. No rank-four assumption or
positive common residual is required to REMOVE these moments; retaining
a common conjugate is the separate 0052 construction.

This precision matters on a solid torus: a curl-free response can have a
nonzero period and need not annihilate every compact solenoidal field with
nontrivial toroidal flux. No such unproved extension is used. The selected
fields and all their projection corrections already possess the required
compact potentials, so this changes none of their physical jets or rank.

All eta_a vanish near the physical jets. Consequently Pi preserves their
exact affine deformation. Every compact divergence-free generator inside
D also has integral_D xi=0, so both the centroid displacement and, after
Pi, its momentum variation vanish. This is the material slice of 0052,
not an energy subtraction followed by an unchanged symplectic form.
Recompute H and Omega on Pi xi. The common K pairing is zero exactly.
No further symplectic orthogonality is presumed.

## 3. A full-rank affine bond construction on this same field

Choose nine directions

    n in {e1,e2,e3,(ei+ej)/sqrt(2),(ei-ej)/sqrt(2): i<j}.

For tracefree symmetric E their actual fractional affine bond extensions
are `b_n(E)=n.E n`. The exact frame identity is

    sum_n b_n(E)²=(3/2) sum_i E_ii²+2 sum_i<j E_ij²,
    ||E||_F² <= sum_n b_n(E)² <= (3/2)||E||_F².

Thus these specified geometric measurements resolve all five shear
components. There is no fitted target tensor or assumed rank.

Place nine small disjoint cage balls inside D, away from the core jets and
the fixed response support, in a region where omega_z has a positive lower
bound. Smoothness and nonzero core vorticity supply such a region. In each
use the actual-field negative-helicity generator of 0045,

    zeta_n=-curl(phi_n p_n)/k_n,
    p_n=(cos(k_n z),sin(k_n z),0), k_n/lambda>0.

The carrier amplitude phi_n has displacement units. Its physical coefficient
is tied to the fractional extension of its specified reference bond. For a
slow displacement u(X), that extension is
`n.[u(X+d_n n/2)-u(X-d_n n/2)]/d_n=n.E n+O(d_n² grad³u)`.
An antisymmetric displacement gradient makes every bond extension zero.

The exact constrained affine family is therefore

    S_E=Pi[R_E+sum_n b_n(E) zeta_n].

This is a kinematic premise, not an elastic energy inserted into Euler:
the designated cage displacement follows its reference bond extension,
with independent conjugate fluid variables subsequently varied in the
same action. Its stiffness is calculated using the complete H above.
Changing the bond/cutoff geometry changes the model and its predictions;
the premise is independently falsifiable. At E=0 there is no added
deformation, and the original common/core-angle action is unchanged.

## 4. Finite positive matrix, including projection and relaxation crosses

Let {E_a}, a=1,...,5, be Frobenius-orthonormal symmetric tracefree matrices,
and set S_a=S_Ea. The finite actual-field matrix is H_ab=H(S_a,S_b).
The cage estimate from 0045 gives

    H(zeta_n,zeta_n)
       >=rho[(1+|k_n|/|lambda|) A_n-C_n],
    A_n=integral (phi_n omega_z)²>0.

Disjoint cage supports eliminate their local helicity cross, not their
whole-space projected kinetic cross. The latter is retained and bounded by
their uniform L² velocity norms. The complete cage matrix therefore has
the positive leading diagonal and a carrier-independent bounded remainder.
The four moments F(zeta_n) are uniformly bounded, and in fact decay as
O(|k_n|^-2) by integration by parts against the fixed smooth f_a. Hence
Pi adds bounded multiples of FIXED eta_a, with no loss of the leading
positive term. Moving curl to a fixed core/response field bounds its
energy cross with any cage independently of carrier frequency. No
high-frequency projection or zero eigenvalue is assumed.

For clarity, an explicit sufficient bound can be assembled as follows.
Use one common sufficiently large carrier magnitude kappa. Let B be the
nine-by-five bond matrix and A_min=min A_n. Its norm satisfies ||B||²<=3/2.
Let C_c bound in operator norm the full cage-matrix remainder divided by
rho, using the 0045 projection-error bounds and the pairwise kinetic norms.
Let V_n bound ||v_zeta_n|| independently of kappa and set
`M_c=||B|| (sum V_n²)^(1/2)`. Write
`S_a=sum_n B_na zeta_n+R_a`; the R_a are bounded combinations of the
fixed core and response generators. Choose uniform

    R0²>=sum_a ||v_Ra||²,
    R1²>=sum_a ||curl v_Ra||²,
    C_s=||B||² C_c+R0²+R0 R1/|lambda|
                            +2(R0+R1/|lambda|) M_c.

Then the COMPLETE shear Hessian obeys

    H_shear >=rho[(1+kappa/|lambda|) A_min-C_s] I_5.

Let p be all selected compact conjugate variables of the joint action,
with its already computed positive momentum block P, and retain the exact
cross matrix N_pa=H(p_p,S_a). For fixed conjugate generators put

    F_p²=sum_p (||v_p||+||curl v_p||/|lambda|)²,
    N_star=rho F_p(M_c+R0).

Then ||N||<=N_star independently of kappa. Their complete static elimination
gives the shear stiffness

    C=H_shear-N^T P^-1 N,
    C >= {rho[(1+kappa/|lambda|) A_min-C_s]
                              -N_star²/lambda_min(P)} I_5.

A finite kappa making the braces positive exists and is explicitly
determined by the displayed SAME-field integrals/norms. One can include all
already positive compact angle variables in the base block instead when
proving positivity of the entire joint elastic matrix. Every projected
KKS entry is recomputed; any resulting strain-rate gradient inertia and
gyroscopic connection are retained in the full Schur action. Time reversal
is averaged only after independently eliminating each realization's
conjugates, as in 0048. This proof neither supplies nor duplicates a mass.

## 5. Isotropic shear coefficient and spin separation

Rotate the complete microscopic realization, bond geometry, physical
core jets, projected generators and its Euler field together by SO(3).
For a declared isotropic distribution with cell number density n_cell,
the symmetric tracefree strain representation has dimension five and is
irreducible. Its exact Haar average is

    average_R C_R=(tr C/5) I_STF,
    W_shear=(n_cell tr C/10) ||E||_F²=mu ||E||_F²,
    mu=n_cell tr C/10>0.

Equivalently n_cell=1/V_cell for a specified single-cell normalization.
The density-weighted sum applies to a distribution of cell geometries.
This is the usual incompressible shear coefficient; no bulk modulus or
compressional Euler mode is inferred from the five tracefree directions.

A local spin angle transforms in the vector (angular momentum one)
representation, whereas E is the angular momentum two representation.
There is no invariant linear map between them. In Cartesian form, an
isotropic rank-three cross tensor is proportional to epsilon_ijk, whose
contraction with E_ij is zero. Hence the ensemble has no local
spin/symmetric-strain potential cross. Spatial reflection pairing also
removes parity-odd strain/wryness cross terms, where present, after the
full action is reduced; it does not erase legitimate even gradient inertia.

## Route result

Established as stated: this same actual smooth EPS field supports a
finite-core, exact constrained affine shear action with a positive derived
isotropic shear scalar, at a finite analytic carrier threshold. The
physical returns, bond attachments, projection, complete Hessian, exterior
motion and static momentum Schur complement are explicit. The parent still
owns construction of its full joint centered momentum sector, kinetic
normal form, isotropic continuum joining and individual claim review.
No unrestricted Euler invariant-manifold claim is used or needed here.
