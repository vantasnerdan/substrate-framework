# Actual compact reaction columns with exactly vanishing acoustic cross rows

## 1. Actual full-phase cross moments

Use the full constrained material momentum of 0154. For a Kelvin
generator xi its initial cotangent column is
pi_xi=−rho P[(D xi)^T u]. For the common macro velocity V its initial
column is (0,rho V). Integration by parts gives

    Ω(V,xi)=−rho〈V·xi〉,
    Ω(D,xi)=rho〈omega·(D×xi)〉
           =rho〈D·P(xi×omega)〉.                         (1)

D is the actual solenoidal macro displacement, not an optical phase.
The products use the real/conjugate Bloch pairing and the SAME full
fluid volume normalization. Thus the cross rows are determined by the
complete optical displacement and Euler-velocity means m,n, respectively.
Tagged spin cannot be substituted for either mean.

## 2. Fifteen constructive potential moments

Use the compact smooth vector-potential version of 0147's finite
toroidal preparation, supported inside a fixed sub-tube. The truncation
is applied to the potential, as in that source. For a macro Bloch
parameter K define the actual initial generator

    xi_K=curl_K A=curl A+i K×A.

Its physical field is exp(iK·x)xi_K and its initial Euler velocity is
the full P_K(xi_K×omega). The potential is fixed as K varies. Periodic
integration by parts gives the EXACT mean rows

    m(K)=i K×〈A〉,
    n(K)=P_Kslow[d+i B K],
    d_j=〈A·partial_j omega〉,
    B_ij=〈A_i omega_j〉.                               (2)

For example (K×A)×omega=A(K·omega)−K(A·omega); the last term is
removed by the actual slow pressure, not by omitting pressure.
Consequently the fifteen linear conditions

    ∫A_i=0,  ∫A_i omega_j=0,  ∫A·partial_j omega=0       (3)

make BOTH m(K) and n(K) vanish exactly for every K. They suffice
uniformly as identities, without differentiating a large-cell inverse.
Some backgrounds have dependencies between these conditions; the
selected multi-CK geometry supplies an independent control family.

## 3. Rank on an exterior region is proved from the actual field

Choose a fixed open control region disjoint from the full support of
the optical displacement, not merely disjoint from the material tag.
Let chi be smooth, nonnegative, compactly supported there and positive
on an open ball. For the actual multi-CK target of 0147 define the
fifteen vector test fields

    F={e_i, e_i omega_j, partial_j omega}.

These are independent on that ball. Indeed a relation would be

    a0+B omega+(a·∇)omega=0.                            (4)

The entire analyticity of the finite CK sum extends (4) to all space.
Its Fourier support has no zero frequency, so a0=0. Each nonzero
Fourier polarization at g obeys

    [B+i(a·g)I] omega_g=0.

The actual CK Fourier support contains circular cones with two distinct
nonzero positive axial wave numbers k1,k2 and their negative partners.
The coefficients on these cones are nonzero in 0147's construction.
If the horizontal component of a were nonzero, a·g would run over
a continuum on one cone, giving continuously many roots of the
degree-three characteristic polynomial of B. This is impossible.
If a is nonzero vertical, the four distinct projections ±a_z k1,
±a_z k2 give four distinct eigenvalues of a three-dimensional matrix,
also impossible. Hence a=0. The helical polarizations around a
nondegenerate cone span C³, so B=0. This proves the independence.

The Gram matrix G_lm=∫chi F_l·F_m is therefore positive definite.
The local C^r periodic approximation of this fixed field preserves
its strictly positive finite minimum eigenvalue by continuity. No
numerical eigenvalue supplies its sign. All control supports, the
reference Gram margin and norm constants are fixed before this later
approximation accuracy is selected.

Set A_l=chi F_l for the actual approximating background. Their
fifteen moments are precisely G. For each optical potential A_a,
solve the finite real system

    c_a=−G^(-1) moments(A_a),
    A_a,repaired=A_a+sum_l c_(a,l) A_l.                    (5)

The repaired generators are smooth compact curls and (3) holds
exactly. This is an actual exterior Euler initial-data construction,
not an assigned phase matrix. On the original axisymmetric reference,
the large angular harmonic band excludes the low angular orders in
F; its corresponding moments vanish. The axisymmetric transverse
cutoff preserves this angular exclusion. Consequently the coefficients
in (5) are small under the controlled periodic/field approximation.

## 4. Physical tag and action changes have their own scale

The initial displacement on the tag is unchanged EXACTLY, because the
added potential and its curl have disjoint support. Thus its initial
central angle, centroid and displacement moments are unchanged. The
induced velocity includes a genuine pressure tail on the tag; its
spin and angle-rate change is NOT set to zero by exterior support.

If the original fifteen moment residuals have norm epsilon, the
correction has H^r norm at most C_control epsilon. The full Euler
finite-time energy and material-observation estimates bound its
velocity, pressure, angle and spin effects by the corresponding
finite constants times epsilon. Constants include the actual small
quadrupole denominator, packet normalization and observation derivatives.
Choose epsilon below these already fixed relative margins. Existing
spin/angle matching then persists with an explicit error, or the
nearby marker moment IFT can retune its exact finite matching conditions.

The KKS cross between the original packet and an exterior return is
zero by disjoint generator support. Thus the initial optical KKS
changes only quadratically in the return amplitude. With Vcell the
actual cell volume,

    |delta beta_density|≤C epsilon²/Vcell,
    beta_density=beta_packet/Vcell>0.                    (6)

The relative error is C epsilon²/beta_packet, NOT an absolute error
divided by an unrelated bulk energy. Tagged spin-density and angle
errors obey the same per-packet-to-density cancellation of Vcell.
Action/velocity cross terms at later times need not vanish and are
bounded by the actual linear Euler comparison, not by support alone.

For the selected finite cell, sufficiently small K preserves the
nondegenerate optical phase form. Its complete real Bloch encoding is
retained; it need not be a scalar beta(K)J before polarization reduction.
The acoustic/optical cross form, however, is exactly zero by (1)-(3)
and remains zero for all time under the actual joint Euler evolution.
This supplies a genuine block-orthogonal INITIAL phase family. It does
not set the later physical observation cross rows to zero, and does
not assert 0155's unconstructed fixed-cell optical spatial curvature.
