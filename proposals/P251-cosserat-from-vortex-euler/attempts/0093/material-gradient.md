# Positive material curvature after its added inertia

## 1. Same-action inputs and exact finite margin

Use the material Jacobi forms of0084/0090, including density:

    M(V,W)=rho integral V.W,
    K(V,W)=integral V.Hess(p)W-rho (u.grad V).(u.grad W).

The background is stationary smooth Euler. These forms are local in the
ACTUAL divergence-free material displacements; they are not the coadjoint
Hessian with a Leray tail. The time-reversed backgrounds have identical M,K
and opposite gyroscopic form. Their constrained material shapes are varied
coherently, with pressure, mean gauge and individual circulation variables
retained as in0091. This proposition adds prescribed gradient attachments
to that finite Cauchy--Born family; it does not freely relax the attachments.

For0091's radial rotation B_L, write

    M_B=m L^5, A_B=a L^5, d=a-m>0, |K_B|<=C L^3,
    m>0, a>m, D=d L^5.

Here A(V)=rho integral (e cross r).V is the complete-fluid spin row.
Explicitly a,m,d are rho times the integrals with factors chi, chi²,
chi(1-chi); C=(3rho/2)||u||_infinity² integral ||D B_1||_F².
All are computed geometric or Euler quantities. Choose a disjoint compact
positive material cage X from0090, with m_X>0, k_X>0, a_X=A(X), and

    t=(a_X+sqrt(a_X²+4m_X d L^5))/(2m_X),
    Xi=B_L+t X,
    j=M(Xi)=m L^5+m_X t²,
    kappa=K(Xi)=K_B+k_X t².

The exact root enforces A(Xi)=M(Xi), not a target frequency. Put r_X=k_X/m_X.
The important cancellation, valid at each finite L, is

    r_X-kappa/j=(r_X m L^5-K_B)/j.                 (1)

No a_X or t survives in the numerator. Thus a cage with positive stiffness
alone is not enough, but THIS same cage eventually has a strictly higher
stiffness/mass ratio than the completed physical-spin mode.

Here is a finite sufficient choice, without a limiting sign test. Require

    |a_X t| <= min(a,d/2) L^5,
    C L^3 < min(r_X m/2, r_X d/2) L^5.

Such L exists because t=O(L^(5/2)); explicitly
t<=|a_X|/m_X+sqrt(d/m_X)L^(5/2), which provides a polynomial bound for the
first inequality. Then kappa>0, j<=2a L^5 and

    r_X-kappa/j >= gamma := r_X m/(4a)>0.         (2)

For translated cages on u_E+epsilon b, choose their positions at sufficiently
far periodic-copy sites after L is selected. 0090's C² continuity estimate
gives a uniform lower bound k_X>=k_*>0 and convergence of every such ratio
to the same positive periodic value. This also gives a finite, non-circular
uniform choice of gamma using k_*/m_X. Additional copies can therefore be
chosen with ratio at least kappa/j+gamma/2. No equality of their stiffnesses
on the nonperiodic EPS component is assumed. Their masses and geometric
spin moments, unlike their stiffnesses, are exactly translation invariant.

## 2. Paired gradient cages and exact mean moments

For each attachment use two translated copies of the SAME compact smooth
divergence-free profile, at disjoint supports away from all existing spin,
strain and core supports:

    Y(y)=X(y-a_plus)-X(y-a_minus).

For every compact divergence-free X, integration by parts gives

    integral X_i=0,
    integral(y_j X_i+y_i X_j)=0.

The first moment is independent of the translation of its support, because
the zeroth moment is zero. Consequently Y has BOTH moments zero:

    integral Y_i=0, integral y_j Y_i=0.           (3)

In particular A(Y)=0, every affine kinetic pairing integral (c+Hy).Y
vanishes, and its linear complete-fluid shape and spin rows vanish.
Locality and disjointness give

    m_Y=2m_X,
    k_Y=k_plus+k_minus,
    k_Y-(kappa/j)m_Y >= (gamma/2)m_Y>0.          (4)

The corresponding microscopic displacement at a marked cell is

    tau [n.q(X_cell+h/2)-n.q(X_cell-h/2)] Y(y).

The bracket is a cell-coordinate amplitude, constant on each microscopic
profile, so material incompressibility is exact; multiplying X by a varying
physical-space envelope and omitting its solenoidal return is not used.
The profile has no core support and therefore leaves the core-angle
observable unchanged. A smooth long-wave field supplies the amplitudes,
and the exact difference symbol is 2i sin(k.h/2).

Equation (3) also treats the actual mean, not just isolated cells. Under the
stationary marked assembly, the coherent Fourier mean of Y is O(|k|²), by
Taylor's formula with its finite second moment. Its gradient amplitude is
O(|k|), so the added mean displacement and velocity are O(|k|³). A bounded
Helmholtz projection in the volume-preserving GLM mean does not lower this
order. Thus mean-centering changes no added coefficient through order k².
This statement concerns the NEW profiles;0091's pre-existing mean,
circulation and core-frame lift still has to be the same full material one.
The spin/covariance current corrections from these zero first moments also
start beyond the retained leading curl-of-spin mean map. Higher moment
currents are retained in the exact discrete observable, not set to zero.
For example, with the Fourier density convention exp(-ik.y), the added
spin row itself need NOT vanish at order k². Its exact profile factor is

    rho tau 2i sin(k.h/2) integral [y cross Y(y)] exp(-ik.y) dy,

whose leading coefficient is

    rho tau (k.h) k_j integral y_j [y cross Y(y)] dy.

This is an explicit O(k²) physical-spin filter, linear in tau; taking its
curl in the momentum improvement makes it O(k³). It is not the quadratic
gradient mass tau² m_g, and the two are not silently identified. The
canonical angular momentum contains the gradient mass; their difference
is the retained higher material moment/current improvement. The actual
spin observable follows this Fourier profile and the field map in section4.

## 3. Full kinetic and stiffness additions

Take three orthogonal bond directions h_l of equal length h, with distinct
disjoint paired supports, and average their complete configurations over
proper rotations and reflection. An independent isotropic axis mark n
has E[n n^T]=I/3. Equal masses are not needed for positivity: a uniform lower
bound on (4) gives a positive definite bond/axis tensor. In the identical
coefficient special case, with marked intensity nu,

    Delta M_PhiPhi = tau² m_g |k|² I,
    Delta K_PhiPhi = tau² k_g |k|² I,
    m_g=nu h² m_Y/3, k_g=nu h² k_Y/3,
    k_g-(kappa/j)m_g >0.                        (5)

For nonidentical actual-EPS copies replace these by their Palm expectations;
(4) survives averaging and gives the same strict tensor lower bound. The
angle q is relative to the macro frame. Replacing q by Phi-curl U/2 makes
the new mixed entries start at degree three and the new U entries at degree
four. Those exact pullback entries exist; they do not contribute to the
degree-two matrix jet.

All local material cross terms with previously selected compact internal
profiles vanish by disjoint support. The bare macro translation has
K(c,Y)=0 by integrating the pressure Hessian against div Y=0. Its affine
potential cross with a gradient angle is a parity-odd rank-four polar/axial
tensor and vanishes in the stated reflection-paired ensemble. Its affine
kinetic cross already vanishes exactly by (3). Thus the NEW degree-two
translation blocks and degree-one mixed blocks are zero, not uncomputed.
The shear-only family0094 contributes its actual U-gradient mass and shear
stiffness separately; disjointness removes its local cross with Y.

No isolated-cell inverse is taken. If additional finite material reaction
coordinates have been retained, choose their frozen profiles before the
new supports. Their cross blocks with Y vanish for the same local-form
reason. An exact algebraic Schur reduction of those blocks consequently
leaves (5) unchanged. This is a statement about the declared constrained
finite family, not an assertion that an unrestricted infinite Euler
reaction space is support-orthogonal. A new freely varying reaction that
overlaps Y requires its actual Schur contribution; it cannot be introduced
and then ignored. Cauchy--Born excludes that nonaffine relaxation. The
mean-pressure gauge contribution, which is not an arbitrary free shape,
is dealt with independently by the O(k³) mean estimate above.

## 4. Complete normal form, not stiffness alone

Let the ALREADY COMPUTED reflection-even transverse action before adding
the cages have complete jets (curl helicity sigma=+/-1)

    M=[[rho+m_U k², b sigma k],[b sigma k,j+m_Phi k²]],
    K=[[A k²,g sigma k],[g sigma k,kappa+C k²]].

This includes the pre-existing mean/current and any legitimate reaction
Schur terms; it does not assign b, g or rho from an unrelated action.
The physical-to-normal field map and complete pulled-back curvature are

    U_phys=(1-m_U k²/(2rho)) U_N-(b/rho) sigma k Phi_N,
    Phi_phys=(1-(m_Phi-b²/rho)k²/(2j)) Phi_N,
    C_N=C-2gb/rho-(kappa/j)(m_Phi-b²/rho).

For longitudinal spin, C_L,N=C_L-(kappa/j)m_Phi,L; incompressible Euler has
no longitudinal displacement. Adding (5) changes BOTH m_Phi and C, so

    Delta C_T,N=Delta C_L,N
        =tau² [k_g-(kappa/j)m_g] >0.             (6)

The leading b, g, j, kappa, rho and shear coefficient are unchanged. Thus
all fixed finite pre-existing normal-form corrections can be overcome:
if delta=k_g-(kappa/j)m_g>0, any

    tau² > max(0,-C_T,N,-C_L,N)/delta

gives strictly positive transverse and longitudinal curvature. In the
nonidentical tensor case use its positive minimum quadratic-form bound.
This is a finite geometric attachment strength chosen by structural
positivity, not by matching a comparator or supplying a modulus.

The added spin-gradient mass remains in the physical map. It is incorrect
to report Phi_phys=Phi_N unchanged or to delete the mass before evaluating
the optical slope. The latter changes by delta/j, whereas acoustic k²
stiffness is unchanged. Likewise0094's U-gradient mass belongs to the map,
although it enters acoustic dispersion only at k⁴. Boundary terms are those
of the same mean/current map and periodic or compact slow-field action;
this proof does not replace action positivity by pointwise density equality.

## Verdict and parent interface

Established: given0091's complete finite material/GLM lift, gradient-only
compact material attachments exist with a finite strict positive gain in
BOTH normal-form curvature eigenvalues, after their actual added inertia,
mean-centering, fixed reaction blocks and physical field map are retained.
The proof gives the mechanism and finite bounds; a positive K cage alone
would not establish it. The original smooth-Euler parent additionally uses
0091's same-action mean/core/Kelvin identification and0094's actual shear.
This result supplies the gradient construction, not an independent proof
of those inputs or permission to replace them by old coadjoint coefficients.
