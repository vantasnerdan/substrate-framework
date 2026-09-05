# Full joint translation–angle jet, its normal form, and its physical response

## 1. Exact statement and derivative counting

N3 asks for an objective isotropic affine energy in the displacement gradient
and independent microrotation, with slow fields; N4 also specifies the kinetic
action. This is a second-gradient continuum target, not an arbitrary finite-k
Euler invariant-subspace target. Nevertheless a coordinate transformation
cannot create a physical centroid response absent from the original action.
The following identities distinguish these two statements precisely.

Use Fourier convention exp(i k.x), and a transverse curl helicity h=+/-1,
so curl has eigenvalue h k. The displacement is polar and the physical angle
is axial. Time-reversal pairing removes the odd temporal terms; reflection
pairing removes the chiral spatial terms. Through total spatial derivative
degree two the general real Hermitian transverse matrices are

    M = [[rho+m_U k², b h k], [b h k, j+m_Phi k²]],
    K = [[A k², g h k], [g h k, kappa+C k²]],

with rho,j,kappa positive. Coefficients here are the complete joint-action
coefficients after all fluid reactions have been varied. No origin or sign
of g,b,C is inferred from the display. In the physical objective Cosserat
coordinate, g=-kappa/2 and A=mu+kappa/4, with kappa=4 alpha.

In three dimensions replace scalar transverse coefficients by transverse and
longitudinal projectors. The polar/axial off-diagonal block is b curl in M
and g curl in K; it has no longitudinal part. The remaining longitudinal
entries are rho+m_U,L k², j+m_Phi,L k², A_L k² and kappa+C_L k².
For incompressible Euler the longitudinal displacement is absent, while the
longitudinal spin remains. These exhaust the reflection-even isotropic jets
at this degree. In particular a parity-even constant scalar mixed mass or
stiffness between a polar and an axial vector is unavailable.

## 2. General mass normalization pulls back the potential too

Put d=m_Phi-b²/rho. The map from normal fields to physical fields is

    U_phys=(1-m_U k²/(2rho)) U_normal-(b/rho) h k Phi_normal,
    Phi_phys=(1-d k²/(2j)) Phi_normal.

Exact congruence followed by coefficient extraction through k² gives

    M_normal=diag(rho,j),
    K_normal=[[A k²,g h k],[g h k,kappa+C_eff k²]],
    C_eff=C-2 g b/rho-kappa*(m_Phi-b²/rho)/j.

Thus 0053's special formula is correct when its stipulated g=-2 alpha is
used. A different computed mixed potential g must enter the general formula;
replacing it by -2 alpha without an objective field-map derivation would
change the spin-gradient coefficient. The longitudinal spin coefficient is
C_L,eff=C_L-kappa*m_Phi,L/j. The longitudinal displacement map only normalizes
its mass at this degree, and leaves A_L unchanged.

The map contains curls and second derivatives, so it preserves uniform rigid
translations and common rigid frame rotations. It does change nonuniform
centroid and spin observables. Its inverse and the corresponding boundary
current improvement are part of the result, rather than a license to rename
the new displacement the unchanged mass centroid.

Two invariant root expansions expose incorrect omitted terms:

    omega_ac² = (A-g²/kappa) k²/rho+O(k⁴),
    omega_op² = kappa/j
       +[C/j-kappa*m_Phi/j²+(g-kappa*b/j)²/(rho*kappa)] k²+O(k⁴).

These follow either from the complete determinant or from the same pulled-back
normal matrices. The term m_U does not enter either displayed coefficient;
it first multiplies an acoustic frequency already of order k², or enters the
optical translational response at a higher power. Its absence here is derived,
not an omitted mass term.

## 3. The full relative-angle pullback

Write q=Phi-s curl U, where the physical affine frame has s=1/2. Allow the
most general relative-coordinate jets before specializing:

    M_rel=[[rho+m_0 k²,d_0 h k],[d_0 h k,j+n k²]],
    K_rel=[[a k²,r h k],[r h k,kappa+c k²]].

The exact coordinate matrix T=[[1,0],[-s h k,1]] maps (U,Phi) to (U,q).
Congruence T^T M_rel T and T^T K_rel T gives the full entries

    M_UU=rho+(m_0-2s d_0+s²j)k²+s²n k⁴,
    M_UPhi=(d_0-sj)h k-s n h k³,
    M_PhiPhi=j+n k²,

    K_UU=(a-2s r+s²kappa)k²+s²c k⁴,
    K_UPhi=(r-s kappa)h k-s c h k³,
    K_PhiPhi=kappa+c k².

Thus the curvature of q does generate extra potential terms; their total
spatial degrees are THREE and FOUR, not two. Dropping them in a k² matrix
jet cannot change the displayed acoustic or optical k² dispersion. Conversely
using the truncated matrix to claim the exact finite-k factorization would be
incorrect. Existing 0061 explicitly uses the full congruence for that identity.

For the separable action a=mu, r=d_0=m_0=n=0, s=1/2, the exact determinant is

    det(K-omega² M)
      =(mu k²-rho omega²)(kappa+c k²-j omega²).

Its k² matrix jet is equivalent to the standard Cosserat-shaped normal form
with alpha=kappa/4, A=mu+alpha and C_eff=c-alpha*j/rho. The normal-form optical
slope is C_eff/j+alpha/rho=c/j, just as the exact separable determinant requires.
This is a legitimate formal second-gradient equivalence; exact factorization
does not refute that equivalence. Positive target curvature requires the
computed c and all other contributions to give positive C_eff, rather than
assuming c>0 alone suffices.

## 4. Physical coupling is a different, explicitly computable observable

On the optical branch the ORIGINAL physical centroid displacement satisfies

    U_phys/Phi_phys = (j*g/kappa-b) h k/rho+O(k³).

Define l=g-kappa*b/j. This is the leading physical centroid/optical-angle
mixing coefficient in the stated physical coordinates. In relative variables

    l = r-kappa*d_0/j,

independent of the choice of s in the angle map. For the separable action
l=0: the original centroid has no leading optical response. In fact it is
exactly decoupled when all orders of the exact pullback are retained.

The real-field variation makes the mechanism transparent. If E_U and E_q
are the two Euler–Lagrange residuals in the original (U,q) variables, then

    E_Phi=E_q,
    E_U,physical=E_U-s curl E_q.

Self-adjointness of curl here includes the periodic or explicitly retained
boundary term. Once the spin equation holds, the physical translation
equation is exactly its old independent equation. The apparent asymmetric
stress term and the apparent kinetic cross term cancel in that observable;
discarding the latter would manufacture leading translation–spin transfer.

Accordingly the separable construction DOES supply a conditional Cosserat
normal-form action through second gradient when its full coefficient
inequalities hold. It does NOT by itself supply the additionally requested
physical centroid transfer or a new mode hybridization. The normal displacement
contains a spin-dependent curl filter of the original centroid. Reporting
that filter explicitly preserves the strongest true result without replacing
the physical response requested by the parent.

## 5. What joint term can actually repair the physical response?

An arbitrary order-k potential term r q.curl U is not an available objective
repair: q is a relative angle and stays fixed under a common rigid frame
rotation, whereas curl U shifts. Its cross variation would change the energy.
Moreover an isotropic scalar coupling between an axial vector q and the
symmetric affine strain has coefficient proportional to epsilon_ijk, whose
contraction with that symmetric strain vanishes. Thus in the declared
reflection-even isotropic objective affine ensemble the relative-coordinate
order-k stiffness r vanishes. This conclusion does not assert that the
unaveraged cell mixed Hessian is zero; it identifies the representation
component surviving the stipulated ensemble average.

An actual relative kinetic coupling d_0 between qdot and curl Udot IS allowed:
a superposed time-independent rigid frame change does not shift the velocities.
It must be derived from the same Euler material/Routh action, including the
affine fluid reaction. With r=0, it gives l=-kappa*d_0/j, nonzero if d_0 is.
The exact pullback above repairs every affected mass and stiffness coefficient:
b=d_0-j/2, m_U=m_0-d_0+j/4, A=mu+kappa/4, g=-kappa/2,
m_Phi=n and C=c. No separately assigned rotor inertia is needed to state this
repair target, and no value or sign for d_0 has been postulated.

The eleven-affine KKS projection in 0059 legitimately constructs its selected
zero-pairing sector, but that positive theorem does not prove a nonzero d_0.
A next joint construction can retain and compute an affine/reaction block
instead of projecting it away, then take the FULL block Schur complement.
An independently selected nonzero physical momentum response, or the actual
material moving-frame kinetic connection, can supply that block; merely
assigning a potential r forbidden by objectivity cannot. All finite gradient
corrections must then use the general C_eff and full reaction jet, not the
separable specialization.

## Route verdict

Established: complete joint-jet normal form, induced-term order accounting,
and exact distinction between second-gradient Cosserat equivalence and
physical centroid–spin transfer. No same-degree potential term was missing
from the separable k² pullback. The parent physical-transfer construction
requires a computed nonzero l; the separable 0059 specialization has l=0.
This names the needed joint response without an all-k requirement or a new
scope gate. The continuation is the same-Euler affine/reaction kinetic block,
with the exact coefficient repair above ready for its output.
