# Fixed detector versus full current: an explicit boundary action map

This is the registered representative alternative, not an assertion that
a projected spin equals the complete mechanical moment. The physical
cross-section observer is fixed first: the C016 tag chi(psi), its body
X axis, the actual axial spin S_X and its actual normal covariance angle.
The positive whole-field law pushes forward that complete observer.
The finite axial profile in the full-current calculation is a bookkeeping
cut whose material transport and face currents remain explicit.

## The coefficient is measured from the observer difference

Let I_X=rho L_alpha integral chi(Y²+Z²), in the same per-length and cell
normalization as all other moments. For the acoustic rotational rate,
the body-X current is I_X Omega_dot and the complete three-axis trace
is s(t) from full-affine-spin.md. After adding 0231's actual b return,
the projected and full whole-field spins are, respectively,

    S_int=(I_X+cS0)Omega_dot/3,
    S_full=(s(t)+cS0)Omega_dot/3.

Thus their difference is exactly

    Delta S=f(t)Omega_dot,
    f(t)=[s(t)-I_X]/3.                                (1)

The target optical mass j0 and the b coefficient CANCEL from (1).
There is no choice of f after specifying a desired constitutive answer:
it is the full physical response minus the already fixed cross-section
observer. Its nonconstant part is the derived material-cut shear and
correlations, not an assigned microscopic modulus. The observer itself
is stationary on the wrapped cell; the complete finite-cut accounting
knows the chosen cut's age and transport. Confusing those two observables
was precisely the unlicensed complete-spin extension of an axial detector.

For a common slow acoustic field, Omega_dot=curl V/2 at its first
spatial order. Define q(t)=f(t)/2 and the explicit tensor

    Q_ij=q(t) epsilon_ijk V_k.

Then Delta S=div Q at that order. The actual acoustic U_dot differs
from its initial V only at second spatial order, so replacing V by
U_dot in Q changes div Q first at degree three, on the fixed finite
window. Finite smooth preparation and observation remainder constants
are kept in this statement; the ordinary prescribed diagonal must
dominate them. If only a first-jet error bound is available, this map
has precisely that scope and does not upgrade the remainder by naming
it a second jet. Whole O(3) polar-to-axial parity removes the even
spatial coefficient where the actual smooth jet is licensed.

The D-column averaged static spin is zero: it rotates the reference
spin vector, whose whole-field mean is zero. It does not add an
uncomputed time kernel to (1). Initial displacement-current constants
are retained separately, not determined by integrating a spin rate
without its initial value.

## Complete current transformation at the same actual momentum

The coefficient and its initial rows are explicit. With the definitions
in the companion proof,

    q=[2 tr I0+2 Sigma t²-I_X
          +2(C_chi(t)-C_chi(0))
          +2(t C_d'(t)-C_d(t)+C_d(0))]/6,
    q'=[4 Sigma t+2 C_chi'(t)+2t C_d''(t)]/6,
    q(0)=[2 tr I0-I_X]/6,  q'(0)=0,  q''(0)=Sigma/3.

Multiplying the fixed tag by a positive fraction scales I, Sigma,
C_chi, C_d, I_X and q linearly; the normalized bar psi is unchanged.
Equal-probability fractions use their actual weighted q, not a second
normalization by tag mass. The units of q are those of spin inertia
density, after the same cell-density normalization as S.

For the integrated physical current G=G0+integral S, (2) below implies

    G_full-G_int=Delta G0+div integral_0^t Q(s)ds.      (2a)

If an initial spatial charge Q_G0 has actually been supplied, so that
Delta G0=div Q_G0, the accumulated tensor is
Q_G=Q_G0+integral Q. It is generally NOT q epsilon U:

    integral_0^t q(s) U_dot(s)ds
       =q(t)U(t)-q(0)U(0)-integral_0^t q'(s)U(s)ds.   (2b)

The final integral is a physical cut-current memory, including its
initial value. In particular the independent optical supplier uses
S_full,2=S_int,target,2+q curl U_opt,dot and (2a) for its initial/current
row. For U_opt=-beta curl Phi at leading optical order, its addition
is -q beta curl² Phi_dot. This is the extension of the displayed
current/action representative; the acoustic first-jet calculation
alone does not equate a raw optical cross-section detector with it.

Start with the literal hybrid momentum/current and stress representative
of 0232. Define

    S_int=S_full-div Q,
    mu_int=mu_full-Q_dot,
    J_H,int=J_H,      sigma_int=sigma_full.             (2)

The transported spin flux can be kept in the same total flux convention;
at linear order about the prescribed mean state, (2) simply acts on
the total angular-current flux. No background convective term is
dropped if one works with the exact total current instead of splitting
it into names. Equations (2) give the identity

    partial_t S_int-div mu_int+ax sigma_int
      =partial_t S_full-div mu_full+ax sigma_full.      (3)

Most importantly, the physical translational momentum and its U chart
have NOT been changed. The complete angular four-current transforms by

    L_int=L_full-div Q,
    J_ang,int=J_ang,full+Q_dot,
    L=x cross J_H+S.

It is a spacetime superpotential with components B^{0j}=-Q_ij,
B^{j0}=Q_ij. For any fixed comparison region O,

    integral_O L_full=integral_O L_int+integral_boundary Q n,
    J_ang,full n=J_ang,int n-Q_dot n.                 (4)

These are the exact boundary charge and flux to retain. They include
the transverse finite-cut contribution rather than declare it zero.
On the original periodic/compact-variation domain, the integrated
charge difference vanishes. On a cut domain it generally does not.
The microscopic force/torque and moving-face formulae in 0232 compute
the actual full side of (4); arbitrary mechanical boundary tractions
are not equated by convention.

## Same-action realization and unchanged translation momentum

The current transformation has an explicit variational implementation.
For arbitrary smooth physical U,Phi and q(t), set
Q_ij=q epsilon_ijk U_dot,k and add to a local quadratic action density

    Delta L=-(div Q).Phi_dot+Q_dot:grad Phi
           =partial_t[Q:grad Phi]
                       -div[Q^T Phi_dot].             (5)

Thus the full action changes only by its displayed temporal endpoint
and spatial boundary functional. Its Euler-Lagrange equations and
periodic symplectic two-form are unchanged. Its rotational momentum
and couple current shift by -div Q and -Q_dot, respectively, precisely
as in (2). The term Q_dot contains U_ddot; it is part of an exact
boundary identity, not a new higher-time-derivative bulk dynamics.

It is essential to compute the COMPLETE translation momentum of (5),
including its U_ddot and grad U_dot terms. Let
a_ijk=q epsilon_ijk. In components,

    Delta L=-a_ijk partial_j U_dot,k Phi_dot,i
           +(a_dot,ijk U_dot,k+a_ijk U_ddot,k)partial_j Phi_i.

The translation charge from this boundary term is

    partial Delta L/partial U_dot,k
      -partial_t(partial Delta L/partial U_ddot,k)
      -partial_j(partial Delta L/partial(partial_j U_dot,k))=0. (6)

The a_dot pieces cancel, and the remaining opposite derivatives of
Phi_dot cancel. Therefore this is not the earlier momentum-normalizing
field map that altered physical U. It realizes the angular improvement
at the SAME actual translational momentum. Uniform translations and
the measured core-angle observable are untouched.

## Where the memory goes, and the exact remaining interface

In a normalized 0230 local action, choosing the fixed cross-section
internal-spin representative puts the explicit finite-cut current
difference into the boundary map (4)--(5). Conversely the complete
mechanical representative has spin S_int+div Q and couple current
mu_int+Q_dot. Its derived time-dependent cut memory is present there;
it has not been turned into a stationary bulk modulus or deleted from
a physical torque experiment. The autonomous periodic bulk equations
are unchanged by (5).

This proves a legitimate current representative and its explicit full
action transformation. It does NOT, from conservation alone, assert
that the literal pressure-bond mu_full equals a proposed variational
mu_int+Q_dot pointwise: comparison of those actual fluxes, possible
remaining divergence-free localization terms, and the actual supplied
physical histories still belongs to the full constitutive supplier.
The normalized action cannot replace that supplier by a balance identity.
Nor is this a claim of a stationary ensemble of finite-cut material
parcels whose axial variance (5 in the other proof) grows with time.
Stationarity belongs to the fixed wrapped cross-section observer; the
finite cuts are accounted for by their stated boundary/face currents.

The positive result is that the original projected internal-spin observer
can be used without changing actual U or denying the full angular
momentum: its difference is a derived explicit boundary current at
the retained acoustic spatial order. The literal complete-spin route
still requires the exact kernel repair if no such current improvement
is allowed. These are two distinct, fully recorded options under the
original periodic action scope, not an unannounced change of observable.
