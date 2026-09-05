# Same-Euler affine angular reaction and absolute physical-angle inertia

## 1. A changed reaction ensemble, with computed moments

This construction changes the selected microscopic reaction directions; it
does not rename 0059's affine-orthogonal action. Use the same actual smooth
Beltrami good patches, robust core geometry, compact physical angle jet, and
finite-carrier negative-helicity cage. The eleven independent affine moment
functionals and their disjoint KKS-isotropic dual responses are those proved
in 0059. Choose their basis as three translations, five symmetric trace-free
strains, and three rotations

    R_i(y)=e_i cross y,    F_A(xi)=Omega(X_A,xi).

All y are measured about the patch reference center. Let eta_A be the fixed
compact responses, with F_A(eta_B)=delta_AB, supported in eleven disjoint
off-core balls. Their span is KKS-isotropic and disjoint from the raw angle
and cage supports. Put Pi=I-sum eta_A F_A.

For a patch with marked physical axis n, take raw Q0=Q_core+C1 and S0=C2.
The actual Euler integral B=Omega(Q0,S0) is nonzero by the finite-carrier
proof. Define the new reaction pair

    Q=Pi Q0,
    S=Pi S0+B sum_i n_i eta_rot,i.

Every coefficient in this definition is a computed Euler KKS integral or an
inverse response Gram coefficient. There is no assigned inertia, locking
modulus or frequency. The equality of rotational and angle moments is an
explicit selected reaction geometry: the dimensionless ratio is one because
both represent the same physical angle, not because an empirical spectrum
was fitted. It is falsifiable by measuring the selected displacement and
moment maps, independently of the resulting energy coefficients.

Exactly, by the disjoint isotropic response construction,

    Omega(Q,S)=B,
    Omega(R_i,Q)=0,          Omega(R_i,S)=B n_i,
    Omega(T_i,Q)=Omega(T_i,S)=0,
    Omega(E_A,Q)=Omega(E_A,S)=0.

Here T denotes translations and E the five symmetric trace-free affine
strains. The physical core jet of Q is unchanged, and S has zero core-angle
jet. The response corrections cannot feed back into B: every response-response
KKS term vanishes, and every raw-response term vanishes by support separation.
This is the explicit same-action repair of the zero affine rotational
pairing in 0059. Rotational pairing is retained; it is not added as a new
independent global-body coordinate.

## 2. Finite carrier positivity and full spatial return

The old fixed moment corrections are bounded uniformly on the selected good
event. So is the new term B sum n_i eta_rot,i: B has a finite nonzero limit
and a uniform bound, and each inverse Gram and response derivative has a
uniform bound. The principal high-carrier helicity contribution to H on Q,S
is positive of order |k_carrier|. Cross terms with a fixed response are bounded
independently of that carrier: disjoint support kills the leading local
helicity product, while full-Leray velocity products are bounded and the
remaining helicity derivative can be moved to the fixed response by integration
by parts. Thus H retains strict positivity at a sufficiently large FINITE
uniform carrier. This proves positivity for the changed pair; it does not
silently reuse the numerical value of the old Hessian.

For independent stationary patch reaction amplitudes, use the full stochastic
Leray projector and the compact parametrix of 0059 before summing the patches.
Disjoint source supports control their L² amplitude norm; the full projector
is a contraction. The same bounds establish a bounded coercive reaction
operator P on the full reaction Hilbert space. The added fixed responses do
not require isolated-patch kinetic-energy factorization.

The force moment now changes in a precisely known way. For F=S cross omega,
translation and symmetric-tracefree moments give

    integral F=0,
    M_ij=integral y_j F_i=c delta_ij-epsilon_ijm B n_m/(2rho).

Its antisymmetric part is not zero and must be retained. A compact vector
potential A with integral A=B n/(2rho) has precisely this first moment for
curl A. Subtract that COMPACT CURL and a compact gradient of mass -c; the
remainder has zero zeroth and first moments. The parent's 0069 construction
then represents the remainder as a smooth compact double divergence. The
compact curl is fixed by Leray, so all full stationary second-gradient jets
remain defined with its cross terms included. No infrared cutoff or discarded
antisymmetric moment has entered the repair.

## 3. The retained affine pairing belongs to the actual material angle

Let U be the material mean displacement already normalized by the full Euler
mass density rho, and beta=curl U/2 its LOCAL affine rotation. At a patch
center its affine displacement is

    U(X_a)+E_a y+beta_a cross y.

The constant and symmetric-strain symplectic pairings with Q,S vanish by the
computed moments. Its remaining pairing with S is exactly B n.beta_a. Thus
the pullback of the SAME Euler KKS form, in the retained coordinates beta,q,s,
contains

    Omega = B (dq+d(n.beta)) wedge ds,
    L_sympl = B s [qdot+n.betadot].

The physical marked core angle is n.Phi=q+n.beta. Consequently the momentum
couples to the physical ABSOLUTE angle rate n.Phidot. No angular kinetic term
has been appended to the action. The affine field supplies beta as a spatial
derivative of U, and the compact fluid reaction supplies B; both are in the
same Euler pullback.

At finite macroscopic wavelength the orbital translation/momentum is already
part of rho|Udot|²/2 and its full mean material action. There is no independent
uniform global rotation variable B(t), no finite assigned orbital inertia of
all R³, and no cancellation of that infinite-volume rigid-body motion by
relabeling the stationary ensemble. Origin shifts in the local rotation
generator only add a translation, whose pairing is zero exactly.

Time-independent frame rotations transport both the microscopic field and
its marks; objectivity of the Euler energy then makes the local potential
depend on relative q and symmetric strain, not beta alone. This Ward identity
is applied to the full transported affine state, not to a fixed-wall isolated
patch energy. Mixed strain/reaction Hessians can remain before elimination.
They are included in the full matrix; isotropy removes a scalar symmetric
strain/relative axial-angle cross only after the complete reduction/average.

The microscopic angle variables and reaction coordinates need not form an
invariant finite-dimensional Euler subsystem. They define the original
declared affine constrained ensemble, using actual volume-preserving
generator flows and a pullback of the closed Euler KKS form. All pressure,
material-boundary and ambient return terms from the full mean action remain
in the parent assembly. This construction supplies its missing order-one
affine rotational pairing rather than replacing those already retained terms.

## 4. The pairing measures physical angular momentum, including ambient fluid

For any of these compact rearrangements,

    delta u=P(xi cross omega),
    delta omega=curl(xi cross omega),    integral delta omega=0.

The rotational Euler moment is

    delta J_e=Omega(e cross y,xi)
      =-rho/2 integral |y|² e.delta omega.

For a ball centered at the patch and containing the source support, the exact
velocity angular-momentum identity gives

    integral_ball rho y cross delta u
       =-rho/2 integral_ball |y|² delta omega
           +rho/2 R² integral_sphere n cross delta u.

The last integral equals integral_ball delta omega by Stokes, and is ZERO.
Thus the compact angular impulse equals actual velocity angular momentum of
the ambient-inclusive ball, exactly and for every larger radius. It is not
just an arbitrary canonical label for a vortex-center coordinate.

The force-moment decomposition of section 2 also gives decay sufficient for
absolute convergence: after the compact curl/gradient pieces, the double-
divergence remainder has zero moments through first order, so its Leray tail
is O(r^-5). The angular-momentum integral therefore converges absolutely. The
mean momentum is zero and the spin moment is independent of the origin.

For Q the moment is zero; for S it is B n. Thus the fluid momentum response
in the retained sector is B s n, the same quantity that multiplies Phidot in
the pulled-back action. In an arbitrary material partition, rather than a
spherical coherence moment, surface and ambient terms do not separately
vanish: the exact current improvements of 0052/0055 are retained. The special
sphere identity proves the physical content of the full moment; it does not
permit dropping material-face fluxes in a different partition.

## 5. Full reaction-operator elimination and its centroid response

Let q=Phi-beta be the relative-angle vector, E the symmetric affine strain,
and z=(q,E). Retain ALL independent fluid reaction coordinates s in the
stationary reaction Hilbert space. D maps the physical angle vector to the
computed patch KKS couplings. The same Euler quadratic action is

    L_int=<s,D Phidot>-<s,P s>/2-<s,N z>-z.H z/2.

P includes every full-Leray interaction and every new response direction.
N includes angle/reaction AND affine-strain/reaction energy blocks. The mean
material action and any explicitly retained background affine two-form are
added only as their already derived Euler terms. The latter is odd under
time reversal, like the entire KKS form.

In the two time-reversed realizations hold physical U,Phi coherently but vary
s_plus,s_minus independently. Their KKS matrices have opposite sign and their
energy matrices agree. Eliminating them gives

    J=D* P^-1 D>0,
    H_eff=H-N* P^-1 N,
    L_pair=Phidot.J Phidot/2-z.H_eff z/2.

The mixed gyro cancels after this variation. An average of isolated-cell
inverses is not used. Positivity of the angle block follows from the full
positive compact angle/reaction construction; simultaneous positivity of
the retained shear/gradient blocks is checked in the parent's full joint
assembly, using the same Schur matrix rather than reusing old coefficients.

After simultaneous orientation/parity averaging of the full action and its
independent reactions, J=j I and the uniform angle block is kappa I, with
j,kappa positive and alpha=kappa/4. All coefficients are Euler integral/
operator expressions. The physical zero-gradient kinetic and affine locking
terms are therefore

    T=rho |Udot|²/2+j |Phidot|²/2,
    W_lock=kappa |Phi-curl U/2|²/2.

In relative variables the SAME kinetic term is
j|qdot+curl Udot/2|²/2: its genuine kinetic connection is d_0=j/2,
m_0=j/4. Substituting in 0066 gives physical b=d_0-j/2=0 and
m_U=m_0-d_0+j/4=0, while g=-kappa/2. Therefore

    l=g-kappa*b/j=-kappa/2 !=0,
    U_phys/Phi_phys=-j h k/(2rho)+O(k³)

on the optical branch. This is actual leading centroid–spin transfer, absent
from 0059's old sector. It is not produced by deleting a kinetic cross term:
the computed affine/reaction pairing has changed the microscopic action.
Other independently derived gradient-current or affine kinetic contributions
must be retained before applying the general 0066 normal form. They cannot
be silently identified with this displayed leading sector or with the old
separable values.

## Route result

Established as stated: a positive same-Euler microscopic angle/reaction
sector with prescribed physical affine angular moment, absolute-angle
inertia, an ambient-inclusive physical spin moment, and a nonzero leading
physical centroid-transfer coefficient. It is a constructive replacement
for the old fully affine-orthogonal sector. Full joint shear and gradient
coefficients, their common Schur reduction and claim promotion remain the
parent assembly's next work; their positive values are not assigned here.
