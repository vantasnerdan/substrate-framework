# Bounded same-curl insertion: prescribed tubes and a positive bulk response

## 1. Prescribed local objects and the exact Fourier measure

Fix lambda>0 and an elementary periodic positive-curl wave u0. Let uT
be the actual entire constant-curl target built in0145, or0147's
multi-CK extension, with any finite periodic EPS knot insertion already
licensed by0145. Fix the finite collection of compact tube neighborhoods
and a larger comparison ball. All desired tube/twist/packet margins are
fixed before approximation. The target is not inferred from a scalar
eigenvalue or a finite numerical field sample.

These explicit targets have a finite vector Fourier measure on the sphere
|q|=lambda: CK terms have a bounded cone measure; finite periodic Beltrami
fields and u0 have finitely many atoms. Thus uT-u0 has the representation

    v(x)=integral_(S2) exp(i lambda n.x) dmu(n),
    n.dmu=0, i n cross dmu=dmu,
    dmu(-n)=conjugate(dmu(n)).                                (1)

Its total variation M is finite. This representation is a property of
the constructed target, not a universal assertion about every entire
Beltrami field.

Let Pplus(n)=[I-nn^T+i cross(n)]/2. It is the Hermitian orthogonal
projector onto the positive-helicity line. Smooth mu with an antipodally
symmetric probability kernel K_epsilon(n,m) supported at angular distance
epsilon, then project:

    a_epsilon(n)=integral K_epsilon(n,m) Pplus(n) dmu(m).

The field h_epsilon with this smooth angular density is real and has
exact curl lambda. Its angular L1 norm is at most M. For every fixed
spatial derivative order r,

    ||h_epsilon||_C^r <= C_r(lambda) M                       (2)

globally, independently of epsilon. On a ball of radius R its difference
from v in C^r is bounded by C_r(lambda)(1+lambda R)M epsilon:
the exponential, the frequency monomial and Pplus are Lipschitz on
the sphere. No exponentially large global Runge norm occurs.

The smooth-angular field also decays as |x|^-1 with each fixed derivative.
One elementary proof aligns the polar axis with x, integrates first over
the azimuth, and integrates exp(i lambda |x|t) against the resulting
smooth function of t in[-1,1]. One integration by parts, retaining both
endpoint values, gives the stated bound. Its constant may depend on the
fixed angular smoothing; bound(2) does not. Thus u0+h_epsilon approaches
u0 spatially while approximating uT on the prescribed ball. Both fields
are actual stationary Euler because their common curl is lambda and
p=-rho|u|^2/2, including every cross term.

## 2. Periodic quadrature with small normalized bulk defect

Partition the sphere into small antipodally paired patches. In each
positive patch choose a distinct rational unit vector n_j, avoiding the
finite u0 directions, and use the coefficient

    a_j=Pplus(n_j) integral_patch a_epsilon(n) dOmega.

Use its conjugate at -n_j. Rational sphere points are dense by the
stereographic map already checked in0145. Every finite selection has a
common denominator D, hence period P=2pi D/lambda. The exact finite sum
h_N is real, positive-curl and periodic. Its coefficients obey

    sum_j |a_j| <= M,
    max_j |a_j| ->0,
    sum_j |a_j|^2 <= M max_j |a_j| ->0.                      (3)

Here j includes both signs. The second statement is absolute continuity
of the now smooth measure; it would be false for unsmoothed atoms.
The same local C^r approximation estimate follows from patch diameters.
The global C^r bound in(2) persists for every N.

Orthogonality on the full periodic cell gives the exact physical energy
identity (the u0 frequencies were excluded deliberately):

    <|u0+h_N|^2>=<|u0|^2>+sum_j |a_j|^2,
    ||h_N||_(H^r,average)^2
                         =(1+lambda^2)^r sum_j |a_j|^2.      (4)

The displayed H^r convention uses the Fourier multiplier
(1+|q|^2)^(r/2). Therefore local insertion strength, bounded global
derivatives, and arbitrarily small normalized bulk energy can coexist.
They are different norms, not contradictory estimates. For a selected
finite N the inserted tube regions have positive, finite spatial density.
Uniform phase and whole-field SO(3) rotation give a stationary isotropic
positive-energy Euler law. It is not asserted Gaussian or translation
ergodic, and is not substituted for an accepted Gaussian-law hypothesis.

Choosing smoothing first and quadrature second preserves every prescribed
local tube margin. This produces ONE stationary field with the elliptic
optical region and the licensed EPS knot regions. Coexistence still does
not identify an optical degree of freedom with a different distant knot.

## 3. Actual bulk response at a fixed finite wave number and time

Freeze a finite set of real wave vectors K and finite time T. Use periodic
or quasiperiodic Bloch data with their complete Euler pressure projector;
the mean is the actual normalized Fourier coefficient. Define actual
common-V and Kelvin-D preparations from u_N=u0+h_N, as in0151. For the
reference u0 their exact Euler histories w0 have bounded spatial C1 norms
on[0,T], with constants depending on the fixed K,T,lambda and data, but
not on the larger periodic quadrature cell. The explicit one-coordinate
operator of0151 supplies those bounds on its fixed reference period.

The exact difference e=w_N-w0 solves

    e_t=L_(u_N,K)e
        -P_K[(h_N.grad_K)w0+(w0.grad)h_N].                   (5)

P_K has L2 norm at most one at every frequency, including the defined
harmonic/mean sector. Integration by parts bounds the exact homogeneous
Euler propagator by exp(T||grad u_N||_infinity), uniformly in N by(2).
The forcing norm in normalized cell L2 is at most
C_TK ||h_N||_(H1,average). The Kelvin initial difference has the same
bound; the common-V initial difference is zero. Hence

    sup_[0,T] ||w_N-w0||_(L2,average)
                           <= C_TK sqrt(sum_j |a_j|^2).       (6)

The actual Lin histories satisfy the corresponding estimate after
differentiating their transport equation and retaining the full
material-position term. This is not an isolated-cell pressure estimate.
Whole-field phase/rotation averaging does not increase(6); the reference
constants are uniform on the compact rotation group for the fixed
nonzero wave-number magnitudes.

For the Jacobi action, use the exact material expression
rho<|D eta|^2-eta* Hess(p)eta>/2. Its mixed bilinear coefficients between
the actual prepared columns also converge: D eta=w+(grad u)eta,
global coefficient derivatives are uniformly bounded, and differences
of Hess(p_N) from Hess(p0) have normalized L2 norm bounded by
C_M ||h_N||_(H2,average). Pair the latter with the bounded reference
columns, and use Cauchy--Schwarz for the difference columns. Thus the
action comparison follows from the same actual fields, not an appended
wave Lagrangian. The canonical full-phase mass remains exactly rho.

At each fixed K the right side of(6) can be made smaller than any
prescribed fraction of |K|^2 times the reference response scale.
Together with0151 this transfers its finite-window positive mean-wave
content while preserving the actual inserted tube regions. This statement
does NOT establish a uniform C2 Bloch derivative as P grows, nor a
long-wave limit with |K|P small, nor an acoustic-time uniform error.
Those are distinct homogenization/action licenses. In particular small
bulk L2 energy alone cannot be differentiated twice through an inverse
whose lowest nonzero periodic wave number tends to zero.

## 4. Local optical action is controlled on its own scale

If0147 supplies a whole-space localized actual optical history family
in H^(s+1), with two carrier derivatives and finite positive KKS/action,
the large comparison ball above is chosen using that family's uniform
tail tightness. The field difference is small on that ball and globally
bounded by(2). The full-pressure Duhamel comparison of0145 then controls
each optical history and physical material observation relative to its
own fixed packet normalization. It is not inferred from bulk estimate(6).

A finite periodic realization can use a still larger common supercell.
Periodize the smooth solenoidal Kelvin initial generator, retaining its
actual carrier derivatives. For initially localized smooth data, periodic
and whole-space Leray operators agree locally as the cell grows: their
Green functions differ by a smooth regular term whose second derivatives
are O(P^-3) in three dimensions. On a fixed ball this follows directly
from G_P(x)=P^-1 G_1(x/P) after the common Newton singularity is removed.
Truncate a fixed H^s tail first; the global L2 projector bound controls
the discarded part. The regular-kernel estimate and the same uniform
Euler energy bound then give finite-time local and per-cell action
convergence. No artificial pressure wall is imposed on the packet.

The packet action per unit volume is its actual finite cell action
divided by the supercell volume. Its relative error is unchanged by
that division. Its resulting j is strictly positive for a selected
finite construction, but tends to zero if the cell grows with a fixed
packet. This theorem neither fixes j in that limit nor imports a finite
Cosserat density from a single localized packet. An actual dense
homogenization or scaled-family closure remains a separate achievement.

## 5. Exact scope and active continuation

Established construction: a bounded exact stationary constant-curl field
with prescribed finite local tube regions, an exact periodic approximation
with arbitrarily small normalized bulk defect, and actual fixed-finite-K
Euler mean/action comparison. The local optical transfer consumes0147's
explicit finite-action license; its pending completion is not assumed here.

This is a useful same-field join, not yet the full autonomous Cosserat
continuum. Remaining requirements are the physically normalized joint
optical/acoustic phase action, its genuine long-wave/isotropic closure,
and the requested relation between the optical structure and EPS topology.
The full objective remains active. Next candidates are0151's small
fixed-cell elliptic field (which avoids the diverging-period inverse),
the actual phase-action construction0154, and a spatial-density-scaled
insertion family with its own long-wave bounds.
