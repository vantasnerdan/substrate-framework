# An exact stationary elliptic tube near the one-wave acoustic field

This executes candidate B. It is a closed tube on a flat periodic
three-torus, not a claim of a compact Euclidean EPS knot or a completed
intrinsic angle/action normalization.

## 1. Exact constant-curl insertion with a small global norm

In dimensionless coordinates of period 2pi, choose 0<d<1/10 and

    u_d=(cos z+d sin y, sin z, d cos y).

Direct differentiation gives div u_d=0 and curl u_d=−u_d. Consequently
it is exact stationary Euler with p_d=−|u_d|²/2; no superposition of
nonlinear Euler solutions is presumed. Relative to the circular
one-wave u_0=(cos z,sin z,0), its difference is O(d) in every fixed
C^r norm. Restoring dimensions uses u=v u_d(qx), curl u=−q u.

The transverse motion has the exact first integral

    H(y,z)=−cos z−d sin y,
    ydot=H_z=sin z,   zdot=−H_y=d cos y,
    xdot=−H.

Thus y=pi/2,z=0 is an actual periodic core streamline, traversing the
x circle with speed 1+d. Its transverse Hessian is diag(d,1), and
the linearized transverse matrix is [[0,1],[-d,0]]. The physical
Floquet multipliers on one x return are

    exp[±i 2pi sqrt(d)/(1+d)],

strictly elliptic for the declared interval. Small positive levels of
E=H+1+d below the first transverse saddle form nested invariant
circles; crossing them with the x circle gives a finite solid vortex
tube. The boundary is a material vorticity surface because curl u is
parallel to u. This is actual finite-core geometry, not a scalar phase
declared to be a material angle.

There is also a nonzero small-amplitude twist. Put Y=d^(1/4)(y−pi/2),
Z=d^(−1/4)z. This is a symplectic transverse change of coordinates.
With I=(Y²+Z²)/2, the quartic Birkhoff average is

    E=sqrt(d) I−(1+d) I²/16+O(I³).

The transverse frequency divided by the actual x speed has derivative

    d/dE [nu(E)/(1+d−E)] at E=0
       =[8d−(1+d)²]/[8sqrt(d)(1+d)²],

which is nonzero on 0<d<1/10. The exact invariant level tubes already
exist without invoking an external KAM theorem. The twist is recorded
as useful data for subsequent local persistence/embedding work.

## 2. The actual fixed-time acoustic response survives the insertion

Keep d, the periodic cell and a finite time T fixed before macro k.
The same two physical initial columns used in `one-wave-response.md`
are defined from the actual u_d and its curl, together with the common
initial velocity V. The complete Euler operator is

    L_(d,k)w=−P_k[(u_d·∇+ik κ·u_d)w+(w·∇)u_d].

It is not replaced by the one-wave operator. On periodic H^r spaces,
the usual differentiated transport energy estimate follows directly
by integration by parts and the Leray bound. The first two k
derivatives are obtained by differentiating this equation: derivatives
of each nonzero Fourier projector are bounded by constants times
|g|^(−j), j≤2, while the zero Fourier projector is the fixed Pκ.
The minimum nonzero microscopic wave is fixed at q, uniformly under
whole-field rotation. There is no large-cell inverse-limit shortcut.

The operator difference is O(d) from H^(r+1) to H^r. Duhamel's formula,
the transport estimates, and the same estimates for two k derivatives
give a C²-in-k finite-time difference bounded by C_(r,T)d on smooth
prepared data. One extra derivative closes each differentiated source;
all fields here are finite Fourier and hence smooth. Constants are
uniform for 0≤d≤1/10 and for the whole SO(3)/phase law.

For EVERY stationary field in this comparison, the zero spatial jet is
the exact translated/Galilean solution

    w_0(t)=V−[(D+tV)·∇]u_d,
    m_0=V.

The first physical mean spatial jet vanishes for all t: its mean stress
is the sum of a constant times 〈u_d〉=0 and a periodic total derivative
of u_d⊗u_d. Hence the response difference begins at k², rather than
being only a small absolute O(d) error. Together with the uniform
third-spatial-derivative finite-time estimate on these smooth fixed
cells, this proves, in the scales of the previous attachment,

    sup_[0,T]|Xbar_(d),tt+k²(2v²/15)Xbar_d|
      ≤ k²v² C_T [d+|k|/q]
          (|D|+|V|/(qv)).                              (1)

The exact initial isotropic coefficient is instead
(2/15)〈|u_d|²〉=(2v²/15)(1+d²); either positive reference differs from
(1) by a smaller controlled term. Equation (1) licenses the ordered
choice: prescribe a fixed-time relative response accuracy, choose a
strictly positive d small enough, then choose k/q small enough.
Each selected finite d has a nondegenerate elliptic vortex tube; it is
not necessary to take a degenerate tube limit as the constructed object.

## 3. Actual action and remaining same-field joins

The same-field material Jacobi action is

    L_d=rho/2〈|D_(d,k)η|²−η* Hess(p_d)η〉.

Its coefficients, actual physical-current rows, and finite-time
prepared response depend continuously on d in the same fixed-cell
norms. This supplies an action-level comparison with retained corrector
and endpoint data, not a license to drop them. Common physical V still
has mass rho. None of these estimates supplies a nonzero optical tag
inertia, an absolute material angle, or its coupling to X.

The small-d core is strongly strained: its transverse matrix has
rotation magnitude (1+d)/2 and strain magnitude (1−d)/2. Therefore a
near-axisymmetric optical-core theorem cannot be inserted without a
new actual mode calculation or a different localized construction.
Likewise these winding torus cores are not closed curves in Euclidean
R³ before periodic identification. The parent's bounded Herglotz/CK
insertion and global-torus constructions address that distinct spatial
embedding/optical compatibility achievement. This attempt establishes
one concrete nonuniform stationary geometry compatible with the
positive controlled acoustic window and leaves those joins explicit.
