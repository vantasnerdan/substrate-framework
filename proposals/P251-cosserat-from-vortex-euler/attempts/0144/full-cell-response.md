# One-array generic-K Euler material action and physical mean response

This is an active construction attachment. Its exact operator and
observation identities do not yet assert an autonomous isotropic
Cosserat equation. Accepted C-CST-008/009/010 are unchanged.

## 1. The same stationary array and the actual three-dimensional operator

Use the exact smooth triangular periodic array of 0139/0141. Its total
planar vorticity is ζ=q_core−Γ/Acell, its planar velocity is v=J∇ψ with
zero harmonic mean, and its reduced pressure p means physical pressure
divided by rho. The exact Bernoulli lift is (v,W), W=sqrt(2(C−B)),
B=p+|v|²/2. Work in the actual axial Galilean frame removing the mean W:
u=(v,h), h=W−〈W〉. Thus 〈u〉=0 and

    div u=0,       (u·∇)u=−∇p.

All spatial averages below are full-fluid primitive-cell averages,
including compensation and ambient fluid. Set K=εκ, with κ fixed and
arbitrary, and Dε=∂t+u·∇+iε(u·κ). A material displacement Bloch amplitude
η obeys the exact Euler/Lin equations

    w=Dεη−(η·∇)u,
    Dε²η+Hess(p)η+∇εp1=0,       divεη=0.                 (1)

Indeed the Eulerian linear equation for w equals
Dε²η−[(u·∇)∇u+(∇u)²]η; stationary Euler turns the bracket into
−Hess(p). This is the full pressure-constrained material Jacobi equation,
with action

    Lε=rho/2 〈|Dεη|²−η*·Hess(p)η〉.                    (2)

It applies to the full Euler phase initial data. A fixed Kelvin leaf is
an additional initial-data condition, not an identity imposed on every
macro velocity. The underlying Euler trajectories conserve their own
circulations. No constant-curl formula is used for the nonconstant-factor
Bernoulli lift.

## 2. Generic three-dimensional cell corrector, including pressure border

Take transverse mean material displacement U, κ·U=0, and write

    η=U+iεχ+ε²η2+O(ε³),   〈χ〉=〈η2〉=0,
    div χ=0,                    div η2=κ·χ.

The last equation is essential. A periodic zero-mean inverse Laplacian
constructs its gradient part; its remaining solenoidal part remains an
actual higher cell state. Put A0=u·∇, D0=∂t+A0 and a=u·κ. The first cell
response, with its actual initial data, is

    P0[D0²χ+Hess(p)χ] = P0 Fκ[U],
    Fκ[U]=−2a Udot+(κ·∇p)U+κ(U·∇p).                    (3)

P0 projects onto periodic solenoidal zero-mean cell fields. Mean
translation and the harmonic pressure constraint are kept separately.
At this order D0 on the columnar cell amplitudes is ∂t+v·∇xy; h enters
the source a, preparation and physical observation, not an unexamined
axial derivative of the cell function.

For clarity, integration by parts gives the pressure border

    〈U·Hess(p)η2〉=−〈(U·∇p)(κ·χ)〉.

Expanding the *full* action (2), not a guessed oscillator, gives after
the exact time boundary rho ε² d〈a U·χ〉/dt,

    Lε = rho |Udot|²/2
       + rho ε²/2 {〈|D0χ|²−χ·Hess(p)χ〉
                       +〈a²〉|U|²+2〈χ·Fκ[U]〉}
       + O(ε³).                                      (4)

The bare covariance term has the sign of negative elastic stiffness.
The actual cell response can change that sign; deleting it or deleting
the last pressure term changes the action. The exact causal solution
of (3), plus its propagated initial cell state, is the full retained
branch memory. No inverse of a static operator or spectral gap is
assumed by writing that solution. The material first-order formulation
through Euler/Lin supplies its well-posed finite-time evolution on the
fixed smooth periodic background.

## 3. Physical momentum removes part, but not all, of this memory

The full Euler Fourier mean m=〈w〉 is not 〈η〉dot. Exact integration of
the Lin relation gives

    m=〈η〉dot+iK_j〈u_jη−u η_j〉,
    m=Udot−ε²〈aχ−u(κ·χ)〉+O(ε³).                      (5)

The physical mean displacement X is defined by Xdot=m with its declared
initial value; replacing X by U before retaining (5) loses real current.
The mean Jacobi equation and stationary Euler give

    Uddot=ε² Pκ{〈a²〉U+2〈aχdot〉
                  +〈(κ·∇p)χ+∇p(κ·χ)〉}+O(ε³),

    mdot=ε² Pκ{〈a²〉U+〈a D0χ+u(κ·D0χ)〉}+O(ε³).        (6)

Pκ=I−κκᵀ/|κ|² is the actual slow pressure projection. These are generic
three-dimensional formulas on the one array. They include filament and
compensation deformation, pressure response, and the full-fluid mass
rho. They do not average planar modes over a set of measure zero.

For example the *physical* negative constant vorticity of a pure planar
base changes under a solenoidal material map by

    δ(−γ ez)=−γ ∂zη,       γ=Γ/Acell.

It is stationary under planar maps only. Holding it fixed in a bending
calculation omits an actual Euler state. Equations (1)–(6) transport
the complete background, so that omission is impossible here.

## 4. Exact axial mean bending response in the pure-planar limit

First set h=0 and κ=ez. This is the actual pure-planar Euler array,
not a wall-supported approximation. Let U be horizontal. The axial
first cell corrector decouples exactly:

    D0²χz=U·∇p=−U·A0 v.                               (7)

For a displacement-phase preparation with static leading U=U0, the
fixed-Kelvin first-jet condition is

    D0χz=−U0·v.                                      (8)

Its initial value can be a genuine physical cell tilt. Let c be a smooth
streamline-invariant cutoff supported in a contractible vortex cell,
0≤c≤1, and r the actual central Cartesian coordinate there. Choose
χz(0)=−c U0·r. Since A0r=v and A0c=0, the required initial rate is
χz,t(0)=−(1−c)U0·v. It is the Kelvin return, not an externally supplied
inertia or a frozen compensating background. Omitting this return adds
a secular homogeneous term to (7).

Using the physical current (5), the entire axial cell history cancels
from the displacement-phase mean acceleration:

    mdot=ε²〈v D0χz〉=−ε² C_v U0,    C_v=〈v⊗v〉.         (9)

Thus the actual physical axial bending stiffness is rho C_v at this
prepared second-jet scope. It is positive on the horizontal plane for
the nontrivial triangular array. Sixfold symmetry makes
C_v=c_b² I2, where c_b²=〈|v|²〉/2. In the dilute core annulus,
v=Γ eθ/(2πr)+regular field, hence

    c_b²=Γ²/[4π Acell] log(b/a)+O(Γ²/Acell)             (10)

for a≪b comparable to a fixed fraction of the cell scale. The coefficient
is the actual full Euler covariance, not a postulated filament tension;
(10) explains its positive logarithmic scale. A cutoff trial action
would give 〈c(2−c)v⊗v〉 instead of 〈cv⊗v〉. That finite collar discrepancy
is precisely part of the retained response, not an exact identification.

This statement is an exact spatial second-jet identity on each fixed
finite time interval. It alone gives no uniform remainder at times
growing like 1/|K|.

## 5. Full-phase common velocity: rho is restored, but memory is explicit

Candidate C uses actual initial Eulerian velocity w(0)=V exp(iK·x),
κ·V=0, together with η(0)=0. Then ηdot(0)=V and the microscopic kinetic
quadratic term is exactly rho |V|²/2. No projection onto an SH fraction
changes that mass. These initial velocities need not all belong to one
common Kelvin leaf; each resulting Euler solution conserves its own
circulation. This is the full phase boundary authorized by the parent
and already present in C-CST-009's common macro V construction.

For axial K, h=0 and horizontal V, the leading U(t)=tV and the exact
first cell initial data are χz(0)=χz,t(0)=0. Let

    R_v(t)=〈v(x)⊗v(X_(−t)(x))〉,

where X_t is the actual planar material flow; equivalently the second
factor is exp(−t A0)v. Direct integration of (7), not a frequency average,
gives

    D0χz=−t V·v+∫_0^t exp(−sA0)(V·v) ds,

    m(t)=V−ε² [C_v t²/2−∫_0^t(t−s)R_v(s)ds]V+O(ε³),

    X(t)=tV−ε² [C_v t³/6
                    −∫_0^t (t−s)² R_v(s)/2 ds]V+O(ε³).       (11)

The same full-fluid density rho multiplies both terms. At t=0 the
correlation R_v(0)=C_v exposes a real initial departure from an immediate
autonomous wave equation. A field-dependent preparation can move that
initial response into a microscopic slip; ignoring it cannot.

On a closed-streamline region, v=A0r is an exact bounded cohomological
identity. It suggests integration-by-parts control of the correlation
integrals in (11), rather than assuming a spectral gap at a separatrix.
The already registered 0146 continuation owns that actual acoustic-time
normal-form/control problem. It must control the true Euler remainder,
not identify a finite-time K² Taylor jet with a result uniform at 1/K.

## 6. Generic-K in-plane transverse polarization

For h=0 split κ=κh+κz ez and choose U horizontal with κh·U=0, i.e.
U parallel to ez×κ. The horizontal part of (3) is exactly the planar
cell equation with κh: it contains neither κz nor χz. Its initial data
still need matching: 0139's compact core-translation generators are not
silently identified with a uniform material-U chart. In particular a
circulation impulse removed by the planar harmonic pressure projection
can have a nonzero projection when κz is nonzero.

There is an explicit full-phase repair directly in physical Euler
initial data. Keep the actual planar SH horizontal initial velocity
w_h,planar(Kh). For the displacement phase add

    w_z(0)=−i Kz U·v,
    δw_h(0)=−Kz² D_h,Kh (Delta_h,Kh)^−1(U·v).          (12a)

The displayed inverse is on the nonzero microscopic Fourier modes;
U·v has zero mean. This is an O(Kz²) bounded cell correction and makes
the complete initial velocity exactly divergence free. Its gradient
energy and changed circulation are retained in the same full phase
action. It is not relabeled as an unchanged fixed-Kelvin preparation.

For a pure-planar base, the full linear pressure source uses only w_h
and horizontal derivatives. Its nonzero-cell Fourier denominator is
|g+Kh|²+Kz². Thus the horizontal evolution has no new first Kz derivative;
(12a) changes it only at second order. The SH projection of its mean
stress needs the horizontal first jet only. The mean harmonic pressure
projection acts identically on SH in the planar and three-dimensional
problems. Consequently the actual planar contribution is retained.

The added axial first row solves D0 Z_z=U·∇p and has Z_z(0)=−U·v;
for the static leading displacement it is exactly Z_z=−U·v. Its physical
mean-stress contribution is therefore −Kz² C_v U. This proves the
matched-state displacement-phase addition without borrowing a planar
mode on a set of measure zero. The restoring second jet is the sum of
the actual planar response and

    −κz² C_v U.

At 0139's point-action/controlled smooth-history boundary this gives
the positive SH benchmark

    c_SH(n,Khat)²=c_p²[1−(n·Khat)²]+c_b²(n·Khat)²,
    c_p²=Γ²/(16π Acell),      c_b²=〈|v|²〉/2.           (12)

The planar coefficient has exactly the finite-core/ordered-limit
license proved in 0139, not a new fixed-core eigenvalue theorem.
Both full initial phases carry their actual circulation/current rows;
the common-V velocity phase in (11) displays the extra correlation
term explicitly. Equation (12) is not permission to discard that term
or to assert uniform acoustic-time propagation.

The lifted h is small in every fixed derivative norm when the Bernoulli
constant is large. Equations (3)–(6) retain its contributions directly.
The finite-time cell coefficients and initial Kelvin/Lin rows depend
continuously on h with controlled O(C_T||h||) errors at fixed microgeometry.
This is a coefficient-level perturbation: zero-mean cell projection
is fixed before differentiating K. It is not an assertion that an
absolute velocity error is o(K²). The actual h-dependent physical
coefficients accompany the comparison.

## 7. Whole-field SO(3)/TR law and the optical interface

Rotate the complete array, ambient compensation, pressure, material tags
and initial data by one R∈SO(3), then average R with Haar law. Time
reversal likewise acts on complete initial Euler states. No interaction
between independently rotated columns is invented. This operation
makes the exact ensemble response rotationally covariant; it does not
automatically remove its temporal cell memory.

If only the SH history is populated, its actual projection is
P_SH(n,Khat) onto n×Khat. The geometric identities are

    E P_SH=P_K/2,
    E[(n·Khat)² P_SH]=P_K/6.

Hence the reconstructed SH scalar speed is
(2c_p²+c_b²)/3. But assigning conditional velocity 2P_SH V to make the
*raw physical mean* equal V gives kinetic energy rho |V|², i.e.
inertia 2rho. Using P_SH V instead gives raw mean V/2. Neither may be
silently called common mass rho. Candidate C restores common full
solenoidal V and mass rho while retaining the complementary
polarization and its response.

One exposing complementary case is h=0, κz=0, V=Vz ez. It is the exact
passive axial perturbation

    w_z(t,x)=Vz exp[−iε κh·(x−X_(−t)(x))].             (13)

Its mean second jet is the actual material-displacement covariance,
not (12). The corresponding displacement-phase fixed-Kelvin generator
parallel to ez has ξ×ζez=0, so it supplies no axial restoring force.
These are rank/initial-data facts about this array, not a parent no-go.

0141 supplies actual optical carrier histories on this *same* base.
Because the base is independent of z, distinct actual axial Fourier
sectors have zero full-fluid quadratic cross-action. A low-K acoustic
sector and optical carriers ±k* do not acquire an energetic cross-term
just by averaging the whole field's orientation. A material tag/hybrid
observation may nevertheless have an actual derivative coupling. It
must retain

    J_E−J_H=−iK_j C_ij−K_jK_l T_ijl/2+O(K³),
    C_ij=I_dot,ij/2−ε_ijm S_m/2.

The optical spin of 0138/0141 enters this identity, with its shape and
quadrupole terms. Changing from the physical Euler mean to a hybrid
centroid coordinate changes the observed transfer and the kinetic
cross-current together. It does not create a missing full-fluid
cross-action. The exact constitutive comparison must keep that map.

## Active continuation and oracle

The full generic-K operator, action, physical mean identity, exact axial
displacement-phase bending stiffness and common-V correlation formula
are the positive results above. The complete two-polarization
whole-field acoustic-time closure, optical/current joining and global
stationary geometry remain distinct active parent obligations.

`jacobi_join_verify.py` evaluates the full action on a genuine
nonconstant periodic Euler field using canonical Fourier algebra. It
exposes the nonzero divergence-pressure term, wrong omission, actual
three-dimensional compensation deformation and axial-sector
orthogonality. Its initial aborted equality was solely SymPy's
unpropagated conjugation of derivatives of declared real functions;
the diagnostic and repaired output are both preserved. Further exact
mean/correlation and orientation checks accompany this attachment.
