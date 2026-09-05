# Full-poloidal spectral isolation and the closed-ring continuation operator

## 1. Result and exact scope

The next positive step is not merely axial quantization. The0181 nodal
column carrier can be chosen so its positive laboratory-action eigenvalue
is isolated in the FULL poloidal-channel active-vorticity operator. This
property survives the actual smooth compact-core taper. A closed ordinary
vortex ring has an exact transport-plus-compact-pressure decomposition in
each toroidal harmonic, with an explicit transport-band resolvent. This
supplies the correct spectral continuation space and contour test.

No mode-bearing closed ring or EPS pole is inferred from these statements.
Section6 records the exact remaining same-background construction and a
different global generalized-Beltrami source candidate. The transport
space includes every excited poloidal channel and the complete exterior
velocity, not an artificial impermeable wall.

## 2. Full straight-column operator with compact active vorticity

Let u=r O(r)e_theta, Z(r)=2O+rO', with Z smooth, nonnegative,
nonincreasing and supported in r<=c. It equals2Omega in the inner core;
the taper is thin enough that O>=Omega_min>3Omega/4 on the active disk.
Fix axial p>0. The phase is exp(i m theta+i p z-i omega t).

The Hilbert space is

    H_p={w in L2(R2;C3): supp w subset disk_c,
                           div_p w=0 distributionally}.

Velocity is the actual complete-space inverse curl

    v=B_p w=curl_p(-Delta_xy+p²)^(-1)w.             (1)

It has its actual exponentially decaying exterior tail. The bounded
active support does not truncate that velocity. At fixed p, B_p maps
L2 to H1, hence restriction to disk_c is compact. This is a local
Rellich statement on the vorticity source space, not compactness of
the whole-space Leray projector.

Linear vorticity evolution is

    L_p w=T w+K_p w,
    T w=-[u,w], K_p w=[Z e_z,B_p w]
                       =i p Z B_p w-Z'(B_p w)_r e_z. (2)

It preserves H_p: Lie transport preserves divergence and the invariant
active support, while the last bracket is supported where the background
vorticity or its derivative is nonzero. Every exterior pressure effect
is included through(1). The operator K_p is compact for each fixed smooth
taper; no taper-uniform bound on Z' is asserted.

In cylindrical COMPONENTS the rotation of the basis cancels the solid
rotation part of Du, giving

    T=-O(r) partial_theta+N(r),
    Nw=r O'(r) w_r e_theta, N²=0.                  (3)

For a spectral parameter z outside the full transport bands

    S_T=closure union_{m in Z}{-i m O(r):0<=r<=c},

the m-th resolvent is exactly

    (z-T_m)^(-1)=d_m^(-1) I+d_m^(-2)N,
    d_m=z+i m O(r).                                (4)

Its norm is bounded by distance(z,S_T)^(-1) plus
||N|| distance(z,S_T)^(-2), and its angular derivative is bounded:
large |m| gives |d_m|>=|m|Omega_min-|z|. The natural domain is
{w in H_p:partial_theta w in L2}; the resolvent preserves H_p because
it is the resolvent of actual divergence-preserving Lie transport.

Analytic Fredholm applied to I-K_p(z-T)^(-1) now proves discreteness
and finite algebraic multiplicity of all full-channel Euler eigenvalues
off S_T. A right-half-plane resolvent exists from the bounded-gradient
transport group and bounded K_p, so the analytic Fredholm alternative
is the meromorphic one, not the identically singular possibility.

In particular the0181 frequency interval

    -Omega/4 < omega < 0

lies between m=0 and m=-1 transport bands, away from either endpoint.
Its selected fixed negative frequency therefore has a finite-rank Riesz
cluster on the FULL space H_p. This already repairs the missing full-
poloidal isolation license; it is not isolation only in the m=1 ODE.

## 3. Choosing a simple positive cluster, without a hidden Krein collision

For a Rankine reference write j_M(l)=l J_M'(l)/J_M(l),
k_M(x)=x K_M'(x)/K_M(x), x=p a, M>=0. For negative core Doppler
sigma=-tau Omega, tau in(0,2), l=x sqrt(4-tau²)/tau, the determinant is

    D_M=j_M(l)+(l²/x²)k_M(x)+M sqrt(1+l²/x²).       (5)

The radial Sturm identities used in0181 hold for every integer M:
j_M'<0, k_M'<0, k_M<-M (also k_0<0). Thus D_l<0 in every interval
between pressure poles. At fixed tau the derivative with respect to x is

    h j_M'(hx)+h² k_M'(x)<0,
    h=sqrt(4-tau²)/tau.

Since l_tau<0, D_tau>0. Every fixed nodal branch tau_{M,j}(x) is
strictly increasing and analytic. It tends to2 as x tends to infinity:
its l root stays in its finite Bessel interval while x grows. The
surface branch when present has the same large-x limit. For M=0 the
first interval starting at l=0 has no positive root; the nodal intervals
between positive pressure zeros still have exactly one each.

Choose0181's m=1 nodal branch tau_1=1+d with0<d<1/4. The only
possible competing channels at omega=-d Omega are m=0,1,-1,-2.
Here is a useful exclusion that ALSO holds after monotone smoothing.
The actual radial system for any m has the form

    f'=-(1/r+q)f+bP, P'=c0 f+qP,
    q=2m O/(r s), b=(m²/r²+p²)/s²,
    c0=s²-2O Z, s=omega-m O.

Consequently

    [r Re(f conjugate(P))]'=r b |P|²+r c0 |f|².    (6)

For m>=2, s²>4O²>=2O Z on the core, since O'<=0.
For m<=-3 the same holds because |omega|<Omega_min and
|s|>2O there. The exterior potential matching gives at r=c

    P=i s phi, f=i phi'/s,
    Re(f conjugate(P))=Re(phi' conjugate(phi))<0,   (7)

where phi is proportional to the decaying K_|m|(p r). The center term
in(6) vanishes. Equations(6),(7) contradict a nonzero mode. This argument
uses only the exterior boundary value; it does not integrate a singular
material displacement through any exterior particle critical radius.

At the Rankine reference the remaining coincidence equations are

    m=0:  tau_1-tau_{0,j}=1,
    m=-1: tau_1+tau_{1,j}=2,
    m=-2: tau_1+tau_{2,j}=3.                        (8)

None is an analytic identity: their left sides tend respectively to0,4,4
at large x. On a compact carrier interval away from omega=0, only finitely
many radial branches can enter (8), because tau has a positive lower
bound there and l<=2x/tau. Their coincidence sets are discrete. Choose
one carrier in0181's open positive-curvature interval avoiding them.
The m=1 root itself is simple by D_l<0. The full eigenvalue is therefore
simple and has positive laboratory action, not merely one positive
vector inside an indefinite coincident cluster.

For the finite remaining m channels,0181/0137's full first-order radial
transfer has no Z' coefficient and converges through the thin annulus.
The finitely many relevant real-frequency determinants therefore retain
their noncoincidence for sufficiently small smoothing width. The other
m channels cannot coincide at that real frequency by(6). Fix this
smooth taper FIRST. Section2's Fredholm theorem then supplies a genuine
full-operator contour isolating its simple eigenvalue; it need not be a
contour uniform in a vanishing taper. In particular(6) was not used as
a claim excluding every nearby complex high-m frequency. This gives an
actual smooth column with a FULL-poloidal simple positive-Krein
eigenvalue. Its carrier second derivatives and0181's tag are unchanged.

For analytic carrier comparison one may work on the fixed all-vector
L2 source space with the same extended formula(2). Divergence of a
nonzero eigenvector satisfies (lambda+u.grad)div_p w=0, including any
distributional boundary trace. Off the scalar transport bands this
has only zero solution, by Fourier division by the nowhere-vanishing
lambda+i m O(r). Thus the isolated eigenvectors are automatically in
the physical solenoidal space H_p. This avoids silently identifying
different p-dependent divergence constraints when comparing contours.

## 4. The exact closed-ring operator, not a local-wall replacement

Now let u be an actual smooth axisymmetric no-swirl stationary Euler ring
in its fixed Galilean translating frame. Its vorticity is supported on
a bounded solid torus D, away from the symmetry axis; its meridional
streamlines in D are nested regular closed curves about a nondegenerate
elliptic center. Use the actual volume form r dr dz dphi.

The meridional action is the enclosed r dr dz area divided by2pi.
The angle is normalized time along the actual closed orbit. This gives
volume coordinates (I,theta,phi) with

    u=Omega(I) partial_theta,
    omega=F(I) partial_phi.                         (9)

The second identity is the actual steady no-swirl law omega_phi/r=F(psi).
The coordinate map is regular in polar form at the elliptic center;
use a smooth disk chart there, not an unweighted singular dI basis.
Assume Omega has a positive lower bound on the complete active core.
This hypothesis is satisfied by a sufficiently C2-close thin-ring
continuation of the fixed near-Rankine column; mere nested topology
without a nondegenerate center would not by itself supply the bound.

For a fixed nonzero toroidal Fourier harmonic n, let H_n consist of
global divergence-free L2 vorticities supported in D and equivariant
with that harmonic under axial rotations. Its induced velocity is

    v=B_n w=curl(-Delta_R3)^(-1)w,                  (10)

with its complete free-space tail. The local H1 bound follows from the
Fourier multiplier and compact support: the high-frequency gradient is
Calderon–Zygmund bounded, while the low-frequency velocity integral is
finite in3D using ||w||_1<=|D|^(1/2)||w||_2. Local Rellich makes
B_n:H_n->L2(D) compact. Full cylindrical-vector harmonics, including
the n+1 and n-1 Cartesian components, are part of this inverse.

In coordinate-vector form the exact Euler generator is

    L_n w=-Omega partial_theta w+Omega' w^I partial_theta
                       +i n F B_nw-F'(B_nw)^I partial_phi. (11)

This is obtained directly from -[u,w]+[omega,v]. Both brackets preserve
divergence and the active support. It includes poloidal stretching,
every poloidal harmonic and the pressure reaction from all of R3.

The last two terms are compact; the first two have precisely the
resolvent(4), with the harmless norm-equivalence factors of the physical
disk chart. Hence every point outside

    closure union_m {-i m Omega(I)}                 (12)

is Fredholm of index zero, and the spectrum there is discrete unless
absent. This is a theorem for the actual closed-ring operator under
the stated stationary-background hypotheses. Toroidal quantization
n alone did not prove it: bounded active vorticity, the full inverse
curl and the explicit transport resolvent did.

## 5. What actually transfers an eigenvalue and its action

Let a family of these GLOBAL stationary rings have, after physical
pullback to one meridional disk, C2 velocity/profile and action-angle
convergence to the fixed smooth column of Section3. Take integers n_R
with n_R/R->p0. To turn the Fredholm step into a pole continuation one
needs the complete resolvent comparison, for example on the fixed
positive-mode contour C:

    sup_C ||(L_R-L_0)(z-L_0)^(-1)|| < 1.            (13)

Then the ordinary resolvent identity and Riesz integral preserve the
simple projection. A Hamiltonian perturbation with continuous nonzero
KKS on that eigenspace preserves its definite energy and imaginary
generator eigenvalue. This last conclusion uses positivity of the
whole isolated cluster proved in Section3, not a generic claim that
all vortex-ring perturbations are stable.

The actual Kelvin leaf can be reconstructed on the active core rather
than assumed from an arbitrary vorticity eigenvector. Write ad_u xi=[u,xi]
and let lambda be the mode generator eigenvalue outside the transport
bands. Solve (lambda+ad_u)xi=v on D. Commutation [u,omega]=0 gives

    (lambda+ad_u)[omega,xi]=[omega,v].

The eigenmode w satisfies the same equation, so w=[omega,xi]. Its
divergence vanishes by the same transport inverse, hence xi has zero
total boundary flux and admits a smooth divergence-free exterior
extension. Because omega vanishes outside D, any such extension gives
the same orbit tangent and compact KKS. The actual exterior material
displacement may need time-dependent relabeling to satisfy full Lin;
its velocity and pressure are still the exact Euler normal mode. A
monochromatic L2 material label displacement in the uniform exterior
is not part of this conclusion.

The ring's translating frame is a genuine Galilean frame. Its stationary
phase Hamiltonian includes the corresponding Euclidean momentum
subtraction, not an assigned change of internal frequency. The material
angle and spin are centered observables and are invariant under that
common translation. At the straight-core limit the same-frame action
therefore tends to0181's action. An action-continuity use of(13) retains
this frame/momentum term explicitly.

The local scalar Green comparison entering(13) has the actual
meridional operators

    A_{R,n}=-partial_x²-(R+x)^(-1)partial_x
                           -partial_z²+n²/(R+x)²,
    A_infty=-Delta_xz+p0².                          (14)

The vector inverse uses n,n+1,n-1, not just one scalar potential.
On a fixed compact neighborhood their coefficients converge. To earn
operator-norm rather than pointwise convergence, a growing buffer and
the coercive n²/r² barrier between that buffer and r near0 or infinity
must control the complete inverse. The companion
toroidal-kernel-transfer.md executes that full inverse-curl estimate
directly, including two discrete carrier jets. Section4 alone does not
prove(13); the remaining premise is the global stationary same-profile
background and coordinate convergence.

Similarly, discrete values p=n/R do not define a continuous physical
Bloch parameter on one closed ring. A second-jet conclusion needs the
differentiated meridional resolvent family (or controlled discrete
second differences) and the corresponding physical observable map.
No O(1/R) localization error is divided silently by a gradient O(1/R²).

The actual material tag must also be transported on THIS ring. No-swirl
gives stationary positive weights chi(I,phi), including a finite arc.
This supplies an actual finite stationary reference shape. For its mode
xi, the literal rows are

    G=rho integral chi r_centroid cross xi,
    S=G_t+2rho integral chi xi cross u.              (15)

They include the moving-position term. Equation(15) is exact for an
invariant density; it is not automatically the0181 rigid-core relation.
The new mode and the positive tag moment map must continue together.
Its homogeneous covariance-angle normalization and the full KKS remain
measured quantities, not a supplied rotor mass.

## 6. Source licenses and the next actual geometry candidate

[Cao–Lai–Qin–Zhan–Zou, arXiv2206.10165, Proposition4.1](https://arxiv.org/pdf/2206.10165)
provides global no-swirl thin rings for nonnegative nondecreasing bounded
vorticity laws vanishing on the negative half-line. Smooth flat laws
meet those hypotheses and elliptic regularity supplies smooth velocity.
The proposition gives concentration and diameter control; it is not a
norm-resolvent theorem for0181's chosen radial profile. Its nonlinear
stability conclusions for other specified laws are not imported here.

Thus Candidate A has ESTABLISHED the full-channel positive column
isolation and the exact closed-ring Fredholm/domain step. The next
construction at this stage was a same-profile stationary thin-ring
family with(13), then the stationary material-moment continuation(15).
The companion global-ring-and-tag.md executes those steps by a bordered
full-space Green equation and an actual prepared time-reversal law.

Candidate B has a distinct useful source:
[Abe, arXiv2008.09345, Theorems1.1/2.1 and Section2.4](https://arxiv.org/pdf/2008.09345)
constructs global generalized-Beltrami rings in a nonzero uniform
far-field frame, with compact vorticity and nested invariant tori. The
factor is nonconstant; its displayed power law yields finite regularity
depending on the chosen exponent, not an automatic C-infinity limit.
The theorem does not specify a near-Rankine optical spectrum or thin-
core branch. It is nevertheless a genuine stationary closed-tube Euler
ambient, unlike a local toroidal Dirichlet field.

For a smooth integrable ring with swirl, straightening the periodic
toroidal drift gives u=Omega(I)partial_theta+U(I)partial_phi. The
transport bands become m Omega(I)+n U(I). Its vorticity is tangent to
the same surfaces. The term omega_theta partial_theta B_n is order zero
but relatively compact against the transport resolvent: commute
partial_theta through B_n, use compactness of B_n and of its first-order
commutator, and boundedness of partial_theta(z-T)^(-1). The commutator
kernel gains one difference of the smooth coordinate-vector coefficient,
so it retains inverse-curl order minus one. This supplies the same
Fredholm construction if Omega_min>0 and the chosen contour avoids ALL
shifted bands. Neither condition is supplied just by the word Beltrami.

An analytic noncompact constant-lambda EPS field is a DIFFERENT ambient.
There omega is nonzero outside a selected torus. The full tail v then
creates exterior vorticity via [omega,v], so H_n supported in that torus
is not invariant in general. Local approximation cannot import the
compact-support Fredholm theorem. Its exterior spectral response or a
genuinely reducing symmetry space is a remaining construction, not a
refutation of the EPS route or the parent objective.
