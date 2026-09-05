# Same-profile global ring, actual optical pole, and prepared physical current

This continues the two preceding operator files. The global ordinary
Euler ring is constructed here by a bordered integral equation, rather
than attributed to a source theorem with a different core profile.
It is not claimed to be an EPS constant-curl field.
Lengths in logarithms use the fixed seed core length as their unit;
restoring that unit replaces log R by log(R/a_ref) and shifts A only.

## 1. Fixed smooth radial seed and its full-plane border

Use0181's smooth near-Rankine vorticity Z(s): exactly2Omega for s<=a,
smooth monotone taper on [a,c], c=a+epsilon, and zero for s>=c.
Set V(s)=s^(-1) integral_0^s t Z(t)dt and

    phi0(s)=integral_s^c V(t)dt,
    f(phi0(s))=Z(s), f(t)=0 for t<=0.

The taper is flat at its endpoints. Consequently f is smooth,
nondecreasing and bounded, with a constant plateau at large positive
arguments. Its exact full-plane equation is -Delta phi0=f(phi0).
Write Gamma=integral_R2 f(phi0)>0 and
K0 g=-(2pi)^(-1)integral log|x-y| g(y)dy. Then
phi0=K0 f(phi0)+A0 for a uniquely chosen constant A0.

All nonlinear unknowns below are restricted to one disk D_b with b>c,
where phi0(b)<0. The global extension is constructed from the integral,
not an imposed wall at b. Work in C^(2,alpha)(D_b), even in z.

The radial linearized equation is Delta h+Q(s)h=0,
Q=f'(phi0)=-Z'/V>=0. For angular m>=1, V is the positive m=1 solution.
The ground-state identity gives only the two translation kernels at
m=1 and no decaying harmonic kernel at m>=2. The z-even restriction
retains only partial_x phi0 from those two kernels.

The radial regular solution h0(0)=1,h0'(0)=0 has a nonzero logarithmic
coefficient outside c. This is proved directly in the thin taper:

    (s h0')'=-s Q h0,
    integral_a^c s Q(s)ds ->2,
    sup_[a,c]|h0-1|=O(epsilon).

Thus its exterior coefficient s h0' tends to -2, and is nonzero for a
fixed sufficiently thin smooth taper. This is an analytic mass-border
margin, not a sampled soft eigenvalue. Imposing integral Qh=0 removes
this radial homogeneous solution; the remaining constant row fixes A.

The translation cokernel of I-K0 Q is Q partial_x phi0. Its pairing
with the drift source x is

    integral x Q partial_x phi0
                   =integral x partial_x f(phi0)=-Gamma !=0. (1)

The same identity shows that the center row integral x Qh detects the
translation kernel. Therefore the complete bordered linear map

    (h,delta A,delta U) ->
    (h-K0 Qh-delta A+delta U x,
                         integral Qh, integral x Qh)        (2)

is invertible. To see surjectivity, I-K0Q is Fredholm of index zero on
the stated bounded disk; K0Q gains two derivatives. The angular kernel
classification and the nonzero radial/translation borders just proved
give a zero kernel for the index-zero bordered map. This is a global
logarithmic matching argument, not the unrelated Dirichlet invertibility
of0136's local torus.

## 2. Actual whole-space stationary ring equation

Let r=R+x and r'=R+y_x. For a toroidal vorticity
omega_phi=(r/R) f(phi), the exact axisymmetric vector-potential formula
gives the streamfunction psi=R phi plus its uniform-flow/background terms.
Its integral kernel on D_b is

    K_R(x,y)=r r'^2/(4pi R²) integral_-pi^pi
                   cos t / sqrt(|x-y|²+4rr' sin²(t/2)) dt.   (3)

The meridional velocity of psi is u_r=-psi_z/r,u_z=psi_r/r.
This normalization gives omega_phi/r=f(phi)/R, a genuine first integral.

Put L_R=log(8R)-2 and Ktilde_R=K_R-L_R/(2pi).
The complete elliptic-integral expression for the angular integral is

    4/D [(2/k²-1) K(k²)-2 E(k²)/k²],
    D=sqrt(|x-y|²+4rr'), k²=4rr'/D².

Its near-diagonal expansion, or the same near/far integral split, gives

    K_R=(2pi)^(-1) (r/R)^(1/2)(r'/R)^(3/2)
                        [log(8sqrt(rr')/|x-y|)-2]
          +O(|x-y|² R^(-2) log(R/|x-y|)).                    (4)

In particular Ktilde_R->K0 with O(log R/R) local elliptic-operator norm.
This means the maps from a bounded C^(alpha) source to C^(2,alpha)
potential on D_b, not a uniform pointwise assertion about a differentiated
logarithmic kernel. The leading logarithmic singular coefficient is
smooth and converges; the remainder has the usual local elliptic bounds.

Solve the following actual core integral equation and two finite borders:

    phi=Ktilde_R f(phi)+A-U x-U x²/(2R),
    integral_D_b f(phi)=Gamma,
    integral_D_b x f(phi)=0.                               (5)

The small-parameter IFT with the invertible map(2) supplies a solution
(phi_R,A_R,U_R) converging to(phi0,A0,0) for every sufficiently large
finite R. Smooth dependence is in the unknowns; continuous small
parameters log R/R suffice, and finite-R differentiation is ordinary.
The leading x coefficient of(4) and (1) give

    U_R=Gamma log R/(4pi R)+O(1/R)>0.                        (6)

Indeed the log R/R part of K_R is
(4pi R)^(-1)log R [x integral f+3 integral y_x f].
The borders in(5) leave precisely Gamma x log R/(4pi R).
Its projection onto the translation cokernel fixes(6). Nonlinear
quadratic errors are O(log² R/R²)=o(1/R), so they do not change this
sign. The construction determines the actual translation speed, not
a frequency fit.

Define the GLOBAL function by the same complete kernel and

    phi_global=K_R f(phi_R)-mu_R-U_R r²/(2R),
    mu_R=L_R Gamma/(2pi)-A_R-U_R R/2
                      =3 Gamma log R/(8pi)+O(1)>0.          (7)

It agrees with phi_R on D_b. It is negative on the meridional boundary
of that disk, on the symmetry axis and at infinity. Outside its source
it solves Delta*phi_global=0, where
Delta*=partial_r²-r^(-1)partial_r+partial_z² and Delta* r²=0.
The exterior maximum principle therefore makes it nonpositive there.
There is no extra positive island or unregistered source. Hence f of
the global function is exactly the source in(5).

The physical streamfunction psi=R phi_global yields

    -Delta*psi=(r²/R) f(psi/R),
    omega_phi=(r/R)f(psi/R), u_phi=0,
    u -> -U_R e_z at infinity.                              (8)

The exact identity omega_phi/r=f(psi/R)/R proves steady Euler. Its
pressure is the actual Euler pressure. The smooth flat f makes the
velocity and vorticity smooth across the core boundary, and the
vorticity has compact toroidal support. Elliptic regularity upgrades
the fixed-taper convergence to every required finite derivative order.
The nested streamlines and nondegenerate core persist by C2 convergence;
their action-angle map and positive period bound converge as well.
This constructs the same-profile family needed in the operator proof.

The constant value psi=-R mu_R on the axis is a streamfunction gauge,
not a singular physical vector potential. For a regular global vector
potential use psi_phys=psi+R mu_R (zero on the axis), and write the
vorticity law as f(psi_phys/R-mu_R). The velocity is defined by the
displayed derivatives and the full-space source potential, so there
is no added vortex filament on the symmetry axis.

The background has a uniform far-field velocity in its steady frame.
In a frame with fluid at rest at infinity it is a traveling vortex ring.
No global finite kinetic energy for the added uniform velocity is
claimed; perturbation velocity and the phase action are finite.
Use the rest-at-infinity velocity for its finite energy E and the
compact-vorticity translation moment I_z=(rho/2)integral(x cross omega)_z.
The actual steady-frame Hamiltonian is E-U_R I_z. Indeed, for a compact
generator, delta I_z=rho integral(xi cross omega)_z, which is precisely
the translation KKS row. Thus the frame subtraction is the Euler
translation moment, not an imported internal inertia. The centered
tag spin and covariance angle are unchanged by the common Galilean
translation.

## 3. The actual closed-ring optical mode

Apply Sections2--5 of full-poloidal-operator.md and the full kernel
estimate in toroidal-kernel-transfer.md to this fixed smooth family.
The positive simple column eigenvalue continues for sufficiently large
finite R and integers n with n/R near the selected carrier. Its actual
Euler velocity and pressure are global normal modes, not packets whose
residual has been called a pole. The compact vorticity source supports
the full KKS, with the translating-frame momentum term retained.

Mode comparison uses one fixed local amplitude, or equivalently the
cross-sectional norm per toroidal length. The full finite ring action
then has its actual2pi R factor. A finite-arc tag need not have a spin
overlap bounded away from zero as R tends to infinity; positivity is
claimed at the selected finite R, with no length factor discarded.

The full positive-mode contour is fixed before R. Thus eigenvectors,
KKS and phase energy are continuous. The three actual adjacent harmonic
operators have the controlled discrete carrier jets of the kernel proof.
Their positive action and finite gap persist. The squared-frequency
curvature inherits the column sign in this discrete jet sense. This
does not identify those three ring harmonics with a generic continuous
laboratory Bloch wavevector in a many-ring medium.

## 4. Literal stationary current from an actual prepared TR pair

Choose a finite stationary positive material fraction on the ring,
chi(I,phi), supported in the rigid-seed tag's continued annular regions
and a finite arc. Its reference centroid and covariance are the actual
Euclidean ones. A simple distinguished covariance axis defines one
physical tilt component theta with exactly unit rigid-rotation response.
The chart is nondegenerate by convergence to0181's elongated tag.

For one actual ring eigenmode let nu>0 be the magnitude of its frequency.
After choosing the phase by the measured angle, its material displacement
on the active core has the two real columns

    xi_+(t)=xi_q theta(t)+xi_p theta_dot(t),
    theta_tt+nu² theta=0.                                  (9)

Take the actual time-reversed stationary field -u_R, with the same tag
and the reflected phase preparation

    xi_-(t)=xi_q theta(t)-xi_p theta_dot(t).                 (10)

These are registered common INITIAL angle/rate data, not an arbitrary
instantaneous identification of two off-shell fluid paths. For each
initial pair(theta0,v0), both actual histories have exactly the same
theta(t)=cos(nu t)theta0+sin(nu t)v0/nu. Their conjugate fluid phase
coordinates are related by time reversal, not all tied identically.

Let rho be the actual density and project cross products on the chosen
physical tilt axis. Set

    G_q=rho integral chi r_centroid cross xi_q,
    G_p=rho integral chi r_centroid cross xi_p,
    C_q=rho integral chi xi_q cross u_R,
    C_p=rho integral chi xi_p cross u_R.

The exact invariant-tag material identity gives

    G_+=G_q theta+G_p theta_dot,
    S_+=G_q theta_dot+G_p theta_tt+2C_q theta+2C_p theta_dot,
    G_-=G_q theta-G_p theta_dot,
    S_-=G_q theta_dot-G_p theta_tt-2C_q theta+2C_p theta_dot.

Consequently the equal-weight physical law has

    G_average=G_q theta,
    S_average=(G_q+2C_p)theta_dot.                          (11)

At the0181 column reference, C_p is a nonzero smooth normalization
times its radial axial-stretch moment R2. Its positive three-control
annular Jacobian therefore controls C_p and both carrier derivatives.
The carrier-to-radial-parameter derivative is nonzero here, not assumed:
for M=1 put k=-1-x K0(x)/K1(x). The modified-Bessel Riccati equation
gives x k'-2k=2+x²[1-(K0/K1)²]>2, since the positive integral
representations give0<K0<K1. At fixed l,
D_x=(l²/x³)[x k'-2k-1/sqrt(1+l²/x²)]>0. With D_l<0 this implies
l'(x)>0. Thus the three carrier rows are an invertible triangular
change of the0181(F,F_l,F_ll) rows, including their tag variations.
For the actual ring, use the THREE values at n-1,n,n+1, or equivalently
their scaled divided differences. The Jacobian converges to that same
invertible three-row matrix. The pressure and eigenvector convergence
through discrete jets, together with smooth physical moment maps, makes
the finite-dimensional positive-tag IFT applicable. It adjusts the two
annulus centers and one positive ratio and sets C_p=0 at all three
actual harmonics with ONE stationary tag. All fractions remain positive;
the centroid/covariance chart and tilt response remain nonzero.

Since G_q tends to Delta>0, (11) yields at each of the three modes

    G_average=J_tag theta,
    S_average=J_tag theta_dot,
    G_average(0)+integral_0^t S_average=J_tag theta(t),
    J_tag=G_q>0.                                          (12)

This is exact for every time of each linear normal mode on its actual
prepared TR law. Each individual realization retains its G_p and C_q
rows; they have been averaged explicitly, not declared zero separately.

## 5. Same initial phase form and energy, before elimination

Let the actual positive ring-mode two-form on the initial physical data
be M dtheta0 wedge dv0. The negative-background KKS reverses sign,
while its conjugate phase column in(10) also reverses sign. Thus its
pullback to the SAME(theta0,v0) is again M dtheta0 wedge dv0.
The phase energy is M(v0²+nu² theta0²)/2 in both preparations.

Average these actual full initial forms and energies with weights1/2
before eliminating a phase coordinate. They remain exactly those same
forms, giving

    L=M(theta_dot²-nu² theta²)/2, M>0.                      (13)

No first-order symplectic form was averaged to zero by tying both
conjugate fluid coordinates. No geometric rotor inertia was added.
The measured overlap J_tag/M is positive and is not required to equal1.
Equations(12),(13) supply a genuine closed-ring physical current/action
input to the parent's coherent preparation route. A common continuous
laboratory K, full many-body mean reaction, and EPS strong-Beltrami
geometry remain separate joins. This ordinary no-swirl ring is neither
a constant-factor nor a generalized-force-free field. The two-field
BACKGROUND time-reversal law is stationary in time; its prepared
oscillatory histories are not asserted to be a stationary process.
One selected finite ring has not thereby become a spatially homogeneous
positive-density ensemble.
