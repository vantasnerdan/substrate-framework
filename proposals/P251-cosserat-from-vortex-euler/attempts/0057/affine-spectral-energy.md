# Exact same-ensemble affine shear and its retained reaction correction

Let F be a constant deformation with det F=1. Transport the actual
stationary vorticity field by `omega_F(x)=F omega0(F^-1 x)`. This is an
exact volume-preserving isovortical deformation, not a supplied elastic
energy. It remains a stationary random field (its perturbed state need
not itself be a stationary Euler solution). Its velocity is obtained from
the full-space inverse curl with the ensemble mean fixed separately.

For an original spectral wave q put p=F^-T q. In material coordinates
x=F y its velocity multiplier on an arbitrary divergence-free input is

    M_F(q) u(q)=-[p]_cross F [q]_cross u(q)/|p|²
              =P_p F^-T u(q).

The equality is the cofactor cross-product identity and det F=1. Thus
`e(F)=rho E|M_F u0|²/2` includes the exact projection; no boundary plane
wave or finite-container substitution has been made.

For the Beltrami base on |q|=|lambda|, the same energy can be written

    e(F)=rho/2 integral |F omega_hat(q)|²/|F^-T q|² dmu(q).

After the declared Haar mixture the real transverse covariance on direction
n is proportional to `I-n n^T`; its imaginary helical part drops out of
this real symmetric energy. If `E[u0 tensor u0]=U_*² I`, then with
`C=F^T F`, `e0=3 rho U_*²/2`,

    e(F)/e0=average_n [(tr C-n.C n)/2]/(n.C^-1 n).

Take F=exp(t E), E symmetric tracefree. With a=n.E n, b=n.E² n and
T=tr E², the ratio before averaging is

    1+t a+t²(T-3 b+2 a²)+O(t³).

The exact sphere moments give `<a>=0`, `<b>=T/3`, `<a²>=2T/15`, hence

    e(exp(tE))=e0+(2 rho U_*²/5) t² ||E||_F²+O(t³),
    mu_affine=2 rho U_*²/5>0.

For a rotation F, C=I and the energy is exactly e0. This is the same
random EPS-containing field and its actual affine Euler energy, not the
periodic comparison coefficient in 0043.

## Full compact reaction, not an averaged zero-cross shortcut

Compose this affine map with the compact stationary orbit family. Its
exact joint energy is `rho E|M_F u_q|²/2`. Let p denote all compact
conjugate amplitudes, with their full positive stationary momentum
operator P from stationary-assembly.md. The mixed affine/conjugate Hessian
N is retained before any independent realization's reaction is eliminated.
The corrected shear matrix and isotropic coefficient are

    C_eff=C_affine-N^* P^-1 N,
    mu_eff=mu_affine-tr_STF(N^* P^-1 N)/10.

The trace includes the ensemble/Palm inner product in P and N. The mean
of a cross may vanish by isotropy while the mean of its squared Schur
correction does not. The formula above retains that nonnegative reduction.

It has a finite explicit bound. For tracefree h and n=q/|q|,

    T_h u := d/dt M_exp(th) u at t=0
           =-h^T u+n[n.(h+h^T)u], n.u=0,
    ||T_h|| <=3 ||h||_F.

Consequently

    N(h,p)=rho E[v_p.(T_h+T_h^*)u0],
    |N(h,p)|<=6 rho ||h||_F sqrt(3) U_* ||v_p||_stationary.

If `||v_p||<=C_F ||p||_Palm` (the bounded good-patch coefficients and
stationary Leray contraction provide C_F), then
`||N||<=N_*=6 sqrt(3) rho U_* C_F`, independent of carrier magnitude.
With the previously derived coercive bound `P>=p_min I`,

    mu_eff >=mu_affine-N_*²/(2 p_min)>0

at a finite carrier threshold. Here p_min grows linearly with the carrier
magnitude, whereas N_* is fixed. Every number is a physical input,
actual-field integral or explicit norm bound, not a fit to the desired
modulus. Including the compact angle coordinates in the positive block
gives the analogous positivity proof for the full joint affine/angle
elastic matrix. Local spin/shear crosses vanish only after the complete
reduced ensemble action is formed and its irreducible representations
are averaged.

This is a simpler positive shear route than adding nine auxiliary affine
cages. The latter remains valid branch evidence; its coefficient is not
silently substituted for this newly reconstructed spectral one.

## Intensity-aware cross bound and the distinct removal limits

On transverse inputs the mixed formula simplifies further:

    N(E,p)=-2 rho E[v_p.E u0]
           =-2 rho E[F_p.P0(E u0)].

Since u0 is monochromatic, its projected affine response is the LOCAL
finite-derivative expression

    P0(E u0)=E u0+lambda^-2 grad div(E u0).

It is uniformly bounded by `C_u ||E||_F` on a good patch using the selected
u0 and second-derivative bounds. If
`integral_patch |F_p|<=C_F1 |p_a|`, the exact intensity identity yields

    |N(E,p)|<=2 rho n C_u C_F1 ||E||_F E_Palm |p_a|
             <=2 rho sqrt(n) C_u C_F1 ||E||_F ||p||_Palm,n,
    ||p||_Palm,n²=n E_Palm |p_a|².

Thus the shear Schur correction is O(n/p_min), not an intensity-independent
background modulus. This refined bound is again derived before averaging
the reaction square.

Removing selected core/cage coherence (`n->0`, or its explicitly defined
selected length density `L_selected=n E_Palm length`) at fixed Gaussian
amplitude U_* removes the angle/curvature sector and this reaction
correction, but leaves `mu_affine=2 rho U_*²/5`. Its remaining continuum
is the divergence-free Navier--Cauchy background, not an empty Euler fluid.
The Gaussian background vorticity has NOT disappeared in that limit.

Removing the entire vortical background additionally requires `U_*->0`.
Then the affine shear energy vanishes as well. At exactly U_*=0 the KKS
form loses rank, so that limit is taken on the unreduced material-fluid
action; a quotient such as B²/H is not retained as a spurious spin inertia
after its coordinate chart has degenerated. In both distinct operations
the total material density remains rho. No single ambiguous L_v symbol
is used for both removals.
