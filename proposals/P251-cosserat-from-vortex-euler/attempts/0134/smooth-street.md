# Executed smooth-core continuation and its dynamical boundary

We now construct a genuine C-infinity finite-core traveling Euler street,
not merely a patch analogue or a point-vortex assertion. The proof uses
the radial bordered inverse explicitly proved in attempt0036,
smooth-polygon-construction.md sections1--2 and radial-core-gap.md.
The signed, periodic geometry below is new; it is not inherited from the
same-sign polygon. This proof asserts steady existence and physical
impulse response. It does not identify the steady family with a single
Kelvin leaf or transfer a point-vortex spectral theorem.

## 1. Exact periodic equation

Write J(x,y)=(-y,x), u=J grad psi. The period-a Green function is

    G_a(x,y)=1/(4pi) log[2cosh(2pi y/a)-2cos(2pi x/a)].

It has Delta G_a equal to the periodic point source. Near zero,
G_a=N+R_a, N=log|x|/(2pi), and R_a is smooth. Its linear growth at
infinity is harmless for a neutral pair; their velocities decay there.

Fix a,b,Gamma>0, centres c_+=(0,b/2), c_-=(a/2,-b/2), and reflection
S=diag(1,-1). Use the positive smooth radial core U_0 and flat monotone
F_eta from0036, with integral F_eta(U_0)=Gamma. On a fixed disk B_L,
seek U even in y_1, with exact mass Gamma and vertical centroid zero.
Set

    omega(c_++epsilon y)=epsilon^-2 F_eta(U(y)),
    omega(c_-+epsilon S y)=-epsilon^-2 F_eta(U(y)),

and extend periodically and by zero outside the disjoint disks. Put
D=c_+-c_-=(-a/2,b). The following equation includes *all* periodic and
mutual interactions:

    U(y)-mu+N*F_eta(U)(y)
      + integral [R_a(epsilon(y-z))-R_a(0)]F_eta(U(z)) dz
      - integral [G_a(D+epsilon(y-Sz))-G_a(D)]F_eta(U(z)) dz
      + epsilon V y_2 + lambda y_2 =0,
    integral F_eta(U)=Gamma,
    integral y_2 F_eta(U)=0.

At epsilon=0 its derivative in (U,mu,lambda) is precisely the0036
bordered inverse, with the two coordinate axes interchanged. Evenness in
y_1 removes the horizontal translation kernel; the vertical translation
is removed by the centroid border. The speed V is initially a parameter.
All external kernels are smooth for sufficiently small epsilon, all
composition maps are smooth on C^(1,alpha), and the negative collar
persists. The implicit function theorem therefore solves (U,mu,lambda)
smoothly in (epsilon,V,a,b), on a neighborhood of each fixed a,b>0.
Reflection y_1 -> -y_1 preserves the periodic equation because D_1=-a/2.

The first-order external gradient is -Gamma grad G_a(D), with

    partial_1 G_a(D)=0,
    Gamma partial_2 G_a(D)=Gamma tanh(pi b/a)/(2a)=V_0(a,b).

Thus U has no first-order variation and

    lambda(epsilon,V)=epsilon(V_0-V)+O(epsilon^2).

The quotient lambda/epsilon extends smoothly to zero by Hadamard's
formula. Its V derivative is -1 there. A second implicit function
theorem gives an exact speed V_epsilon(a,b) with lambda=0. In particular

    V_epsilon -> V_0 in C^1 on each fixed small parameter neighborhood,
    partial_b V_epsilon >0 for sufficiently small positive epsilon.

No external force or prescribed speed remains. Elliptic bootstrapping
and the flat negative collar give globally C-infinity vorticity.

## 2. Actual Euler and the attempted global scalar relation

Let psi=G_a*omega with its glide-antisymmetric normalization. The two
cores and reflection give

    psi(x+a/2,-y)=-psi(x,y).

At lambda=0, the equation is exactly U=mu_phys-[psi+Vy] in the positive
core, after absorbing the subtracted constants and logarithmic self
constant into mu_phys. The negative core gives U=mu_phys+[psi+Vy].
For sufficiently small epsilon, mu_phys=(Gamma/(2pi))log epsilon+O(1)<0.
Writing s=psi+Vy, the tempting global extension would be

    omega=f_epsilon(s),
    f_epsilon(s)=epsilon^-2[F_eta(mu_phys-s)-F_eta(mu_phys+s)].

But that extension fails: outside the supports s is harmonic, and it
tends to +/-infinity as y tends to +/-infinity when V>0. Both ends cross
the active core threshold values, so **the displayed globally
single-valued formula is not licensed there**.
Only the two-core local relation is needed for exact Euler: vorticity is
zero in the exterior and flat across its boundaries. Consequently

    (J grad psi-V e_x).grad omega
       =J grad(psi+Vy).grad omega=0

everywhere, and omega(x-Vt e_x) is an actual smooth Euler solution.

The far-field streamfunction issue is substantive for an attempted
global generalized-force-free lift: the moving-frame scalar s is
unbounded although f_epsilon is nonzero on some of its values inside
the cores. The simple global scalar formula would create spurious
far-field vorticity. A separately proved exterior branch/first-integral
construction is needed; a local scalar relation on each disjoint core
is insufficient for the parent's0136 global single-valued lift.

## 3. Finite-core physical momentum and what still needs construction

Neutrality makes u decay exponentially at large |y|. Periodicity and
integration by parts give the exact full-fluid momentum per period

    P_x=rho integral u_x=rho integral y omega=rho Gamma b.

The two mass and centroid constraints establish the final equality for
every constructed finite epsilon, not just in the point limit. Hence
the actual steady-family response is

    dP_x/dV_epsilon=rho Gamma/(partial_b V_epsilon)>0,
    dP_x/dV_epsilon -> 2rho a^2 cosh^2(pi b/a)/pi.

This is a useful positive **smooth-Euler impulse/speed response**.
The family U_epsilon(b) may change its full vorticity distribution, even
though its total circulation is fixed. The semilinear construction does
not fix every Casimir. Therefore this derivative is not asserted to be
the exact dynamic inertia obtained by eliminating a same-Kelvin shape
momentum. Likewise, steady desingularization alone does not transfer the
long-wave lattice spectrum to finite cores. An isovortical core branch
and its full Euler response operator are the next dynamical construction.

This precise boundary avoids a false promotion while preserving an
actual finite-core stationary Euler object and its measurable positive
momentum response. It does not close the stronger parent continuum.

## 4. Executed Bernoulli repair of the failed global-F route

Parent0136 supplied a materially different global construction that does
not require a single-valued vorticity/streamfunction relation. Apply it
directly to the smooth stationary moving-frame velocity v=u-Ve_x above.
Its globally defined periodic pressure p exists because the closed steady
Euler pressure one-form has zero period on a far-field horizontal loop.
Let B=p+|v|^2/2 (kinematic pressure) and choose a finite C>sup B. Both B
and its derivatives are bounded; B tends to the same constant at the two
far ends. Set

    W=sqrt(2(C-B)),  u_3=(v_x,v_y,W),  partial_z u_3=0.

Then v.grad B=0, so v.grad W=0; the full three-dimensional Euler equation
is exactly the planar equation plus a zero axial acceleration. Moreover
its full Bernoulli function is constant. With our convention v=J grad s,
the exact planar identity is grad B=omega grad s, and therefore

    grad W=-(omega/W)grad s,
    curl u_3=(omega/W)u_3.

This is a genuine global smooth generalized force-free Euler street.
Outside the compact cores omega=0, hence grad W=0, and the flat core
collars join smoothly. It requires neither magnetic dynamics nor a
global F(s); it repairs precisely that route's identified obstruction.
The z-independent planar Euler subsystem, including its pressure and
mean response, is exact: the axial W is merely transported under a
planar area-preserving material map. Its integrated kinetic contribution
is constant on that class of maps, not a hidden planar restoring term.

It remains a street extruded along an infinite axis, not a closed EPS
knot or a periodic isotropic continuum. The finite-core dynamical and
nonzero-K complete-momentum boundaries in sections3 and street-action.md
are unchanged by this exact three-dimensional stationary lift.
