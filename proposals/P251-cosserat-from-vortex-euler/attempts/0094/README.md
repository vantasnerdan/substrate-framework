# 0094 — positive shear from the actual material action

Owner `/root`. Parent P251 remains active. Frozen continuation of0090/0091:
derive affine shear from the NEW material Jacobi action rather than copying
the old vorticity-pushforward coefficient. Candidate A is a separate
strain-attached positive material cage with no independent angle; candidate B
is full Kelvin/relabeling reduction of the affine background. This record
executes A at the conditional material-action scope. No empirical comparator,
fitted modulus, or all-wave-number Euler-invariant ansatz is introduced.

## Exact input and the uncorrected background

Use the actual stationary isotropic mean-zero Euler field law with
E[u_i u_j]=U_*^2 delta_ij and the material Jacobi functional of0084.
For a deterministic slow divergence-free displacement U(x), stationarity
gives E[Hess p]=0. Thus the uncorrected material stiffness is

    K_bg(U)=-rho U_*^2 integral |grad U|^2.

This is NOT the frozen-vorticity affine coefficient of0057. For compactly
supported divergence-free U, integration by parts gives
integral |grad U|^2=2 integral |sym grad U|^2, so its shear coefficient
in the potential K/2 is -rho U_*^2. This negative contribution is retained
in the calculation, not renamed as a positive modulus.

Use a positive compact material cage X on the SAME field from0090,
with k_X=K(X)>0 and m_X=rho integral |X|^2>0. Sample its complete field and
finite geometry through a stationary isotropic marked ensemble. Select a
separate population of disjoint cage supports with finite positive intensity
nu_s. The actual good-patch reconstruction includes uniform positive lower
bounds and bounded norms; it is the same support/Poisson construction as0057,
with these material rather than coadjoint functionals. Material marks are
held fixed under variations.

## Strain attachment and its derived coefficient

Let T0 be a unit symmetric tracefree tensor and let T=R T0 R^T be its
proper-rotation mark, with the background and cage rotated simultaneously.
Define the scalar attached amplitude as a=t_s T:E, E=sym grad U, and the
actual microscopic material displacement as a X at each cage. This is a
declared affine kinematic attachment, not a prescribed elastic energy.
The irreducible five-dimensional tracefree representation gives exactly

    E_R[(T:E)^2]=|dev E|^2/5.

Consequently the self term in the averaged potential is

    Delta W_self=nu_s t_s^2 E_Palm[k_X] |dev E|^2/10.

The full cross with the bare macro lift is retained. Its translation part
vanishes exactly: for any constant c and compact divergence-free X,

    K(c,X)=integral c_i p_ij X_j=0,

by integrating in j and using div X=0. The affine cross is finite and
linear in t_s. Isotropy makes its symmetric-tracefree restriction a scalar
multiple of E:E. Let b_s denote that actual signed scalar (computed from
the full field/cage functional), not its absolute-value bound. Then

    mu_total=-rho U_*^2+b_s t_s+a_s t_s^2,
    a_s=nu_s E_Palm[k_X]/10>0.

For any fixed b_s and finite lower-order contributions from the remaining
fixed geometry, this quadratic is positive at a finite t_s. For example
if c_s bounds the absolute sum of those nonquadratic shear terms, any
t_s> (|b_s|+sqrt(b_s^2+4a_s c_s))/(2a_s) suffices. These are analytic
sign conditions on a declared material geometry, not comparison to target
elastic data. All coefficients are actual action integrals.

## Keep the added inertia and physical currents

The same substitution adds

    Delta T_self=nu_s t_s^2 E_Palm[m_X] |dev Edot|^2/10.

It is a strain-rate-gradient inertia term and remains in the second-gradient
mass matrix. It is not the unchanged-kinetic construction of0065. The
compact divergence-free identity integral X=0 removes the uniform-translation
cross. Its first moment is antisymmetric, so contraction with symmetric
Edot removes that bare affine kinetic cross. Other retained mean/current
and spin-coordinate terms follow the SAME full action and are not dropped
by this observation. An isotropic strain-rate-to-axial-spin linear response
vanishes by tensor type; this does not erase the cage's microscopic spin.

For the acoustic branch omega^2=O(k^2), this additional O(k^2) mass first
changes dispersion at O(k^4), not the leading shear wave speed. For an
optical branch whose macro translation amplitude is O(k), its pure U-block
gradient mass similarly enters only beyond the retained optical k^2 term.
Both statements require the full coupled mode scaling; they are not a
license to delete this mass from the action or observation map.

## Distinct removal limits and scope

Keep this affine-only population while removing the separate population of
independently rotating coherent tubes. The resulting positive continuum has
no independent microrotation and is in the incompressible Navier–Cauchy
sector with the ACTUALLY derived mu_total. This is a sector/field-limit
statement, not equality to the old five-scalar filament modulus. If all
vortical background and attached material structure are removed together,
their stiffnesses vanish and the neutral incompressible Euler limit remains.
The exact population/scaling definition must be part of the parent contract;
there is no suggestion that deleting all material correlations at fixed
U_* preserves a positive shear coefficient.

Established conditional proposition: positive material affine shear can be
constructed with its entire same-action inertia and cross terms retained.
This does not resolve the remaining0091 mean/core-frame/Kelvin reduction.
The parent continuum action, physical interpretation and promotion remain
active. No accepted claim or earlier attempt is rewritten by this route.
