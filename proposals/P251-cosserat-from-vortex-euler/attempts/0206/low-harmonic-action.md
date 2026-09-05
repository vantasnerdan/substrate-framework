# Actual rigid-motion branch and the positive nonrigid bending form

The root activated this attempt after central schema success with 268
accepted claims. The sources opened afterward were the exact 0186
toroidal kernel/full-poloidal operator and the 0195 global force-free
family/full current. All conclusions below distinguish an actual Euler
solution from a quadratic form restricted to geometric tangents.

## 1. An exact Euler translation/orientation response, not a fake gap

Let v_*(x-Ut e_z) be the laboratory velocity of the actual isolated
traveling ring. Its relative steady field is u_*=v_*-Ue_z. For a fixed
orthogonal rotation Q and initial center a,

    v_Q,a(t,x)=Q v_*(Q^T(x-a)-Ut e_z)              (1)

is an EXACT Euler solution. In the original translating frame its
center is a+Ut(Qe_z-e_z), and its geometric axis is N=Qe_z.
Differentiation in a transverse small rotation epsilon gives

    delta X=a+Ut delta N,  delta N=epsilon cross e_z.               (2)

The actual covariance tensor of any full invariant tube fraction is
axisymmetric at the base. Its perpendicular/parallel eigenvalues are
different for a thin ring, so the principal normal is an actual smooth
Euclidean orientation observable with unit rigid-rotation response.
Rotating the actual material tube and its actual velocity in (1) gives
precisely (2), including the moving domain and all ambient fluid.

The translation momentum on the Euler orbit is the Kelvin impulse

    I=(rho/2) integral x cross omega dx=I0 N, I0>0.

It is not silently equated to the finite tube's mass times velocity.
Its KKS moment identity pairs a translation with a rotation as

    Omega(translation(a), rotation(epsilon))
       =a dot delta I=I0 a dot(epsilon cross e_z),                 (3)

up to the single consistent sign convention for the symplectic form.
The rotation/rotation pairing ALSO retains the actual axial angular
impulse J0: Omega(K_x,K_y)=J0. For this axisymmetric compact-swirl field
J0 is the complete axial swirl moment; it need not equal the spin of
only a selected inner tube. On the Euclidean orbit the laboratory
energy is rotation invariant.
The ACTUAL translating-frame Hamiltonian is E-U I_z. Its transverse
second variation is U I0 |delta N|²/2. Writing p=I0 delta N therefore
gives the exact quadratic Euclidean suborbit action

    p dot X_t + J0/(2I0²)(p_x p_{y,t}-p_y p_{x,t})
                  -|p|²/(2M_eff), M_eff=I0/U>0.   (4)

The retained momentum-space gyro term is essential to the full phase
form, although it drops from (2) because p_t=0 on these histories.
Its equations therefore recover (2). This is an actual positive translational
phase metric supplied by the same complete Euler fluid, with its
ambient impulse reaction. No rigid-body mass has been appended.
The orientation is constant on these histories; the zero eigenvalue
is a translation/tilt Jordan pair, not a positive optical gap.

The literal tube mass M_D gives P_D=M_D U N in the laboratory frame,
whereas I0 need not equal M_D U. Likewise its base intrinsic spin is
S_D=S0 N for the swirling ring, and delta S_D=S0 delta N. This is a
rotated pre-existing spin, not a newly generated S=j epsilon_t. Under
the registered time reversal that static spin reverses sign. These
physical rows are kept separately from the canonical impulse in (3).

## 2. The same-fluid leading nonrigid bending energy

Parameterize the perturbed centerline by

    X(phi)=(R+q(phi))e_r(phi)+z(phi)e_z.

Use smooth low fixed harmonics and actual volume-preserving tubular
deformations of the base vorticity. A tubular coordinate map transports
the cross-section and corrects its area by its longitudinal Jacobian;
a smooth local volume correction completes it without changing the
centerline. Thus these are genuine Euler coadjoint tangents with their
full inverse-curl velocity, not a postulated filament modulus.

The complete Biot--Savart energy has the local logarithmic part

    E_log=Tau_R length(X),
    Tau_R=rho Gamma² log(R/a)/(4pi).               (5)

Here a is the fixed microscopic core scale; changing a by a fixed factor
changes only the bounded finite part. To derive (5), split the exact
double Green integral into |s-s'|<a, a<|s-s'|<cR and its complement.
In the intermediate range the integrated axial vorticity is the actual
circulation Gamma and the leading kernel is Gamma²/|s-s'|. Integrating
the two sides supplies its logarithm and the factor in (5). The compact
poloidal vorticity has zero cross-sectional vector integral, so it has
no extra logarithmic monopole. Its energy and all nonlocal return terms
belong to the finite part; they are not set to zero.

On a fixed smooth dimensionless neighborhood of the circular centerline,
subtracting this explicit local singularity leaves a C² shape functional
with bounded second derivative. Near the diagonal, subtract tangent and
Jacobian at the same arclength before differentiating: each remainder
gains a power of separation, making its first two shape derivatives
integrable. In the finite core, rescale the two cross-section variables
by the FIXED smooth profile. In the far part use their positive geometric
separation. Restoring physical amplitude q,z gives an O(rho Gamma²/R)
bound for that finite second form on each fixed harmonic space. The
bound includes the actual swirl, volume return and profile deformation.
This is a C² energy estimate, not an asserted dynamical spectral limit.

The same full source moment gives

    I_z=(rho Gamma/2) integral (R+q)² dphi
                      + finite-core moment corrections.          (6)

Those corrections have no logarithmic enhancement. The actual speed
U=Gamma log(R/a)/(4pi R)+O(1/R) from 0195 therefore makes the first
variation of E-U I_z vanish, and its leading SECOND variation is

    H_bend=Tau_R/(2R) integral(q_phi²+z_phi²-q²)dphi
                      +O(rho Gamma²/R)||q,z||².   (7)

The -q² term is the actual translating-frame impulse reaction. Omitting
it would give a false positive restoring force for rigid translation.
The leading KKS form follows by integrating actual vorticity over each
cross-section:

    Omega_bend=-rho Gamma R integral dq wedge dz
                          +O(rho Gamma a)||q,z||².                (8)

For q=Q cos(n phi), z=Z cos(n phi), the leading Hessian and symplectic
coefficient are

    H=pi Tau_R diag(n²-1,n²)/R,
    B=-pi rho Gamma R.                            (9)

For every fixed n>=2 the complete restricted finite-core H is positive
and B nonzero at sufficiently large finite R: its positive logarithmic
margin dominates the derived bounded finite part. The corresponding
geometric TWO-COORDINATE action would have

    omega_bend²=[Gamma log(R/a)/(4pi R²)]² n²(n²-1).                (10)

Equation (10) is the leading restricted-action frequency, not yet a
full Euler pole. The next companion retains the exact complement before
deciding when these two coordinates follow it.

For n=1, (9) is semidefinite. Its zero radial translation coordinate
and conjugate tilted plane reproduce the exact Jordan structure (2).
No n=1 negative stability conclusion or optical gap is inferred.

## 3. Exact global tensor selection

For any axisymmetric invariant fraction, a physical cylindrical-vector
perturbation proportional to exp(i n phi) has a global vector moment
only at n=0,1: transforming its components into Cartesian vectors adds
only frequencies 0,+1,-1. This applies to the COMPLETE linear spin
integrand, including moving-domain, centroid and background-velocity
terms, since those expressions are equivariant vectors. It is not
specific to a chosen local profile.

An actual rank-two covariance has harmonics through |n|=2. Its normal
tilt uses n=1 because Q_xz,Q_yz transform as planar vectors. Its planar
anisotropy Q_xx-Q_yy and 2Q_xy uses n=2. Hence the positive n=2 bending
form has a real global quadrupole response, but zero global linear
spin vector for an axisymmetric whole-tube fraction.

The circular background has Q_xx=Q_yy: an in-plane angle obtained by
dividing by that zero gap is not a linear physical orientation. The
nonrigid quadrupole components themselves remain valid observables.
For instance, at the leading thin-ring order and uniform material
line measure M_D/(2pi), q=Q cos(2phi) yields

    delta(Q_xx-Q_yy)=M_D R Q,
    delta Q_xy=0.                                 (11)

The actual finite-core tensor rows have their explicitly small geometric
corrections. Introducing an actual nonaxisymmetric stationary structure
or a transported nonnegative tag can create a nonzero angle gap; it
also changes the current/transport problem and is not done by silently
renaming (11) as microrotation.

## 4. Primary comparison, with its actual license

[Jerrard--Seis, *On the vortex filament conjecture for Euler flows*,
Theorems 1--2](https://arxiv.org/pdf/1603.00227) gives a quantified
instantaneous Euler/filament correspondence and a finite-time conclusion
conditional on persistence of vorticity concentration and an energy
excess bound. It does not supply that persistence for this nonaxisymmetric
swirling-ring perturbation. Its Hamiltonian/Poisson viewpoint supports
the representation used here, but is not imported as an unconditional
mode or optical-window theorem. The exact same-field rigid histories
in Section 1 require no filament approximation at all.
