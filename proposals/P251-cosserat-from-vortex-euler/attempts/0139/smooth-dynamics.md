# Smooth stationary compensated cores and controlled actual Euler histories

This is a constructive transfer, not a claim that steady desingularization
alone supplies finite-core eigenmodes. All constants below refer to a
fixed finite supercell, fixed nonzero acoustic k and fixed finite time T.
The small core radius is chosen after those quantities.

## 1. Exact C-infinity stationary base

On the primitive triangular torus of area A use Delta G=delta-1/A and
G=N+R near zero, N=log|x|/(2pi). The regular part R is C-infinity and
sixfold invariant. Use the smooth monotone flat core F_eta,U_0 of0036,
with integral F_eta(U_0)=Gamma; eta here is fixed once, not the Gaussian
point-energy regularizer in triangular-action.md.

On a fixed disk B_L, restrict U to sixfold invariant functions and solve

    U-mu+N*F_eta(U)
       +integral [R(epsilon(y-z))-R(0)]F_eta(U(z)) dz=0,
    integral F_eta(U)=Gamma.                                    (1)

At epsilon=0 this is the radial Green equation. The only bounded kernel
without the mass border consists of translations, by the explicit0036
radial gap proof; sixfold invariance excludes them. The mass border
removes the radial logarithmic direction exactly as in0036. Therefore
the bordered derivative in (U,mu) is invertible. Ordinary smooth IFT
gives U_epsilon=U_0+O(epsilon²) in C^(1,alpha), since grad R(0)=0.
The flat negative collar persists and elliptic bootstrap gives every
fixed finite regularity norm. Define the *total* physical vorticity

    omega_epsilon(x)=epsilon^-2 F_eta(U_epsilon(x/epsilon))-Gamma/A

on the core, extending the first term by zero elsewhere. It is globally
C-infinity. Its excess circulation is exactly Gamma, so the total mean
vorticity is zero. Its streamfunction is periodic and (1) says the excess
is a function of psi on the core. Outside, total vorticity is constant.
Thus u_epsilon=J grad psi is an exact smooth stationary Euler field.
Copying it into the finite supercell gives N identical stationary cores.

This construction fixes actual profiles and all their Casimirs before
the dynamical preparation. It need not vary the steady profile when
constructing a wave: each core is subsequently transported by an actual
area-preserving map.

## 2. Exact Kelvin-prepared perturbations and their complete action

Choose a real supercell Bloch mode with nonzero k and zero sum of centre
displacements. For small amplitude delta, translate each core by delta
xi_j using a smooth Hamiltonian material map equal to that translation
on a neighborhood of its support and smoothly returned in disjoint
collars. Such a map is the time-one flow of a compactly supported stream
function that is linear on the core neighborhood. These maps preserve
area, every vorticity Casimir and the uniform compensating background.

Use the actual coadjoint push-forward of the stationary velocity one-form,
including its torus harmonic component, to set the initial velocity.
Because the total circulation-weighted core displacement is zero and
there is no winding, the induced harmonic impulse is zero. Equivalently,
integrating the infinitesimal mean Kelvin variation along these local
translations gives J sum_j Gamma delta z_j/(NA)=0. Thus the preparation
has the fixed zero mean used by the point Green model; it is not imposed
by resetting a nonzero Kelvin mean afterward.

The KKS pairing of two such core translations is exactly

    Omega=rho Gamma sum_j det(xi_j,eta_j).                       (2)

The uniform part contributes a constant times the integral of a Poisson
bracket and hence vanishes. Inside each core the generators are constant,
so its contribution is exactly its circulation times the determinant,
independent of finite core shape. The full kinetic Hessian along these
maps retains every ambient velocity and delta²omega term. Self-energy
is translation independent, while distinct-core interactions converge
in every fixed centre derivative to the torus Green interactions as
epsilon ->0. Indeed, Taylor expansion against the exact centred profiles
gives O(epsilon²) on any fixed noncollision configuration neighborhood.
Consequently its finite-dimensional action converges to the same positive
point-lattice block in triangular-action.md, with a controlled finite-core
error smaller than the k² stiffness once epsilon is sufficiently small.

This action restriction alone is not used to assert Euler reconstruction.
The following actual evolution estimate supplies that missing step for
the specified finite histories.

## 3. Source-licensed localization and the periodic transfer

Marchioro, On the Localization of the Vortices, BUMI1-B (1998),571--584,
http://www.bdim.eu/item?fmt=pdf&id=BUMI_1998_8_1B_3_571_0,
Theorems1.1 and2.1, proves support localization of positive concentrated
vorticity blobs on every fixed pre-collision interval. For each exponent
0<alpha<1/3 and initial support radius epsilon with an L-infinity bound
M epsilon^-gamma0, the support stays within C_(alpha,T)epsilon^alpha of
the corresponding point path. The theorem's proof first allows an
arbitrary bounded, divergence-free, Lipschitz external velocity field.

Here it transfers directly, with an explicit local decomposition. In a
disk around the ith blob, the periodic kernel is the planar singular
Biot--Savart kernel plus J grad R. The latter and all interactions with
other separated blobs are bounded Lipschitz external fields, uniformly
in epsilon. The negative uniform background is exactly unchanged by
transport and is already included in R; it introduces no new fluctuating
source or singular kernel. Cutoffs extend the local smooth external
streamfunction to a bounded Lipschitz field if necessary.

The source proof uses: cancellation of the singular self-induced centroid
velocity; second-moment inequality I_i'<=2L I_i; the radial kernel identity
x dot K(x)=0; and a dyadic cutoff estimate for exterior blob mass. Each
survives unchanged under that external-field decomposition. In particular
the outer-mass recursion has coefficient

    C[epsilon/r³+epsilon²/r^4+1],

and its factorial iteration makes the outer mass smaller than every power
when r>=C epsilon^alpha, alpha<1/3. The remaining radial velocity obeys
|r'|<=Cr+C epsilon/r²+superalgebraic remainder. Gronwall gives the stated
localization. A bootstrap using the fixed minimal point separation keeps
these local disks disjoint throughout T. Thus the proof applies on this
finite flat torus, rather than importing a whole-plane theorem across an
unexamined pressure boundary.

The exact Euler excess blobs q_i=omega_i+its labelled constant subtraction
are most simply defined as transported positive initial core densities;
omega=sum_i q_i-Gamma/A. They remain positive, with conserved mass Gamma.
Their centres Z_i^epsilon and the point paths z_i obey

    sup_[0,T] |Z_i^epsilon-z_i| <= C_T epsilon,
    support q_i subset B(z_i,C_(alpha,T)epsilon^alpha).           (3)

The centroid estimate follows from the same antisymmetry and Lipschitz
Gronwall proof; the smooth torus self-remainder adds only O(epsilon).
All smooth solutions exist for the required time by two-dimensional Euler
well-posedness, and the corresponding three-dimensional z-independent
extension is an exact solution for that same time.

## 4. Full physical coarse velocity, acoustic time and ordered errors

Use the actual k-Fourier component of the full Euler velocity. For k!=0,
the constant background has no Fourier coefficient and inverse curl gives

    uhat(k,t)=-i Jk/(|k|² NA) sum_i integral q_i(x,t)e^-ikx dx.

Comparing each exponential with e^-ik z_i(t) and using (3) yields

    |uhat_Euler(k,t)-uhat_point(k,t)|<=C_(alpha,T)epsilon^alpha.  (4)

This is a complete periodic pressure observable; it retains the ambient
flow and does not replace fluid momentum by a vorticity-weighted centroid.
Physical transverse material displacement is its time integral with a
declared initial value. Its error is at most T times (4).

Choose a fixed sufficiently small nonzero admissible k. The point-lattice
linear mode has frequency omega=c_T|k|+O(|k|³), positive action and actual
mean velocity amplitude of order delta|k|. Choose T equal to a fixed
number of those periods. Smooth finite-dimensional point dynamics near
the separated equilibrium follows its linearized solution with error
C_T delta² over T; its exact nonlinear mean-velocity observation has the
same Taylor control. Now choose delta small, then epsilon still smaller,
so that

    C_T delta² + C_(alpha,T)epsilon^alpha << delta|k|,
    epsilon² << |k|².                                          (5)

These are finite, independently selectable geometric/amplitude scales,
not a fitted sound speed. They give genuine smooth Euler histories with
nonzero transverse coarse velocity and the positive acoustic phase/action
normalization through a full prescribed acoustic observation interval.
Their initial data are exact Kelvin preparations of the constructed smooth
stationary periodic base. All mean and ambient contributions are retained.

This is a controlled finite-time/ordered-limit realization. It does not
prove an exact Bloch eigenbranch for one fixed finite core radius, a
uniform joint k,epsilon ->0 theorem without (5), or a complete isotropic
three-dimensional optical/acoustic Cosserat continuum. Those are distinct
joining obligations; the present construction supplies actual smooth Euler
transverse acoustic content, beyond a projected oscillator.
