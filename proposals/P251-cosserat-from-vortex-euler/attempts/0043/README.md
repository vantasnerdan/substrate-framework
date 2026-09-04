# 0043 — affine translational energy of the same Beltrami fluid

Parent: P251 N3/N4, exact conditional affine Euler coarse-graining. This is
an append-only representation change motivated by 0040's positive internal
sector and its still-missing nonuniform translational action. Main owns this
directory. The object is the exact Biot–Savart energy of a volume-preserving
affine pushforward of the stationary Beltrami field, including its advected
vorticity and deformed periodic lattice. No independent elastic potential is
introduced. The material-action mass from 0037 is a distinct input, not
inferred from circulation-centroid dynamics.

For F in SL(3), pushforward sends Fourier vorticity omega_k to F omega_k
and its wave covector k to F^{-T} k. Incompressibility is preserved exactly.
The cell energy follows by inverting curl in the deformed lattice. We will
derive its full quadratic coefficient under F=exp(E), E symmetric traceless,
then average a declared isotropic distribution of actual tube-cell frames.
Each background realization remains an exact stationary Euler field before
affine perturbation. Their ensemble is an expectation of energies, not an
unsupported superposition of independently stationary fields.

Selection criteria frozen before evaluation: actual Euler pushforward,
objectivity, no fitted coefficient, nonzero positive shear response, and
explicit finite-cell/long-wave boundary. Compare a circular-helicity shell
identity with direct Fourier evaluation of 0040's two-mode tube. This is a
fixed algebraic target with no empirical comparator or numerical remainder.
An affine energy alone licenses only the corresponding constrained elastic
coefficient, not a complete coupled kinetic action or EPS knot embedding.

Status: active. Exact symbolic derivation is the strongest practical oracle;
the incorrect fixed-wavevector transformation will be checked as a mutation.

## Exact affine energy and isotropic coefficient

Let omega_k be the Fourier vorticity of a periodic zero-mean Euler field.
The deformed lattice is part of the material affine map, so

    k_F=F^{-T}k, omega_{F,k}=F omega_k,
    u_{F,k}=i k_F cross omega_{F,k}/|k_F|²,
    E(F)/volume = rho/2 sum_k |F omega_k|²/|F^{-T}k|².

The last equality uses k_F.omega_{F,k}=k.omega_k=0 and det F=1.
These are the physical curl inverse and Euler kinetic energy, not the
kinetic energy of a merely transported velocity vector. The latter would
give a different coefficient. Mean momentum can be restored independently
with its exact Galilean mass term as in 0037/0040.

Each Fourier pair of a real constant-eigenvalue Beltrami field is circularly
polarized. If n=k/|k| and C=F^T F, its energy ratio is exactly

    R(F,n) = [tr C - n.C.n]/[2 n.C^{-1}.n].

For F=exp(t E), E=E^T, tr E=0, define s=n.E.n, u=n.E².n,
and T=tr E². Direct expansion, including the reciprocal denominator, gives

    R=1+t s+t²(2s²+T-3u)+O(t³).

The declared isotropic distribution is the rotation of the entire actual
cell and all its generators together, not independent randomization of
the internal locking angle. Thus every shell direction is uniformly
distributed while the core/cage correlation remains. The existing exact
sphere moments give <s>=0, <u>=T/3, <s²>=2T/15. Consequently

    <E(F)>/volume = E0 + (4 E0/15) tr E² + O(E³),
    mu_affine = 4 E0/15 > 0.

Here E0 is the actual unstrained Euler energy density. For the exact tube
u0=(-b sin y,a sin x,a cos x+b cos y), E0=rho(a²+b²)/2, hence

    mu_affine = 2 rho(a²+b²)/15.

This scalar has pressure units and contains no fitted elastic constant.
At the geometry a=2b chosen in 0040 it is 2 rho b²/3. This positive shear
is a coefficient of the same fluid ensemble that supplies the internal
angle/shape action, rather than a separately assumed P242 filament spring.
The calculation is an exact quadratic affine energy identity for a smooth
finite-scale field; no thin-core limit or singular cutoff enters it.

Common spatial rotation sends F to QF and rotates both omega and wave
covectors together, leaving the energy exactly unchanged. The shear
coefficient therefore does not assign a spurious spring to free rigid
rotation. Under incompressibility, volumetric strain is not an admissible
variation: its multiplier is pressure. This result does not claim an
independent compressible Lamé lambda or an incompressible longitudinal wave.

The isotropic average is an expectation over exact backgrounds, not a claim
that their superposition remains Beltrami or that independently oriented
periodic cells can be glued without an interface construction. Correlated
mean/internal perturbations and the kinetic reduction still need to be
computed from the joint action before adding coefficients into N3/N4.

## Evidence and route verdict

The first scientific run of verify.py exits zero with 12/12 exact checks
(stdout.txt, empty stderr.txt). It checks all five independent traceless
strain components, derives rather than supplies the shear modulus, verifies
the curl inverse directly for an exact finite affine deformation, detects a
fixed-wavevector mutation, and tests rigid rotation and zero-amplitude limits.
The sphere-moment API is used as conditional unpromoted P242 infrastructure;
its accepted-status limitations are unchanged.

Route verdict: established for the exact isotropic affine Biot–Savart shear
energy of a declared stationary Beltrami-cell ensemble. This advances the
translational energy dependency without closing the full coupled action.

Reusable extraction: euler_affine.py computes the full finite affine energy,
including solenoidality and determinant-one input checks. Four independent
API tests pass (pytest-first.stdout.txt). The extraction replay adds a direct
API/curl-inverse comparison: 13/13 checks pass in api-replay.stdout.txt. The
original 12-check receipt is preserved. No scientific equation changed.
