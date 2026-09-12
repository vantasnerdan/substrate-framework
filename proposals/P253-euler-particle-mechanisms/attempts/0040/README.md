# P253/0040: nonaxisymmetric rotating-wave branch from a stable Cao ring

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0040/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0040 artifact. No derivation body, verifier,
numerical run, source edit, accepted API, claim, or receipt is produced. Root
owns central files, commits, activation, review, and promotion. This worker
owns only eventual append-only artifacts in `attempts/0040`.

## Parent objective and bounded positive intent

The P253 objective remains an actual robust localized Euler mechanism for the
electron first and the neutrino subsequently. Attempt 0040 addresses one
LP2/P4 candidate: construct one exact, compact-vorticity, finite-kinetic-energy
Euler rotating wave which bifurcates from an orbitally stable classical Cao
vortex ring and carries a nonzero classical rotational action. Success would
provide a time-dependent, nonaxisymmetric, same-carrier collective sector. It
would not identify that sector with a particle or supply quantum mechanics.

Fix `kappa>0`, `W>0`, an integer `p>=2`, and one sufficiently small
cross-section parameter `epsilon=epsilon_0` in Cao--Lai--Qin--Zhan--Zou,
arXiv:2206.10165. Its no-swirl potential vorticity is

    zeta_epsilon=omega_epsilon^theta/r
      =epsilon^-2(psi_epsilon-(W/2)r^2 log(1/epsilon)
                    -mu_epsilon)_+^p,                       (1)

with circulation `kappa`, translating speed
`c_epsilon=W log(1/epsilon)`, compact vorticity, and limiting ring radius
`r_*=kappa/(4 pi W)`. The cited theorem gives `C^3` velocity and `C^2`
vorticity for `p>=2`; this attempt does not silently upgrade that statement to
`C-infinity`. If a bifurcation theorem needs more derivatives, the regularity
earned from the chosen larger integer `p` and the free-boundary composition
`s_+^p` must be proved explicitly.

| Frozen field | 0040 value |
| --- | --- |
| `positive_intent` | One exact nonaxisymmetric rotating-wave branch on the coadjoint leaf of a fixed stable Cao ring |
| `primary_object` | A compact-vorticity finite-energy Euler solution with axial translation and physical `SO(2)` rotation, bifurcating from an isolated azimuthal Kelvin/shape mode |
| `carrier_scope` | The polynomial-profile no-swirl family in arXiv:2206.10165, not the unrelated swirling family in arXiv:2009.13210 |
| `required_structure` | Exact full 3D Euler operator, dynamically accessible leaf, translation/rotation/isotropy reduction, isolated mode or an executed embedded-mode replacement, physical period, and nonzero KKS/angular-impulse pullback |
| `pass_licenses` | A classical same-carrier relative equilibrium/relative-periodic mechanism and its classical rotational action |
| `does_not_license` | Nonlinear stability of the new branch; half-integer spin; FR sign; `hbar`; Born rule; quantization; relativity; electron or neutrino identity |
| `maximum_verdict` | Existence and local control of the named exact branch, with its physical `SO(2)` period and nonzero classical action coefficient |
| `failure_scope` | Only the selected carrier member, azimuthal sector, and route whose hypotheses were actually tested |

## Primary-access and operator audit before source work

Full source bodies are cached outside the campaign tree, as required:

| Source | Access receipt and exact locations | Imported fact and boundary |
| --- | --- | --- |
| Cao, Lai, Qin, Zhan, Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation* | [arXiv:2206.10165v2](https://arxiv.org/abs/2206.10165v2), cached `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`; v2 dated 2023-12-25; no-swirl reduction (1.3)--(1.6), PDF pp. 1--2; Proposition 1.4 and Theorem 1.6, PDF pp. 5--6; Theorem 1.8, PDF pp. 6--7; Appendix C, PDF pp. 59--62 | Exact polynomial-profile compact-core translating rings, uniqueness up to axial translation, classical regularity, and nonlinear orbital stability **only for nonnegative axisymmetric potential-vorticity perturbations** in the stated weighted `L1 intersect L2` setting. It supplies neither a full 3D nonaxisymmetric stability theorem nor a Kelvin eigenmode. |
| Cao, Qin, Yu, Zhan, Zou, *Existence, uniqueness and stability of steady vortex rings of small cross-section* | [arXiv:2201.08232v4](https://arxiv.org/abs/2201.08232v4), cached `/tmp/primary-source-cache/P253-0040/2201.08232.pdf`, SHA-256 `1c432136e5546299c7d47e86d21d14e231369347fb086da58b03607127a19156`; v4 dated 2023-12-05; Theorems 1.1, 1.2, 1.5, PDF pp. 6--9; elliptic linearization (2.18) and limiting kernel discussion (2.22), PDF pp. 16--17 | Patch-family comparator and the precise origin of the Lyapunov--Schmidt translation kernel. Its surface-delta elliptic operator is not a time-linearized Euler generator and its lower regularity is not the 0040 target. |
| Pocklington, *The complete system of the periods of a hollow vortex ring* | [DOI 10.1098/rsta.1895.0017](https://doi.org/10.1098/rsta.1895.0017), *Phil. Trans. R. Soc. A* **186** (1895), 603--619; Royal Society metadata and a third-party full-text listing are accessible, but no theorem body is imported before activation and a versioned PDF/hash is not yet pinned | Historical thin hollow-core mode calculation only. It is a Route-B model candidate, not evidence that any mode persists for the smooth finite-core Cao operator. |
| Cao--Zhan, *On the steady axisymmetric vortex rings for 3-D incompressible Euler flows* | [arXiv:2009.13210v2](https://arxiv.org/abs/2009.13210v2), existing 0010 cache receipt SHA-256 `d5da2b24ed5822e52a1d7a30f0dc3dd23289db15f91d1c2e133f12807ab02aed`; Theorem 1.3, PDF pp. 6--7 | An exact axisymmetric **swirling** construction whose paper leaves stability open. It is not the orbitally stable no-swirl Cao carrier frozen here. |

The “linearized operator” used in the Cao uniqueness proofs is an elliptic
stream-function operator for a stationary Lyapunov--Schmidt/Pohozaev argument
(and, in the patch paper, contains a free-boundary surface delta). It is not
the full time-dependent three-dimensional Euler linearization. No eigenvalue
or invertibility statement about that elliptic operator is transferred to the
dynamical problem.

## Exact full-ring dynamical operator

Work in the frame translating with the fixed ring. Write

    U_epsilon=(u_epsilon^r(r,z),0,
               u_epsilon^z(r,z)-c_epsilon),
    Omega_epsilon=r zeta_epsilon(r,z) e_theta.               (2)

On the finite-energy solenoidal space `X=L2_sigma(R3)` the velocity
linearization is the closed realization of

    G_epsilon v
      =-P_L[(U_epsilon dot grad)v+(v dot grad)U_epsilon],     (3)

where `P_L` is the whole-space Leray projector. Equivalently, with
`eta=curl v`,

    partial_t eta =-(U_epsilon dot grad)eta
                    -(v dot grad)Omega_epsilon
                    +(Omega_epsilon dot grad)v
                    +(eta dot grad)U_epsilon.                (4)

The domain is the maximal transport/shear domain on which (3) lies in `X`,
intersected with the dynamically accessible closure; an unproved global
`H1` domain is not substituted.

Because the base is axisymmetric, (3) reduces under
`v(r,theta,z)=v_l(r,z)e^(i l theta)` for each integer `l`. In cylindrical
components the divergence constraint is

    r^-1 partial_r(r v_l^r)+i l r^-1 v_l^theta
      +partial_z v_l^z=0,                                   (5)

and all cylindrical connection terms are retained. In particular the theta
component of `(v dot grad)U_epsilon` contains
`U_epsilon^r v_l^theta/r`; the pressure has
`i l pi_l/r` in that component and solves the full Fourier-sector Poisson
equation with

    Delta_l=partial_rr+r^-1 partial_r+partial_zz-l^2/r^2.    (6)

Thus neither a scalar meridional operator nor a filament dispersion relation
is the object whose spectrum is needed.

Let `O_epsilon` be the coadjoint orbit of the fixed Euler momentum. Its
integrable compactly generated tangent is represented at vorticity level by

    delta Omega=curl(xi cross Omega_epsilon),
    div xi=0,  xi in C_c^infinity,                           (7)

with the corresponding Biot--Savart/Leray velocity. The analysis uses the
closure of (7) in the declared energy-plus-graph topology. It removes axial
translation, fixes circulation and hydrodynamic impulse, imposes a center
slice, and quotients the actual stabilizer. It does not claim a general
Clebsch chart or identify Cao's axisymmetric scalar rearrangement class with
the entire three-dimensional coadjoint orbit.

## Frozen eigenmode and exact branch target

Choose one azimuthal number `l>=2`. The primary spectral achievement is a
simple isolated Hamiltonian pair

    G_epsilon q_l=i nu_l q_l,
    G_epsilon conjugate(q_l)=-i nu_l conjugate(q_l),
    nu_l not=0,                                               (8)

on the reduced dynamically accessible sector, with finite nonzero KKS/Krein
signature and the spectral separation/nonresonance required by the selected
infinite-dimensional bifurcation theorem. The `l=0` family tangents and the
`l=1` translation/tilt sectors are not counted as the desired shape mode.

At the critical pattern speed

    Omega_0=-nu_l/l                                            (9)

with sign fixed by the `e^(i l theta)` and rotation conventions, construct a
real branch `(V_a,c_a,Omega_a)`, `0<|a|<a_0`, solving the exact stationary
Euler equation in the translating and rotating frame. In invariant notation,

    F(V,c,Omega)
      =P_L[(V dot grad)V]-c partial_z V-Omega [R,V]=0,
    R(x)=e_z cross x,                                         (10)

where the bracket sign is checked against the physical ansatz rather than
inferred from notation. The resulting laboratory solution is

    u_a(t,x)=R_{Omega_a t}
       V_a(R_{-Omega_a t}(x-c_a t e_z)).                     (11)

It has `C_l` isotropy and no larger accidental stabilizer. Consequently the
full `SO(2)` group circuit is

    T_SO2(a)=2 pi/|Omega_a|,                                  (12)

while the minimal repetition time of the physical field after quotienting
its `C_l` isotropy is

    T_phys(a)=2 pi/(l |Omega_a|).                             (13)

Both periods, the axial displacement `c_a T`, and the exact relative-periodic
identity are verified. A frequency in a modal equation without (11)--(13) is
not a physical rotating wave.

The branch remains on the fixed carrier leaf:

    Omega_a=(g_a)_* Omega_epsilon,
    g_a volume-preserving,
    ||g_a-id||_{C^k}+||V_a-u_epsilon||_{H^s}=O(|a|),          (14)

in explicitly chosen indices supported by the source regularity and the
bifurcation estimates. Hence its vorticity support is a compact deformation
of the same torus and its Biot--Savart velocity has finite kinetic energy.
Equation (14), uniform inverse bounds on the reduced complement, and an
explicit `O(a^2)` nonlinear remainder are the required same-family
perturbation control. Cao's axisymmetric orbital-stability theorem is retained
as a base-family fact; it is not promoted to stability of (11).

## Rotational moment map and KKS pullback

Set the physical mass density to the declared constant `rho_0`. On compact
vorticity, define the axial angular impulse

    J_z(Omega)=-rho_0/2 integral_R3 |x|^2 Omega_z(x) dx.      (15)

Its variational derivative has curl `e_z cross x`, so (15) is the
vorticity-side moment map for physical axial rotations. The equivalence with
a cutoff limit of `rho_0 integral (x cross u)_z dx`, including the far-field
boundary term and origin dependence under nonzero linear impulse, must be
proved in the chosen centered/impulse-fixed slice; an improperly convergent
velocity integral is not used as a definition. The centered no-swirl base has
`J_z=0` because `Omega_epsilon` is purely azimuthal.

Let `theta` be the actual physical rotation angle and let
`iota(a,theta)=(R_theta)_*Omega_a`. Fix the convention

    i_{X_J} Omega_KKS=dJ_z.                                   (16)

The exact moment-map identity then requires

    iota^* Omega_KKS=d theta wedge d J_z(a),
    J_z(a)=j_2 a^2+O(a^3),  j_2 not=0.                       (17)

The vanishing of the linear term follows from angular Fourier orthogonality;
the nonzero coefficient must be calculated from the normalized eigenmode,
its second-order correction, and the actual density/sign conventions. A
one-dimensional orientation circle has zero pulled-back two-form; the
two-dimensional amplitude--orientation surface in (17) is essential. The
classical action around a full group circuit is
`2 pi J_z(a)` (and around the physical `C_l` quotient circuit is
`2 pi J_z(a)/l`, when the quotient descends). No integer, `N=1`, or `hbar`
fitting is permitted, and geometric prequantization is not quantum dynamics.

## Route A: exact reduced Crandall--Rabinowitz/Lyapunov center

Build a local chart on the integrable part of `O_epsilon`, impose the axial
center/impulse slice, and quotient the base `SO(2)` stabilizer. Treat `Omega`
in (10) as the bifurcation parameter: every axisymmetric base solves (10),
and at (9) the `l`-sector acquires the kernel generated by `q_l`. Prove on the
chosen Banach/Hilbert scales:

1. Fredholm index zero and a one-dimensional complex kernel after the real,
   symmetry, and Casimir reductions;
2. a closed range and uniformly bounded inverse on its complement;
3. the CR transversality of `partial_Omega D_V F q_l`, or the corresponding
   Hamiltonian Lyapunov-center nonresonance and finite KKS signature;
4. tame differentiability at the compact-vorticity boundary, with no
   derivative count beyond that earned from (1);
5. exact reconstruction of (11), (14), and the coefficient in (17).

If a direct Eulerian chart loses the coadjoint constraint, switch within this
route to a volume-preserving Lagrangian displacement chart. Failure of a
generic textbook CR theorem on an unbounded-domain transport operator is a
method failure, not evidence that the branch cannot exist.

## Route B: thin-ring Kelvin/Pocklington asymptotic plus finite-core inverse

Let `delta=epsilon/r_*` and derive the rescaled exact sector operator
`A_{delta,l}` from (3)--(6) on the Cao core. Separately derive the hollow,
filament, or constant-core comparison operator `A_l^model` from an opened
primary source, including its density profile, boundary conditions, frame,
time scale, and mode normalization. The bridge requires an estimate of the
form

    A_{delta,l}=s_delta(A_l^model+E_{delta,l}),
    ||E_{delta,l}||_{D -> X}<=r_delta ->0,                   (18)

together with a contour `Gamma_l` enclosing one nonzero imaginary model mode,

    sup_{z in Gamma_l}||(A_l^model-z)^-1||<=C_l,
    C_l r_delta<1,                                           (19)

and a uniform inverse for the exact reduced complement. The physical scale
`s_delta`, core corrections, logarithms, and spectral contour are derived,
not fitted. Only (18)--(19), or an equally strong Evans/Fredholm comparison,
can transfer a Pocklington/Kelvin mode to the smooth finite-core Cao ring.

This route must control the exterior Biot--Savart/Leray field, the core edge,
translation and rotation zero modes, and any essential spectrum. A filament
frequency, formal local induction law, or numerical dispersion curve alone
does not meet the target.

## Route C: embedded mode, reversible/KAM, or contour-patch replacement

If (8) lies in continuous/essential spectrum or the reduced inverse in Route A
is not Fredholm, first compute the exact resonance mechanism. Competing
continuations are:

- a reversible Hamiltonian Lyapunov/KAM or Nash--Moser construction using a
  derived twist, Diophantine exclusions, smoothing estimates, and control of
  the Euler derivative loss on the same smooth Cao leaf;
- a contour-dynamics bifurcation from the Cao patch family of
  arXiv:2201.08232, where the free-boundary operator is genuinely the correct
  evolution object, followed by a separately proved continuation/limit to the
  polynomial-profile carrier if the smooth-family target is retained.

A patch branch is a valuable route-scoped positive object but is not silently
renamed a branch of the classical polynomial family. If the embedded mode
dissolves into phase mixing, the attempt records that mechanism and generates
the next regular-carrier or vortex-tube candidate; it does not infer a global
Euler no-go.

## Analytic boundary, route competition, and completion contract

0040 begins with exact source/operator/spectral calculus. No production
numerics are preregistered or authorized. Symbolic checks after activation may
verify cylindrical signs, moment-map identities, asymptotic recurrences, and
normal-form coefficients, but sampled spectra cannot establish isolation or
the uniform inverse. If an irreducible small spectral splitting or stability
edge survives the analytic ladder, a new append-only numerical design must
first load the small-ratio skill, pin error scales and convergence bounds, and
receive its own authorization.

Routes are compared by exact retention of the Cao carrier and coadjoint leaf,
strength of the spectral/inverse theorem, physical period reconstruction,
nonzero KKS coefficient, finite-energy pressure/tail control, and assumption
cost. Each attempted route receives exactly one route-scoped verdict:
established, refuted with the mechanism named, or blocked by a precise missing
construction. Method repair, representation change, and the next competing
route remain active in-run.

The positive contract is met only by an exact solution (11) satisfying
(12)--(17) on the same compact-vorticity finite-energy carrier leaf, with the
required nonlinear construction estimates. An isolated formal eigenmode,
the Cao base stability theorem, a filament limit, an abstract Hamiltonian
bifurcation API, or a nonzero classical action by itself does not complete the
attempt. After success, the next LP2/P4 dependency is nonlinear orbital or
modulational control of the new branch; quantum and relativistic bridges
remain separate constructions.
