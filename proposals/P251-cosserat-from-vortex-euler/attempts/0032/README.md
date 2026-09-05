# Attempt 0032: closed relative-angle action in a six-vortex molecule

Parent: P251 / issue #198 / PR #199. Active obligations: N3/N4, with N2 and N7
independent. Base release v0.171.0. Owner: delegated `orientation_construction`;
write surface: this directory. No accepted claims or public APIs change.

## Candidate generation and exact scope

Attempt 0030 derived a positive patch angle oscillator but its exterior flow
was prescribed. The next achievement is a closed microscopic action including
reaction. Two finite patches have rigorous relative-equilibrium existence
theorems, but their arbitrary ellipse shapes are not dynamically invariant.
This motivates a different candidate: two concentric regular triangles of
equal point vortices, forming a regular hexagon at equilibrium. Discrete C3
symmetry makes the two-triangle manifold exactly invariant under point-vortex
dynamics. This provides a completely closed, non-postulated angle action.

This is the singular two-dimensional Euler point-vortex model, equivalently
parallel ideal line vortices per unit length. Finite core, unrestricted
three-dimensional Euler dynamics, and smooth stationary EPS tubes remain
distinct constructions. No gap or comparator was used to choose coefficients.
All calculations below are exact algebra, with a separate Cartesian
Biot–Savart Jacobian oracle; no numerical solver or small-ratio estimate occurs.

## Exact microscopic energy and reduction

Let each of the six vortices have circulation Γ>0. Their locations are

    z_(a,j)=r_a exp(iθ_a+2πij/3), a=1,2; j=0,1,2.

At distinct positions the renormalized physical energy per axial length is

    H=−ρΓ²/(2π) sum_(i<j) log|z_i−z_j|.

The omitted self/outer constants are position-independent at fixed circulation.
The Euler point-vortex equations are

    z_i_dot=Γ/(2π) sum_(j≠i) J(z_i−z_j)/|z_i−z_j|²,
    J(x,y)=(−y,x).

The action is `sum p_a θ_a_dot−H`, with p_a=−3ρΓr_a²/2. Pair interactions are
counted once, and differentiating their dependence on θ1−θ2 gives equal and
opposite torques. No exterior field is imposed.

The roots-of-unity product `prod_(j=0)^2(z−w exp(2πij/3))=z³−w³` exactly
reduces all nine cross interactions. Set

    S=r1²+r2², x=(r1²−r2²)/S, χ=θ1−θ2,
    A=3ρΓ²/(4π), P0=3ρΓS/4.

Common-rotation symmetry conserves total angular impulse. At fixed S the
canonical relative momentum is p_χ=−P0 x, up to an irrelevant constant in
another choice of common angle. Apart from S-dependent constants,

    H_red/A=−log(1−x²)
      −log[2+6x²−2(1−x²)^(3/2) cos(3χ)].

The physical domain is |x|<1 with all inter-vortex distances nonzero.
The staggered state x=0, χ=π/3 is the regular hexagon. Its Hessian is

    H_xx=A/2, H_χχ=9A/2, H_xχ=0.

Eliminating the canonical momentum gives the exact linear oscillator action

    L2=I χ_dot²/2−K(χ−π/3)²/2,
    I=2P0²/A=3πρS²/2,
    K=9A/2=27ρΓ²/(8π),
    frequency²=K/I=9Γ²/(4π²S²).

I has units kg m and K has units N per axial length. Both are positive. This
is a stiffness of ANGLE, with conjugate radial redistribution providing the
inertia; it is not a norm of an independently prescribed angular rate.

The script independently differentiates the full Cartesian Biot–Savart system
at the regular hexagon of radius r=sqrt(S/2). After subtracting its common
angular velocity Ω=5Γ/(4πr²), alternating radial/tangential perturbations obey

    radial_dot=9Γ tangential/(4πr²),
    tangential_dot=−Γ radial/(4πr²).

Their frequency matches the action-derived value. The check directly detects
wrong signs or factors in K. The common rotation is a relative equilibrium,
not a stationary six-vortex field in the laboratory frame.

## Restoring the cage: the actual kinetic map

The common mean angle is B=(θ1+θ2)/2. Before fixing its impulse, the complete
Hamiltonian includes `−5A log S`. Let P=p1+p2=−3ρΓS/2; expansion about the
rotating hexagon gives common inertia

    I_B=P²/(5A)=2I/5.

Thus the complete two-angle quadratic kinetic term in the rotating frame is

    T2=I_B B_dot²/2+I χ_dot²/2.

If θ1 is the physical microangle Φ and θ2 is the cage angle β, this is

    T2=I_B(Φ_dot+β_dot)²/8+I(Φ_dot−β_dot)²/2.

Its cross inertia is `∂²T2/(∂Φ_dot∂β_dot)=−9I/10`, which is nonzero. At fixed
total impulse only the relative-rate term remains. Therefore replacing this
action by `j Φ_dot²/2` while assigning `β=curl(u)/2` changes the equations.
The exact action gives a generalized micropolar model with spin/ambient-rate
coupling and gradient inertia; the standard target requires another justified
field/kinetic map or a stated asymptotic restriction. At optical frequencies
the omitted cross inertia is of the same order as the angle spring, so a
long-wavelength argument alone does not remove it uniformly.

This is a concrete condition on the proposed continuum bridge, not an absence
of reaction energy: reaction is already present in the derived action. A
constant-angle change of variables can diagonalize any positive mass matrix,
but its new field must be identified with a physical collective orientation,
and its transformed potential and macro field must be carried along. Renaming
the diagonal coordinate as the original core angle would hide the change.

## Material-frame and spatial-gradient constructions

For an affinely displaced triangular cage with unit directions n_j, the average
small angle is `(1/3) sum (J n_j)·(grad u)n_j = curl(u)/2` in its plane.
The script proves this for a fully general 2×2 displacement gradient. Hence
the elastic-angle map is a geometric consequence of an explicit affine-cage
premise, rather than identification of a velocity-strain eigenframe with a
material angle. Symmetric strain also deforms the cage; the full Hessian must
retain its other irreducible modes for constitutive closure.

For isotropically oriented tube axes n, the conditional relative-angle map
`δχ=n·(Φ−curl(u)/2)` would give `alpha=L_m K/12`, where L_m is the axial length
density of MOLECULES. The six-filament length density is `L_v=6 L_m`; using the
same symbol for both would introduce a factor of six. The formal relative
inertia average is `L_m I/3`, but multiplies the relative rate as derived above.

There is also a concrete log-leading spatial construction: replace each point
by a thin, circulation-pinned filament and vary its azimuth along axial s.
Local-induction self energy of a helix is

    T sum_(a,j) sqrt(1+r_a² θ_a_prime²),
    T=ρΓ² log(R/a_core)/(4π).

At the regular hexagon its quadratic part is

    W_line=3TS B_prime²/2+3TS χ_prime²/8.

Consequently the relative twist modulus is `C_χ=3TS/4`, while common twist has
`C_B=3TS`. If only the inner triangle twists and the cage is fixed, the
coefficient is `3TS/2`; these are different admissible variations. An expression
using the same χ for all six filament angles describes common rotation and
cannot be assigned the relative modulus. The script checks these distinctions.

This spatial construction is leading logarithmic, not the exact finite-core
Biot–Savart interaction. Core radius, outer scale, slenderness, slow axial
variation, and discarded mutual finite terms must be declared. Isotropic
fourth moments can contract its axial derivative sector into wryness, but a
complete three-dimensional triad and all relevant bending modes remain needed.

## Finite-core route and source assessment

Claudia García's Theorem 1.1 constructs genuinely rotating Euler patches near
each vertex of any regular N-gon, including N=6, by solving contour dynamics.
Thus finite-core relative-equilibrium existence is available from a theorem.
It does not by itself identify the above relative-angle spectral subspace,
prove convergence of its constrained energy and symplectic form, or supply an
exact finite-core value of I and K. Those are the next concrete obligations.

The Hamiltonian stability/desingularization work of Wan and the invariant
subspace estimates of Long–Wang–Zeng give routes to that spectral continuation;
their exact hypotheses must be checked for the rotating unbounded-plane family
before their conclusions are imported. Stability of a point configuration
alone is not a substitute for that verification. EPS smooth stationary tubes
are a separate three-dimensional compatibility question, not established by a
two-dimensional rotating-patch theorem.

## Route verdict and continuation

`route_verdict: established` for the exact closed six-point-vortex relative
angle action and its positive canonical inertia/stiffness. `evidence_scope:
EXACT_POINT_VORTEX_INVARIANT_SECTOR`.

The conditional log-leading helical construction is analytic asymptotic
progress. The full campaign stays active. Next executable constructions are
finite-core spectral continuation, the explicit generalized-to-standard kinetic
field map, and the affine three-dimensional ensemble/EPS embedding. None of
these local results promotes the original C-CST claims or closes N3/N4/N7.

Reproduce with `PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python
proposals/P251-cosserat-from-vortex-euler/attempts/0032/hexagon_action.py`.
The first script execution is captured in `stdout.txt` and `stderr.txt`:
24 checks passed, process exit zero, and empty stderr.

## Continuation: covariant collective field and controlled dispersion order

The cross-inertia finding above activates an exact coordinate repair. Write
the kinetic matrix in (Φ,β) as a=c=11I/10, b=−9I/10. The new LOCAL collective
angle and inverse are

    Ψ=(aΦ+bβ)/(a+b)=(11Φ−9β)/2,
    Φ=(2Ψ+9β)/11.

Under a common infinitesimal physical rotation both Ψ and β shift by that
rotation. This is therefore a valid covariant local collective-angle chart,
although Ψ is not the original inner-triangle angle. The half-integer weights
do not define a global map on unlabeled periodic triangle orientations; the
claim here is the linear theory near the selected equilibrium.

The exact transformed action has

    T2=J_Ψ Ψ_dot²/2+J_β β_dot²/2,
    J_Ψ=2I/55, J_β=4I/11,
    V2=K_Ψ(Ψ−β)²/2, K_Ψ=4K/121.

Thus a standard positive microinertia is obtainable by a declared physical
field map, rather than removing a cross term by hand. The retained cage inertia
becomes macro gradient inertia when β=curl(u)/2. The transformed helical
self-energy has `C_Ψ=12Tr²/121=6TS/121` multiplying Ψ_s²/2, as well as explicit
Ψ_s β_s and β_s² terms. These gradient cross terms remain in the exact action.

Conditional on a bulk kinetic term ρ|u_dot|²/2, the cage inertia adds an O(k²)
mass correction. It changes the transverse determinant by a multiple of

    k² z (j z−K0−C0 k²), z=frequency².

At the acoustic root z=O(k²), and at the optical root z=K0/j+O(k²), that
difference is O(k⁴). Both roots are simple as roots in z at k=0, so the
dispersion corrections are O(k⁴). The transformed helical cross term changes
the coupling from O(k) to O(k)+O(k³), also first altering the determinant at
O(k⁴). This supplies a precise long-wave sense in which the leading standard
micropolar dispersion can agree with the exact generalized action on both
branches, instead of an unjustified dismissal of optical-frequency inertia.

One ensemble distinction survives the coordinate repair. The free isolated
hexagon has relative frequency K/I. If the cage is constrained by affine
macroscopic displacement, its spatially uniform optical motion has β=0 and

    K_Ψ/J_Ψ=10K/(11I).

The different frequency follows from different admissible cage motion, not
from an algebraic error. A physical intercell mechanism enforcing affine cages
and carrying their reaction is still needed for that continuum claim. The
positive bulk mass term is also not supplied by the isolated molecule's center
coordinate: a molecule of circulation 6Γ has a first-order translational
symplectic term. A balanced ±Γ ensemble can cancel mean gyro terms, but that
cancellation alone does not derive an invertible positive translational mass.
The bulk Euler mean-flow action and its fluctuation/collective-coordinate
closure must be established without double counting the swirl energy.

Accordingly this continuation establishes an exact local field map and a
conditional asymptotic bridge. It does not replace the user's exact conditional
Euler/EPS objective by a finite-order filament theorem. The now sharper next
construction is a reaction-bearing affine cage/mean-flow closure plus finite
core and EPS compatibility; the local kinetic-map obstruction itself is solved.

`collective_field_map.py` reproduces this continuation. Its first execution is
captured separately in `map-stdout.txt` and `map-stderr.txt`.
It passed 18 checks with process exit zero and empty stderr. Ruff and
`git diff --check` pass for the completed local write surface. The parent's
attempt 0034 independently develops the same collective field map and should
be treated as a related derivation, not another physical source.

## One next executable analytic route

Apply the contour-corrected six-patch family of García to the exact
isovortical Hamiltonian, retaining cage motion. At a rotating state use the
augmented energy `F=H+ΩJ`, where our canonical rotation momentum is `−J`.
For an area-preserving generator χ, derive the full second variation along
`δω={χ,ω}`, including `δ²ω={χ,{χ,ω}}`. Separate the alternating centroid
radial/tangential sector from patch-shape modes and form the constrained
Hamiltonian/symplectic Schur complement. This carries both finite-core
deformation and reaction and has a known exact point-vortex limit provided by
the Cartesian oracle above. Radial-core external logarithmic moments and the
large core-mode frequency scale Γ/a_core² offer analytic estimates for that
limit; pointwise replacement of a finite-core Hessian by the point result does
not.

This route can establish a finite-core Euler optical action and compute its
coefficient corrections without an additional constitutive or cage-locking
assumption. It still does not produce an affine macro cage merely by existing:
the next collective assembly must derive which cage modes are actual
long-wave material modes. The intermediate result is smooth patch-boundary
Euler dynamics, not yet smooth vorticity or stationary three-dimensional EPS
embedding. No terminal or exhaustion conclusion follows from its outcome.

## Primary references (accessed 2026-09-04)

- C. García, *Vortex Patches Choreography for Active Scalar Equations*, J.
  Nonlinear Sci. 31, 75 (2021), [Theorem 1.1 and §2](https://arxiv.org/pdf/2010.07361),
  DOI [10.1007/s00332-021-09729-x](https://doi.org/10.1007/s00332-021-09729-x).
  The displayed Euler kernel fixes normalization; the theorem supplies rotating
  finite-patch polygon existence, not the action coefficients derived here.
- T. Hmidi and J. Mateu, *Corotating and counter-rotating vortex pairs for Euler
  equations*, [DOI 10.5802/jedp.647](https://proceedings.centre-mersenne.org/articles/10.5802/jedp.647/)
  (2016); longer construction [arXiv:1601.02242](https://arxiv.org/pdf/1601.02242).
- L. G. Kurakin and V. I. Yudovich, *The stability of stationary rotation of a
  regular vortex polygon*, Chaos 12, 574–595 (2002),
  [DOI 10.1063/1.1482175](https://doi.org/10.1063/1.1482175).
  Independent context for full point-polygon orbital stability; the script's
  narrower optical-sector check does not rely on this result.
- Y.-H. Wan, *Desingularizations of systems of point vortices*, Physica D 32,
  277–295 (1988), [DOI 10.1016/0167-2789(88)90056-5](https://doi.org/10.1016/0167-2789(88)90056-5).
  Abstract and bibliographic record inspected; its complete hypotheses have not
  been audited here, so stability preservation is a route pointer only.
- Y. Long, Y. Wang, C. Zeng, *Concentrated steady vorticities of the Euler
  equation on 2-d domains and their linear stability*, JDE 266, 6661–6701
  (2019), [DOI 10.1016/j.jde.2018.11.011](https://doi.org/10.1016/j.jde.2018.11.011).
  The invariant slow/fast spectral decomposition is relevant methodology;
  a rotating-plane application remains to be checked.
- D. G. Crowdy, *Exact solutions for rotating vortex arrays with finite-area
  cores*, JFM 469, 209–235 (2002),
  [author manuscript](https://www.ma.ic.ac.uk/~dgcrowdy/PubFiles/Paper12.pdf).
  Exact patch-plus-satellite arrays are a materially different candidate if the
  six-patch continuation does not supply the required continuum map.
