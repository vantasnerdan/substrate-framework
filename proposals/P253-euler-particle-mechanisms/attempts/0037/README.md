# P253/0037: fixed-Casimir threshold dynamics of the exact Euler solitary wave

Root preregisters the next same-field restoring construction after `0034`.
The object is the independently established `0027/0028` exact whole-space
Euler solitary wave on one fixed compact-vorticity column background. The
target is an all-time axisymmetric propagation and modulation theorem in the
physical finite-excess topology, on a declared mixed-Casimir leaf, with the
translation orbit and the Bessel exterior retained. This is a P2 restoring
mechanism prerequisite. It is not yet full three-dimensional persistence or
an electron/neutrino identification.

## Frozen equations, domain, and success statement

Use the exact `0034` state `q=(zeta,xi)`, the full linearized operator

    eta_t=-{Psi_c,eta}-{K eta,zeta_c}
          +(2/r^4) partial_z(xi_c chi),
    chi_t=-{Psi_c,chi}-{K eta,xi_c},

and the physical finite-excess Hodge energy space. The fixed leaf is defined
by zero first variation of every mixed Casimir, equivalently `delta J=0` in
each regular swirl-label chart together with its two-end and flat-exterior
rows. The translation mode remains. The branch-speed derivative is retained
as a Casimir-changing modulation companion rather than inserted into the
fixed-leaf operator.

The positive result sought is one of these exact alternatives:

1. a bounded or explicitly polynomially bounded linear group modulo
   translation on that fixed leaf, followed by a nonlinear modulation
   bootstrap for a natural open neighborhood; or
2. a genuine growing spectral/semigroup mechanism of the same full operator,
   with its physical accessible initial state and observable gain proved.

A finite-time energy estimate, scalar model, formal Hessian, or discarded
continuum does not meet the target.

## Competing analytic routes

**Route A: constrained threshold resolvent.** Conjugate the column operator by
the exact positive `0030` metric, retain the localized solitary coefficients
as an operator-valued axial-frequency coupling, impose the mixed-Casimir and
translation rows, and construct limiting resolvents at the `k=0` Bessel
threshold. The strict `c>c0` phase/group-speed inequality supplies the
off-threshold estimate; its lost factor `|k|` must be resolved rather than
called a spectral gap.

**Route B: derived long-wave normal form.** Project the actual time-dependent
Euler equations onto the `0027` threshold eigenfunction and derive, with
physical coefficients and the Bessel logarithm retained, the leading
Hamiltonian amplitude equation about the exact homoclinic. Prove the
projection and full radial-complement remainder estimates in the same norm.
The stationary identity `A-A''=beta A^2` does not by itself license KdV,
NLS, or a scalar stability theorem.

**Route C: full Leray/Evans representation.** Construct the whole-axisymmetric
first-order spatial dynamical system or Birman--Schwinger family with the
exterior Dirichlet-to-Neumann map and accessible boundary rows, then classify
point spectrum, threshold resonances, and the translation pole. If a growing
mode appears, lift it to the physical Euler velocity and coadjoint tangent.

A route failure activates repair of its representation and then the other
routes. The nonaxisymmetric azimuthal sectors remain the next required
extension after an axisymmetric theorem.

## Selection criteria and oracle

Selection is fixed by: exact same-field Euler dynamics; preservation of the
mixed-Casimir leaf and physical pressure; a closed operator domain; uniform
control of the zero-frequency exterior; explicit translation modulation;
and a norm that controls the finite-excess kinetic perturbation. The strongest
oracle is an exact resolvent/energy derivation plus limiting-case and mutation
checks in an importable operator API where useful.

No production numerics are licensed here. If analysis leaves a soft
eigenvalue, Evans determinant, or small sign as the only remainder, a separate
small-ratio numerical design will freeze its gauge, continuum truncation,
error budget, convergence ladder, and verdict bridge before values are
inspected.

Root owns `0037`, central records, APIs/tests, and commits. A separate worker
will independently review any proposed same-field theorem. No accepted claim,
release, or parent state changes at preregistration.
