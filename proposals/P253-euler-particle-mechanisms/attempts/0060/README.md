# P253/0060: neutral compact-ring composite and controlled effective band

Status: **activated and author-complete, with the bounded 0061 review
correction recorded below.** The preregistered activation boundary was met by
`0060/activation-schema.exit=0` before any body or calculation was opened.
Root owns the attempt and central integration.

## Parent purpose and fixed boundary

P253/0055 proves that bare incompressible Euler has no strict full-field
causal cone and that a generic compact swirl has an `r^-3` pressure
quadrupole. It also identifies the stronger lattice moment threshold: in
three dimensions a coupling `r^-p` is absolutely summable for `p>3`, while a
finite first spatial moment requires `p>4`. This attempt executes the positive
failure-derived candidate rather than stopping at that obstruction.

Construct an exact same-substrate neutral composite cell from a finite
proper-rotation orbit of one accepted compact-velocity stationary
Euler ring, with disjoint supports and matched exterior pressure constants.
The field sum must itself solve Euler. Use its actual quadratic stress moments
to cancel the pressure quadrupole and the next odd spatial moment, derive the
first surviving far tensor, and test whether the resulting cell supports an
exact or controlled Bloch/large-array band with finite first-moment coupling.

This is a propagation and composite-carrier candidate. It is not called an
electron, neutrino, charge, quantum state, or stable particle. A symmetric
ensemble average is insufficient: cancellations must hold for each physical
cell and for the actual invariant perturbation sector used by the band.

## Fixed cell and exact superposition test

Let `(u_*,p_*)` be one smooth compact-velocity stationary Euler ring at fixed
density, with pressure chosen constant outside its compact support. Choose a
finite point group `G` containing rotations that send its oriented axis to
each signed coordinate axis, and choose cell centers so the transformed
supports are disjoint. Define

    u_cell(x)=sum_(g in G/G_*) g_*u_*(x-X_g),
    p_cell(x)=sum_(g in G/G_*) p_*(x-X_g).                (1)

First prove directly that all cross-advection terms vanish and that (1) is an
exact stationary Euler field. Retain each component's positive energy,
material tags, action, and support. Check whether the selected group includes
orientation reversal through proper rotations of the physical vortex axis;
do not introduce signed algebraic weights that are not realized by actual
Euler copies.

## Route A: cellwise multipole cancellation

For `T_ij=u_i u_j`, define its compact moments

    M_ij=integral T_ij dx,
    M_ij,k=integral x_k T_ij dx,
    M_ij,kl=integral x_k x_l T_ij dx.                    (2)

Use the actual group orbit to prove

    sum_cell M_ij=m delta_ij,
    sum_cell M_ij,k=0,                                   (3)

so contraction with the harmonic Newton Hessian cancels the `r^-3` pressure
term and central symmetry cancels the `r^-4` term. Compute, rather than name,
the irreducible cubic part of `M_ij,kl` and the first nonzero pressure and
acceleration tails. The target is pressure `O(r^-5)` and acceleration
`O(r^-6)`, with an exposing configuration showing whether the coefficient is
nonzero. If cubic symmetry cancels that term too, record and use the stronger
decay.

Then linearize the moments under the proposed internal perturbation. The
band sector earns the cancellation only if its tangent representation leaves
the dangerous lower multipoles zero. If generic symmetry-breaking modes
restore them, retain those modes in the complement and quantify leakage.

## Route B: exact periodic/Bloch linear Euler band

Place identical physical cells on a declared Bravais lattice with separation
larger than the support diameter. The base field remains an exact disjoint
stationary Euler assembly. Construct the actual linearized Euler/Leray
operator, its Bloch fibers, and a finite-dimensional cell mode only after an
invariant Riesz projection or a controlled Feshbach complement has been
proved. A prepared C-CST response is not substituted for autonomous band
dynamics.

For the physical intercell kernel `J(R)`, prove

    sum_R ||J(R)||<infinity,
    sum_R |R| ||J(R)||<infinity,                          (4)

and the analogous derivative/complement estimates. Equation (4), not mere
absolute summability, is the target needed for a bounded group-velocity
symbol. Derive `omega_a(k)` and `grad_k omega_a(k)` from the exact Bloch
generator, preserve its KKS/energy sign, and construct wavepackets with a
finite-time leakage estimate. The full elliptic pressure tail remains
explicit, so success is a strong effective band with algebraic leakage rather
than an exact Lorentz cone.

## Route C: isolated composite and scattering sector

If the periodic Bloch representation is physically too restrictive, use two
widely separated exact composite cells. Derive their leading interaction,
center/internal modulation, radiation and deformation from the same Euler
field. Test whether the symmetry-protected subspace persists under encounter
and whether the first spatial moment gives a polynomial scattering cone. This
route succeeds only with an actual persistence interval or nonlinear
modulation estimate; a static multipole table alone is Route-A progress.

## Competing configurations and selection criteria

Compare before any coefficient evaluation:

1. the six signed-axis orbit with octahedral symmetry;
2. an inversion-paired tetrahedral orbit with fewer components; and
3. a concentric or linked neutral molecule, only if exact disjoint support or
   an actual same-field matching theorem exists.

Select by exact Euler superposition, cancellation without negative weights,
positive energy/action, invariant tangent sector, component count, first
moment summability, robustness to admissible perturbations, and compatibility
with the persistent-carrier work. No empirical particle value is a selector.

## Success contract and continuation

0060 earns its full result only by establishing the exact cell, its first
surviving physical far tensor, the invariant/control status of the canceled
moments, and either the Bloch band plus (4) or the isolated-scattering
alternative. A refuted symmetry configuration activates the next registered
orbit or a nonlinear multipole-compensation mechanism in-run.

Even success supplies only a propagation/composite-carrier ingredient. The
physical positive doublet, action selection, detector, exchange character,
charge/current, electron and neutrino mechanisms remain separate P253
obligations, and the parent campaign remains active.

## 0061 bounded orbit/domain correction

The stationary six-site cell and the 24-element proper-octahedral tensor
identities remain exact, but they are different objects. The six signed-axis
sites form `O/C4`; four rotations stabilize each physical component. Thus a
generic nonzero seed need not have a nonzero `C4` average. The neutral
`r^-5` tangent remains a concrete candidate conditional on exhibiting a
nonzero `C4`-invariant compact DA tangent (then summing six cosets), using a
free 24-site orbit, or directly proving a nonzero overlapping-stabilizer
average. The repaired construction also retains translated vorticity moments,
componentwise circulation/center/harmonic rows, and the full pressure/Leray
graph-domain estimate still needed before the scalar `r^-5` lattice bound can
be promoted to a `C^1` Bloch operator.
