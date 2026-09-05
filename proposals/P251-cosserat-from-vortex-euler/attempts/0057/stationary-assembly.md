# A normalized one-fluid Euler action and its compact EPS reconstruction

## 1. Why the representation changes

The decaying EPS construction supplies a valid whole-space relative orbit
and finite compact defect Hessians. It does not, without another measure,
turn one defect into a positive-density bulk fluid. Dividing its full
exterior-inclusive defect energy by a chosen tube volume and independently
adding all ambient mass is not an assembly proof. The alternative below
constructs an actual stationary random Euler field, rather than gluing
different decaying solutions along unverified faces.

Let u0 have the Gaussian Beltrami law in the source audit, rescaled to the
declared nonzero lambda and physical velocity scale. Its sample fields are
smooth, solve curl u0=lambda u0, and are exact stationary Euler fields with
`p0=-rho |u0|²/2+constant`. Its law is translation invariant, has mean zero
and finite pointwise kinetic and derivative moments. In the source units
`E[u0(x) tensor u0(x)]=I`, so `E[rho |u0|²/2]=3 rho/2`.

Randomly rotate the ENTIRE realization and its geometry by an independent
Haar rotation. Include the physically reflected and time-reversed partners
as declared in the parent ensemble. Every member is still an actual
stationary smooth Euler field. The mixture is stationary and isotropic;
we use its normalized ensemble expectation, not an unproved assertion that
each realization has isotropic spatial averages. Reflection changes the
sign of lambda, and formulas transform with that sign.

The probability measure samples fluid at a typical spatial point. Thus
`E[1]=1`: the material density is rho throughout, including all fluid
outside the selected tubes. A selected tube intensity is not a volume
weight for the entire mass.

## 2. One material action, mean velocity, and no double counting

Start with the single material Euler action, before any core reduction,

    S=integral dt E_material [rho |gdot|²/2
                                +p(det D_a g-1)].

Here E_material means normalized material-volume/ensemble expectation.
For stationary incompressible realizations, a change from material to
Eulerian sampling preserves volume. The same normalization counts all
ambient particles. A partition into reference parcels is only bookkeeping:
their material maps are restrictions of ONE global smooth map, so shared
faces agree, and their pressure-traction terms cancel on internal faces.
Any exterior observation-boundary traction remains explicitly in the
balance. No periodic tiling of independent EPS solutions is assumed.

Define V to be the ENSEMBLE mean material velocity, and v' its fluctuation.
Then the exact identity is

    E_material[gdot]=V,  E_material[v']=0,
    rho E_material |gdot|²/2
                    =rho |V|²/2+rho E_material |v'|²/2.

Individual tube or ambient parcel means can fluctuate; their squared
means remain in the second term. They are not subtracted a second time.
At the quadratic level V=Udot for the mean displacement U. The mean
cotangent term is `rho V.Udot`; eliminating V gives `rho |Udot|²/2`.
This is an energy DENSITY of a normalized infinite-fluid ensemble. It is
not a declaration that a uniform boost of all R3 has finite total energy.

For comparison, using each tube's mass centroid as the macro coordinate
would require the additional parcel-mean Gram and connection of 0052/0055.
We retain those as physical observation maps, but do not silently identify
them with U. This representation preserves mean material momentum,
`P_mean=rho V`, and keeps its full fluctuation energy.

At nonzero macro gradients, any induced polar mean must be removed from
v' by its definition and included in V. For a trial tangent w the exact
quadratic mean subtraction is `rho |E_material w|²`. Uniform axial spin
rates cannot induce a polar ensemble mean in the reflection-paired
isotropic ensemble. At gradient order, curl of an axial response can be
polar, and that cross/observable map remains. Nothing here discards the
parent's gradient inertia or boundary-current improvement.

## 3. Positive-density good EPS patches, with a proved response rank

Choose a reference EPS field on a compact ball containing the desired
robust finite-core tube, physical jets, all cutoffs and response balls.
All cutoffs are geometrical input. To guarantee the Euclidean response
rank without borrowing the old far-field theorem, perturb this field by
an arbitrarily small same-lambda ABC field having nonzero helical Fourier
atoms at each of `+/-lambda e1,+/-lambda e2,+/-lambda e3`.
The perturbation can be small in every required C^m norm on this compact
ball and hence preserves the robust EPS structures and strict local bounds.

This seed has no infinitesimal Euclidean symmetry. Indeed, its decaying
EPS part has a smooth density on the Fourier sphere. A nonzero rotational
generator differentiates at least one of the added atoms tangentially;
neither a smooth density nor an undifferentiated atom can cancel that
distributional derivative. Preserving all three axial atom pairs requires
the rotation vector to vanish. A remaining translation has multipliers
`a.e_i`, so their simultaneous vanishing forces a=0.

For each of the six normalized Euclidean generators K_alpha (three
translations and three rotations about the local patch center), set

    v_alpha=rotation_part cross u0-(K_alpha.grad)u0,
    F_alpha(xi)=Omega(K_alpha,xi)
                     =-rho integral v_alpha.xi.

All are analytic. In any fixed open response ball the Gram of
`curl(-rho v_alpha)` is therefore positive definite: a zero combination
would extend analytically to a Euclidean symmetry of the seed. Choose
six mutually disjoint response balls away from all physical core jets and
raw cages. Their six Gram minima and all needed C^m bounds are strict
finite constants. A sufficiently small C^m neighborhood of the seed
preserves them. Source Proposition 3.8 gives this good-patch event strictly
positive probability. No observed soft eigenvalue is used.

For a concrete stationary locally finite selection, take a cubic grid of
spacing L>2R with a uniform random shift in its fundamental cell, independent
of u0; accept a site when its R-ball satisfies that good event. Rotate the
grid with the whole realization when forming the Haar mixture. The accepted
intensity is `n=p_good/L³>0`; supports in different accepted balls are
disjoint. Each accepted patch contains a genuine finite-core vortex
structure of the ACTUAL sampled Beltrami field, not the reference field.
The full source theorem independently guarantees positive-density knotted
invariant tori. No independent-field gluing at grid faces occurs.

## 4. Exact compact moment reconstruction on each actual patch

In response ball j form its local six-by-six Gram and its dual compact
generator eta^j, normalized by

    F_alpha(eta^j)=delta_alpha,j.

Explicitly, if f_alpha=-rho v_alpha and chi_j is its cutoff, use
`eta^j=sum_beta curl(chi_j curl f_beta) (G_j^-1)_beta,j`.
Because the six response supports are disjoint,
`Omega(eta^i,eta^j)=0` exactly. All these fields have compact vector
potentials and uniformly bounded norms on the good event.

The exact projector

    Pi xi=xi-sum_j eta^j F_j(xi)

kills all translation and rotation moments and preserves every physical
core jet. Choose raw body and circular internal cages and physical opposite
core rotations as in 0048, but evaluate them on THIS actual field. Then

    r=b eta^rotation_axis+Pi A_body,
    Q=Pi(Q_R+C1),  S=Pi C2

have zero translation moments, one normalized common rotational moment,
zero body/internal KKS crosses and the unchanged nonzero circular KKS
pairing. These are exact support and dual-response identities, not a
presumption that an arbitrary projection is symplectic. On time-reversed
partners choose b with the corresponding sign so the geometric generators
are the same, while the KKS signs reverse.

The translation constraints imply
`integral (xi cross omega0)=0` for every selected compact internal column.
Their compact-curl construction also gives `integral xi=0`. These facts
remove both uniform displacement and force means per patch. Consequently
their stationary sums have zero mean in each ergodic component, not just
after the global Haar mixture. They are appropriate internal directions
for the exact mean-material split in section 2.

The six-response construction also makes the proposed common KKS density
origin independent: for an axis e and patch center X_a,
`K_e(x)=e cross(x-X_a)+e cross X_a`; the second term pairs to zero by the
translation constraints. Whether its global symmetry produces the needed
independent physical common rotor, rather than merely reindexing a law,
is the separate explicit construction in 0059. The response rank alone
does not answer that physical question.

## 5. Complete stationary orbit Hessian, with all ambient interactions

Let xi be a stationary locally finite sum of the selected compact
generators with bounded second moments of their amplitudes. Its force is
`F=xi cross omega0`. Work in the Hilbert space of stationary vector fields
with inner product `E[f.g]`. The Leray multiplier is the ordinary orthogonal
matrix `P(k)=I-kk^T/|k|²` on nonzero spectral modes. It is an L² contraction.
Set the zero-mode velocity separately by the mean material coordinate; the
internal projection P0 removes that mode. For the uniform internal columns
in section 4 the force has zero mean component already.

The exact stationary isovortical variation is

    delta omega=curl F, delta u=P0 F,
    delta² omega=curl(xi cross curl F).

Stationarity gives expectation integration by parts:
`E[partial_i f]=0` when the derivative is integrable. Repeating the full
second-variation calculation, including delta² omega, gives

    H_dens(xi,zeta)
      =rho E[(P0 F_xi).(P0 F_zeta)
                              -F_xi.curl F_zeta/lambda],
    Omega_dens(xi,zeta)=rho E[omega0.(xi cross zeta)].

These are one-fluid quantities per unit TOTAL volume. The global Leray
term retains all interpatch kinetic crosses and all ambient/exterior
motion. It is NOT replaced by intensity times an isolated-patch energy.
For the local KKS and helicity terms, disjoint supports permit the exact
intensity identity `E[local density]=n E_Palm[integral_patch local density]`.
For the nonlocal kinetic term no such diagonal simplification is made.

The pressure/shared-face reaction is therefore already that of one smooth
global field. A diagnostic finite partition may expose internal tractions,
but they cancel pairwise; it does not create walls or discard the ambient
momentum. Exterior-supported cages belong to this same probability measure
and kinetic expectation, so they have neither uncounted fluid nor an
independently appended external mass.

## 6. Carrier positivity without dropping the nonlocal crosses

The stationary representation actually simplifies the positive bound. The
Leray part of H is nonnegative. A negative-helicity carrier on a good patch
has the exact leading local helicity contribution

    -lambda^-1 integral F.curl F
        =(|k|/|lambda|) integral(phi omega_z)²+O(1).

The O(1) bound is uniform on the good event. It includes derivatives of
omega, cutoff returns and all fixed-response/core attachments. Integration
by parts moves curl onto a fixed attachment in each cross, making that
cross bounded independently of k. Different compact patch supports have
no local helicity cross. No sign is imposed on their kinetic cross: the
full kinetic Gram is kept nonnegative as a whole.

Thus on any selected compact coordinate vector z, and equally on its
stationary square-integrable patch amplitudes,

    H_dens >=rho n[(|k|/|lambda|) A_*-C_*] E_Palm |z|².

The positive A_* and finite C_* are uniform strict good-event bounds.
For affine shear use the exact nine-bond frame in 0054, reconstructed on
each sample; for angle gradients use the specified first-difference cage
maps. Their local full-rank frame bounds supply A_*. A finite analytic
carrier threshold therefore gives a coercive complete compact Hamiltonian.

If the conjugates vary independently by patch, their momentum block is a
bounded positive OPERATOR, not an assumed diagonal array. Its inverse
exists by this coercivity. Eliminate it in the same action: the effective
mass is D^* P^-1 D and the stiffness is the full Schur complement. Their
exact definitions retain every Leray cross; local one-patch numerical
coefficients are not transplanted. Coherently prescribed physical fields
are the slow-affine closure, while conjugate momentum response is varied
before time-reversal averaging. This is not free nonaffine strain relaxation.

## 7. Physical currents and retained common-rotation dependency

For a tagged tube D the 0055 observable identity remains

    L_D=J_coherence+B_surface-J_exterior,
    B_eff=B_surface-J_exterior,
with the actual field and selected response coefficients recomputed here.
This is a tagged-tube observation identity; L_D alone is not asserted to
be the total intrinsic spin of all fluid. For the latter, form S_full from
the complete material partition, INCLUDING ambient parcel intrinsic spin.
If S_can denotes the internal canonical spin, set Delta=S_full-S_can.
The exact current improvement is then
`P_can=P_mean+(curl Delta)/2`, `S_can=S_full-Delta`, with its corresponding
pressure/current flux. Its tagged-tube contribution is the displayed
B_eff; any remaining ambient term stays in Delta rather than disappearing.
A canonical momentum is not silently identified with rho times a tube
centroid speed. The material mean density in section 2 is exactly rho.

Established here: a source-supported normalized stationary smooth-Euler
ensemble, a positive-density actual EPS core/cage reconstruction, exact
six-response compact moment projection, one-action material mean splitting,
and the complete exterior-inclusive compact Hessian density with a finite
positive carrier bound. This is the explicit assembly representation that
the single decaying-defect normalization did not supply. The independent
physical common-rotor reconstruction in 0059 and the parent's full slow
kinetic normal form are separate dependencies, not consequences asserted
from statistical invariance alone.

## 8. Selected simpler affine-relative representation

The parent and 0059 subsequently selected a simpler construction: the
target needs U and Phi, not an additional independent cyclic common rotor.
Keep the actual compact physical relative-angle pair q,S and define
`Phi=beta+q`, `beta=curl U/2`. The eleven affine moment projection of 0059
removes its KKS crosses with the three translations and eight tracefree
affine generators, with its actual core jets unchanged. It is reconstructed
on the same good patches by the same dual-response method; the seed
Fourier-atom independence extends from six Euclidean to all eleven affine
generators. The compact pair's positive I,K then gives

    L_internal=I |Phi_dot-beta_dot|²/2-K |Phi-beta|²/2

before the parent's vector/isotropic normalization. This has a genuine
physical relative core tilt and retains the complete gradient kinetic
cross; it uses no statistical common-rotation reindexing as an inertial
coordinate. The associated exact same-ensemble affine energy and finite
reaction bound are in `affine-spectral-energy.md`. The eleven moment
constraints also yield the C^2 slow-action construction in `slow-locality.md`.

This is an append-only representation change within the SAME U,Phi target,
not a reduced objective. The earlier six-response/global-rotation route is
preserved above as branch evidence, not used to supply an extra body mass.
