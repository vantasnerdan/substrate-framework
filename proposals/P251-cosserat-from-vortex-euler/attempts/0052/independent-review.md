# Independent review: material centering and boundary-current joining

Reviewer: `/root/smooth_core_review`, distinct from the author of 0052.
Date: 2026-09-05. One bounded scientific review under AGENTS.md and the
physics/small-ratio skills. The frozen attachment scope is the material
decomposition, response projection and current identities, with the
remaining centered common-momentum nondegeneracy explicitly separate.

## Decision and supported result

The material/action and boundary-current bridge is established as stated.
The fixed-response projection is also an established construction on its
declared compact-curl class. It removes the exact physical centroid
momentum response, retains the physical jets, and identifies the precise
additional residual needed to construct a centered common-momentum partner.
No load-bearing correction is requested for these statements.

The strict inequality Delta_common>0 and the nondegenerate centered
internal KKS pair are not proved by 0052. This is an explicit construction
frontier, not a hidden assumption in the established projection identities.
Neither the 17-check receipt nor this review supplies those missing
positive data or certifies the parent's full common-action closure.

## Evidence and material action

I read the README, complete `material-joining.md`, verifier and repaired
17/17 receipt. The verifier uses the shared orbit API and preserves the
initial import-path failure separately. Exact material-map differentiation,
response duality and current algebra are the strongest practical oracles.
The existing receipt was inspected and reused without a redundant run.

Writing g=X+h with its mass-weighted mean constraint gives
`integral rho hdot=0`. Expansion of the full material kinetic energy then
removes the entire Xdot--hdot cross term, giving the physical centroid
mass M and internal kinetic term once each. The canonical momenta and
Hamiltonian displayed in the proof are its exact Legendre transform on
that constrained slice. Shared-face matching remains a constraint with
the Euler pressure reaction; dropping it would instead describe separate
free parcels. Variation therefore reproduces the material force and
intrinsic-spin balances of 0051 without an inserted spring or mass.

This canonical material identity does not license appending an additional
centroid mass to an already uncentered orbit Hamiltonian. The proof
explicitly avoids that double counting.

## Moving material mean and its energy Gram

For the stationary invariant torus, the divergence theorem applied to
x_i u0 gives zero integrated velocity because u0 is solenoidal and tangent
to the boundary. The reference centroid velocity is thus exactly zero.
Differentiating its velocity mean along a geometric orbit variation gives
both the volume velocity variation and the boundary displacement term:

```
C_D(xi)=|D|^-1 [integral_D v_xi
                          +integral_boundary u0 (xi.n)].
```

No denominator variation occurs for a volume-preserving deformation.
Rotating the entire field and domain about their mass centroid rotates
their zero integrated velocity into zero. Hence C_D(K)=0, even though the
fixed-domain integral of vK alone need not vanish. Compact directions
inside D instead have no moving-boundary term.

Since the base centroid velocity is zero, differentiating its energy
twice gives exactly `M C_D(xi).C_D(eta)`, a positive Gram of rank at most
three. Subtracting this once removes the finite cell's centroid energy
from the full orbit accounting; it is not a claim that the remainder
contains no exterior fluid. Cauchy--Schwarz gives the stated bound by
`rho integral_D |v_xi|²` for compact directions. The cage velocity norms
are bounded while their complete orbit Hessians grow, so this finite
Gram subtraction has a bounded cost. Its K row is zero by the correct
moving-domain identity.

The proof appropriately does not infer a symplectic reduction from this
energy subtraction alone. Nonzero physical mean response requires either
the actual centroid/shape connection or an actual tangent projection;
the unchanged old KKS matrix is not licensed by the energy calculation.

## Fixed-response projection and topology

Self-adjointness of the ambient Leray projector gives
`integral_D v_xi,i=integral xi.(omega0 cross P(1_D e_i))`.
The response w_i is harmonic and curl free in the interior of D; with
analytic Beltrami omega0, curl f_i is analytic there. These are local
interior statements despite the distributional boundary source of 1_D.

For a nonzero smooth cutoff chi positive on an interior open subset,
integration by parts gives
`c_i(curl(chi curl f_j))=integral chi curl f_i.curl f_j`.
Thus W is the exact positive-semidefinite response Gram. If a coefficient
vector is in its nullspace, the corresponding analytic curl f vanishes
on that subset, then throughout connected D. On the stated class
xi=curl A with compact smooth A, its response is consequently zero by
another integration by parts. Conversely, vanishing on that test class
forces curl f=0. Therefore c(xi) lies in range W and
`c-W W^+ c=0`, without assuming rank three.

The compact-potential class is material: a curl-free field on a solid
torus need not be globally exact, and arbitrary flux-carrying solenoidal
directions would require additional cohomology responses. The proof
expressly confines this construction to its actual compact-curl
generators. The eta corrections preserve core jets because they vanish
near those jets. Their coefficients are bounded by the physical velocity
mean estimate for the selected high-frequency families. Cross energy
against the fixed smooth corrections is bounded by moving curl to the
fixed field, as in the earlier cage construction.

A compact solenoidal displacement has zero volume integral, hence zero
centroid displacement. Combining it with the projected zero centroid
momentum response places the tangent in the material canonical
`delta X=delta P=0` slice. Its canonical centroid block then has no cross
pairing with that slice. This does not establish that the old selected
internal KKS pair is unchanged: H and KKS must be recomputed after this
generally nonsymplectic projection, exactly as the proof states.

Delta_common is the squared norm of the common response after orthogonal
projection off the mean-response span. It is nonnegative, and is strictly
positive precisely when that functional is independent on the selected
class. Nonzero original common moment alone does not prove this; the
verifier's rank-deficient and parallel-functional cases expose that error.

## Physical spin, current improvement and admissible mass scope

The integrated curl identity gives `L_D=J_D+B_D` with the boundary
orientation `n cross u`. For a compact rearrangement inside D, the bulk
impulse variation agrees with the common moment of 0048 about the same
axis and center. The physical spin still includes B_D. A common rotation
of the field and domain leaves its axis component invariant; no finite
cell boundary term is discarded by appealing to whole-space decay.

With `A_ij=epsilon_ijk B_k/2`, direct contraction gives ax(A)=B and
div A=curl B/2. Setting P_can=P_centroid+div A and
S_can=S_physical-B, while adding partial_t A to the force flux, preserves
both balances exactly. The full moving-center version uses total force
and couple fluxes, retaining the convective terms. Product differentiation
gives the displayed angular-current improvement
`x cross P_can+S_can=x cross P_centroid+S_physical+div(x cross A)`.
Mass continuity is unchanged because div curl B=0. The new momentum is
not generically rho times centroid velocity; the explicit inverse current
map is necessary for physical observables.

For compact momentum tangents, the moving-boundary surface variation
vanishes and the row t is their actual boundary-spin response. Under
time reversal t and D reverse sign while P does not. Independent
elimination therefore cancels static odd terms and retains
`t P^-1 D (Bdot,qdot)^T`. This axial response is not removed by spatial
parity. The script's independent plus/minus momenta verify this ordering.

Finally the exact density rho block belongs to a finite material parcel
with exterior reaction and shared-face constraints, or to a specified
coherent energy-per-volume assembly. A uniform boost of the whole R³
EPS relative orbit has infinite kinetic energy and is explicitly not
used. The material decomposition by itself proves no global assembly or
finite spatial return for the parent's slow fields.

## Frozen hashes and disposition

- `material-joining.md`: `23c04af13360d343b551353f25af9836504ea091c02f64809d3e94e7ade05e1d`
- `verify.py`: `58f3d96f14f6895934a179fb52ec8f27efe413d3a48e4d964e68934920cc6515`
- `stdout.txt`: `cf1c1a65a6f4a782c5cca3430d75c93f8a6ffe46dee393a7e539fdc83d020cb0`

All validation is exact identity, finite-norm control or an explicitly
conditional Gram construction. No sampled small eigenvalue is used to
assert positive Delta_common. Acceptance is recommended at the stated
attachment scope, with the new centered-response construction separate.
