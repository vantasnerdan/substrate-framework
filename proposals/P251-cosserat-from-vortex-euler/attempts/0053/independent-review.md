# Independent review: isotropic average and kinetic normal form

Reviewer: `/root/smooth_core_review`, distinct from the author of 0053 and
its APIs. Date: 2026-09-05. One bounded scientific/code review under
AGENTS.md and the physics/small-ratio skills. Only the new tensor average,
normal-form theorem, API additions and five new test cases are reviewed;
the previous rotor-reduction functions are not reopened.

## Decision

Established as stated: the simultaneous SO(3) average gives the displayed
isotropic fourth-order tensor and coercivity conditions. The displayed
near-identity physical field map normalizes the assumed kinetic pencil
through spatial order two and, when applied to the potential as well,
produces the stated curvature correction. The implementation agrees with
these formulas. No load-bearing correction is requested.

This is conditional algebra on complete reduced action coefficients. It
neither proves their microscopic provenance nor establishes assembly of
the local cells into a common Euler continuum. Its transformed fields
coincide with the physical fields for uniform translations and rotations,
but not for general nonuniform configurations. Those premises and limits
are explicit and remain part of the positive theorem.

## Evidence and exact isotropic tensor

I read the README, complete `normal-form.md`, the new API definitions and
new tests, and the saved nine-test passing receipt. The five added pytest
cases include both helicities as distinct cases. The existing receipt was
reused; exact invariant contraction and direct pullback are the strongest
oracles, not a sampled spectral calculation.

The averaged tensor is symmetric in its i,k and j,l pairs. Isotropy
therefore gives precisely the two displayed invariant coefficients.
Contracting i=k and j=l gives T, hence `9A+6B=T`. Contracting i=j and
k=l gives L, hence `3A+12B=L`. Solving yields
`A=(2T-L)/15`, `B=(3L-T)/30`. The simultaneous rotation of C is essential;
averaging its orientation independently from the axis is a different
ensemble. The test's exact spherical moments independently recover all
81 components for the symbolic axisymmetric input and expose the
independent-orientation mutation. The invariant argument, not that test
specialization, establishes the general positive symmetric C case.

For G, decomposition into symmetric and skew parts uses
`tr(G²)=||sym G||²-||skew G||²`. It gives the API's trace, symmetric and
skew moduli without assuming that the trace coefficient alone is positive.
Positive C and a unit axis imply `0<L<T`, so c_s, c_a and
`3 c_tr+c_s` are strictly positive. The last is the trace-sector condition
after decomposing sym G into its trace and deviatoric parts. The test with
negative c_tr verifies that a spurious stronger sign condition is not used.

For periodic fields the integrated identity
`integral tr(G²)=integral (div Phi)²` gives transverse and longitudinal
bulk coefficients `nu A` and `nu(A+2B)`. A bounded domain must retain
the corresponding boundary current. Both coefficients are positive from
the same inequalities, but this replacement is not a pointwise deletion
of the null-Lagrangian boundary term.

The second spherical moment gives `j=nu J_Psi/3`. Since
`curl U-2 Phi=-2(Phi-beta)`, the restoring normalization is
`alpha=nu K_Psi/12`, not a separately selected coefficient. Similarly the
cage kinetic term is `nu J_beta |curl Udot|²/24`. It remains a gradient
inertia term until the next, explicitly displayed change of fields.

## Same-map kinetic and potential pullback

The input assumes the specified time-even, reflection-paired isotropic
two-field pencil and positive rho,j,alpha. Other gradient coefficients
are allowed to be either sign. Positivity of the zero-gradient mass
guarantees invertibility in a sufficiently small neighborhood of k=0;
no claim of all-k positivity follows from these input assumptions.

In a curl helicity let `d=m_Phi-b²/rho`. The off-diagonal map entry
`-b h k/rho` cancels the O(k) kinetic cross term. Its contribution to
the spin diagonal is `-b² k²/rho`; the spin rescaling removes exactly
the remaining d k². The translation rescaling similarly removes its
m_U k². Thus the mass pullback is diag(rho,j) coefficientwise through
k², as the API claims.

The same off-diagonal map applied to the two potential cross terms adds
`4 alpha b k²/rho` to the spin diagonal. The spin rescaling of the
constant locking term adds `-4 alpha d k²/j`. The other potential
entries retain their displayed coefficients through k². Consequently

```
C_eff=C+4 alpha b/rho-4 alpha(m_Phi-b²/rho)/j.
```

Both helicity tests independently multiply the actual mass and potential
matrices by the returned map and inspect every coefficient through k².
The wrong-sign kinetic map and omitted potential-correction mutations
are exposed. Neither normalization of frequencies alone nor an assertion
that gradient inertia vanishes would pass these checks.

The longitudinal spin sector lacks the curl cross term and gives
`C_L,eff=C_L-4 alpha m_Phi,L/j` by the same scalar rescaling. This
statement does not assert that the corrected curvature must be positive;
its sign is additional information about the actual input action.

The map has identity zeroth-order part and a formal low-gradient inverse.
Its differential form uses curl and second derivatives, which vanish
on the constant spin and affine displacement of a uniform rigid rotation,
as well as on a uniform translation. Nonuniform centroid, spin and
boundary observables must instead use the full returned physical map.
The code explicitly rejects known j=0 and directs that structure-free
branch to the unreduced action; no singular normal-form limit is claimed.

## API boundaries, hashes and disposition

The additions reject known nonfinite/nonreal input coefficients, invalid
axes, nonpositive density and inadmissible helicity. Unknown symbolic
signs retain the documented caller hypotheses. The argument has no
numerical soft-mode or discretization remainder: exact tensor identities,
matrix coefficients and strict input bounds satisfy its applicable
small-ratio verification requirements.

- `normal-form.md`: `9defd9d4d5c6180e0286548c02a287af67abca868050ee8a12391e915d7f49e6`
- `src/substrate_framework/euler_orbit.py`: `3c54b2e0c4b48c61cb45ad721b2298dd55026a1436da99308c1825374c758155`
- `tests/test_euler_orbit.py`: `34db4c92ea61f4c7adcd1b0b3aa81c9be8c923301acdd43a53c77022747cdc5d`

Acceptance of this conditional exact isotropic and second-gradient
normal-form result is recommended. The raw microscopic coefficient and
parent common-action closure obligations remain distinct from this review.
