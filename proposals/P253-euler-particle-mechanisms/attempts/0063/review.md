# P253/0063 joined independent review of 0058 and 0062

## Frozen transaction and correction provenance

This transaction independently reviews target-owned attempts `0058` and
`0062`, both registered to `particle-foundations`. The reviewer authored or
implemented neither target, its formulas, verifiers, corrections, nor
receipts. The frozen README has SHA-256
`dd14e704c9e48d1789541bfdb97efee390d9fc6820b9002bb73cd69fc1c1357b`.
Central activation exited exactly `0`; the command, stdout, empty stderr, and
exit hashes are respectively
`0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The substantive pass found one bounded 0062 correction. Its final receipt has
SHA-256
`6a025c24d1f5c84fafae8a216a3d7a1e829df68b815687f933684e5eef56c169`.
The correction check was limited to the cylindrical Hodge kernel, collar
topology, finite rows, and their direct claim edges. No unchanged oracle was
rerun and no second substantive pass was performed.

The final corrected 0062 hashes are:

| Artifact | Corrected SHA-256 |
|---|---|
| `0062/derivation.md` | `f8e62434a0b52e37ff89bd5d57bcfdecc4859e7089ab6a64d8e9a4b90101b449` |
| `0062/result.yaml` | `d559a4c3b559c245be2994725aff57a6b36307daa6ad140ce379b9c9b5a82d9d` |
| `0062/source-audit.md` | `a6119cfe92d1d2ee05e639cd0938720a1082224716260425aadaac86c2ee738b` |
| `0062/validation.md` | `57cc0ed8e6f47894cdda0a8157d32258d34d2d7b1f3f29b4d89b46c0e6551f77` |

The 0058 target remained unchanged at its frozen hashes, including
`a19f0ac5b84cf1f2eaaf0e0db6757ca911b6e28bb883bd641db7a05f41e8cbdd`
for `derivation.md` and
`4409eaabee0e99d66a1cfd5b88b03b77408ddbae94165ccd4ef1a9fb9f110357`
for `result.yaml`.

## Source applicability

The cached Cao--Lai--Qin--Zhan--Zou paper, arXiv:2206.10165v2, has SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.
It supplies the fixed smooth axisymmetric compact-core translating ring and
its printed positive-part vorticity structure. It supplies no nonaxisymmetric
DA graph theorem, limiting absorption, distorted adjoint, mode, or rotating
branch.

Gallay--Smets arXiv:1805.05064v3 has SHA-256
`081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`.
Its straight infinite column confirms that eliminating velocity leaves a
coupled nonlocal critical-layer problem. Its compactness, source space,
spectral theorem, and boundary conditions do not transfer to the curved Cao
ring. The helical-patch comparator arXiv:2411.02055v1, SHA-256
`d6fbac9b2aa7f2394684807717224f31e3ae0afef8fed092dc7fa4def2223e4e`,
establishes no smooth same-Cao-leaf branch.

The controlling independent 0057 review keeps the local Cao cells, Piola and
pointwise Hodge identities, but leaves `HJ2`, `GR`, common DA graph closure,
and Riesz transfer unproved. The right-reduced KKS convention imported from
`0005/derivation.md`, SHA-256
`90e17874b928e6ccda19b9e745346836f1c4a2140ced8a3e434fbc44eb70abc2`,
has the sign used below. No source supplies the missing 0062 global operator
theorem.

Neither target adds an importable production API. The nearby frozen
`collective_coordinates.py` and `euler_twisted_carrier.py` modules concern
different conditional reductions and do not support the carrier-map,
critical-layer, or collar claims reviewed here. Attempts `0064` and `0065`
were not opened or used as evidence.

## Unit A: 0058 carrier map and exact route boundary

### Whole-space rotating equation and carrier map

For `W=v-c e_z-Omega R`, compact vorticity and whole-space decay give

    curl(W cross omega)=0  iff  [W,omega]=0.

On simply connected `R^3`, the curl-free compact field `W cross omega` is an
exact gradient. The velocity/Leray equation follows with the exterior
circulation row retained in any core/exterior decomposition. This is an exact
Euler representation, not a core-only pressure solve.

For a volume-preserving same-leaf map `omega_g=g_*omega_0`, naturality of the
bracket gives

    v_g-c e_z-Omega R=g_*Y,    [Y,omega_0]=0.

This carrier-map/centralizer equation is exact. With
`omega_0=zeta partial_theta`, `zeta>0`, the nonzero toroidal harmonic
centralizer matrix has determinant `i n^3 zeta^3`. Therefore its smooth
positive-core centralizer has no nonaxisymmetric component. At zero harmonic,
the poloidal streamfunction is locally `F(zeta)` on each connected regular
level region; boundary, zero-core, and exterior centralizers require separate
matching.

### Fourier obstruction and conditional pendulum

The exact two-wave calculation gives transverse output numerator `(1,-1,0)`
at a raw transport divisor equal to zero. It therefore refutes universal
polynomial divisibility of the Euler quadratic tensor by that divisor. It does
not evaluate the Cao source-specific numerator or refute transparency after a
full-Hodge distorted reduction.

The DA principal symbol is onto the transverse plane when
`k dot omega_0` is nonzero, but the full linearized operator also contains the
nonlocal Hodge polarization. Symbolic DA richness is not a bounded right
inverse or an adjoint trace. The physical reverser, rotation sign,
`J(a)=J(0)+(ell sigma/2)a^2+O(a^3)`, and the distinction between full SO(2)
and isotropy-shortened periods are exact kinematics.

The displayed pendulum width and frequency are correctly conditional on a
full-Hodge generalized adjoint, physical KKS measure, and defined nonzero
coefficient. They construct no Euler inner layer, rotating branch, or
restoring neighborhood.

**Unit A verdict: established at the stated carrier-map, positive-core
centralizer, Fourier nondivisibility, DA-symbol, and kinematic scope. The
pendulum continuation and rotating branch remain blocked by the named
full-Hodge and nonlinear matching constructions.**

## Unit B: 0062 full-Hodge coordinates and topology

### Coordinate block, Hodge kernel, and positive-core inverse

Direct coordinate-basis evaluation of
`A_0 eta=-[W_0,eta]-[B_R3 eta,omega_0]` gives all three displayed components
with the correct Hodge signs. The bounded correction fixes the explicit
whole-space kernel. If

    D(a)=diag(1,a,1),

then source coordinate components are converted by `R_vartheta D(r')`, and
the Cartesian result at `(r,0,z)` is converted back by `D(r)^(-1)`. For
`r=0`, this means the regular-axis limit of the physical harmonic field, not
uncontrolled division by `r`. The abstract Cartesian Biot--Savart operator,
coordinate Lie-bracket block, and KKS calculation are unchanged.

For fixed `n!=0` and `zeta>0`, the three-component map
`C_0 xi=-[xi,omega_0]` has the stated pointwise inverse, including the
`grad zeta` correction in the theta component. The divergence identity shows
that divergence-free vorticity recovers a divergence-free displacement. This
is a two-sided local result only where `zeta>=zeta_min>0`, with constants
depending on `1/(|n| zeta_min)`. It is not an unweighted free-boundary
isomorphism.

Jacobi and `[W_0,omega_0]=0` give the exact intertwiner

    A_0 C_0=C_0 M_0,
    M_0 xi=B_R3 C_0 xi+[xi,W_0].

A compact representative is obtained by cutting off outside the vorticity
support and repairing the divergence in the exterior collar. Character
orthogonality removes the rigid `0,+/-1` rows for fixed `|n|>=2`; other finite
same-character rows require finite-dimensional profile corrections.
Bogovskii repairs only divergence and is chosen to preserve already imposed
rows.

### Exact collar counterexample and topology separation

At integer `s=4`, `p>=6`, choose a fixed `|n|>=2` and an `h`-by-`h`
meridional boundary patch. A streamfunction of amplitude `h^5` produces an
exact divergence-free displacement of point size `h^4` with
`H^5` norm comparable to one. Positive-core injectivity, harmonic
orthogonality, finite-row removal, and support separation from exterior
stabilizers give a uniform lower bound in `H^5/ker C_0`.

The exact positive-part representation gives
`|D^j zeta|<=C_j d^(p-j)` for `0<=j<=5`; no derivative is inferred merely
from `zeta=O(d^p)`. The scaled product ledger then gives

    ||C_0 xi_h||_(H^3)=O(h^(p+1)),
    ||C_0 xi_h||_(H^4)=O(h^p).

Because `C_0 xi_h=curl(xi_h cross omega_0)`, the exact identity
`B_R3 C_0 xi_h=P_L(xi_h cross omega_0)` supplies the Hodge estimate and
resolves the compact-curl zero-mean/low-frequency row. It does not invoke a
false unconditional inhomogeneous `H^3`-to-`H^4` Biot--Savart bound. Smooth
boundary tangency gives `W_0^normal=O(h)`, and direct substitution yields

    ||B_R3 C_0 xi_h||_(H^4)=O(h^(p+1)),
    ||A_0 C_0 xi_h||_(H^3)=O(h^p).

Thus both the ordinary ambient-vorticity lower bound and the ordinary
generator-graph lower bound fail to control the displacement quotient. This
is an actual counterexample at the stated integer regularity, not an expected
estimate. It leaves alive three distinct propositions: the interior inverse
with `zeta_min>0`, the explicit weighted orbit norm, and a source-specific
ambient-vorticity sandwiched trace.

### KKS weight, centralizer topology, and layer exclusions

The right-reduced KKS form gives

    Omega_KKS(C_0 xi_n,C_0 chi_-n)
      =2 pi rho_0 integral r^2 zeta
         (xi_n^z chi_-n^r-xi_n^r chi_-n^z) dr dz.

The sign, density, toroidal Jacobian, and weight are correct. This measure
does not create or normalize a spectral mode.

On each connected regular positive-core level annulus, a smooth centralizer
has poloidal streamfunction `F(zeta)`, preserves every corresponding torus,
and has an integrable rigid return on its angular label. Smooth carrier-map
pushforward preserves orbit and separatrix incidence. Therefore a spatial
pendulum layer that creates a center, saddle, and separatrix at such a regular
level is incompatible unless transparency holds, the allowed centralizer
already has matching critical topology, or the construction leaves the smooth
DA class. This does not refute a coefficient-space normal form, a compatible
critical-level construction, or a smooth integrable nonaxisymmetric branch.

For a nontrivial self-similar layer of width `h=sqrt(epsilon)` and physical
amplitude `Theta(h)`, the displayed `H^s`, `H^(s+1)`, and induced
`H^(s-1)` exponents are correct and diverge for `s>7/2`. The verdict applies
only to that nontrivial amplitude-normalized self-similar family. It does not
exclude zero, `o(h)`, global smooth, differently scaled, or separately
licensed weaker-regularity layers.

### Full-Hodge trace boundary

The Sokhotski formula, slope-weighted raw delta trace, and formal distorted
Lippmann--Schwinger adjoint equation identify the correct next operator
problem. They do not prove a limiting-absorption boundary value or construct
its adjoint. Both `q_*` and `V_*` are therefore undefined, not zero.

The next theorem can start with one explicitly constructed smooth DA seed in
one fixed nonzero-`n` block and its first source-bearing resonance. It needs
only the source-specific sandwiched core/interface/exterior Hodge or Grushin
trace and finite rows for that witness. A monolithic all-sector theorem is not
a prerequisite for this first trace, although later nonlinear convergence may
require broader control.

**Unit B verdict: established for the corrected coordinate differential
block and kernel representation, positive-core inverse, intertwiner, KKS
weight, connected-regular-torus topology, exact unweighted collar
counterexample, and scoped layer exclusions. The global graph/LAP adjoint,
`q_*`, `V_*`, and rotating branch remain blocked by the named construction.**

## Oracle scope and final joined verdict

The frozen 0058 verifier and stdout have SHA-256 values
`8bb56725c86568277b0f1f57aceb54ff04f978ccecd7953eec41dc8948357881`
and `3715973336f135eda2bd115adb9eaa68654b12e28d5c4a0f09ebd421c2762124`.
The frozen 0062 verifier and stdout have SHA-256 values
`cd3c2fcc4ab612f73d7b393b4885bc523316e977184aa9b7e51efae0194d62ad`
and `44bf4ba1c7fc71bb27dd55bc678d7f1b083987557eebbc360b7fa4aef421b14a`.
Both recorded executions exited `0`.

Those exact predicates support the Fourier numerator, centralizer
determinant, coordinate Lie brackets, positive-core inverse, KKS density,
volume-generator divergence, and self-similar powers. They do not test the
cylindrical kernel conversion, collar quotient/graph estimates, global
operator domain, limiting absorption, or a physical mode. The corrected
claims are instead supported by the explicit analytic derivations above, so
the unchanged oracles were correctly not rerun.

The single bounded correction fully closes the review findings. The joined
result is established as an exact same-field representation and topology
boundary: the carrier map and local Hodge/DA/KKS algebra are exact; universal
raw divisibility, an unweighted free-boundary inverse/graph equivalence, and a
generic topology-creating smooth pendulum repair are unavailable for the
specified mechanisms. This is not a Cao branch no-go. The weighted-orbit and
source-specific sandwiched fixed-harmonic routes remain live.

The next dependency is one smooth fixed-`n` physical DA seed and its first
source-bearing resonance, followed by the corresponding sandwiched whole-space
Hodge/Grushin adjoint trace and KKS-normalized coefficient. No rotating branch,
stability theorem, P2 completion, particle, quantum, or relativistic conclusion
is licensed.
