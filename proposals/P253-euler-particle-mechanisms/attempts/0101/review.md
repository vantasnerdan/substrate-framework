# P253/0101 independent final review of P253/0100

Reviewer: `particle-balance-review`

Target owner: `root`

Activated review README SHA-256: `a3325beedc03c12f12fcc2d98811ddc68acbe1ccb70e3af68427e3fc1a2145cb`

Frozen pre-review target manifest SHA-256: `dff6ddccd1166019720f6a545344d246195075f9cf462833ea68886c6af8d54f`

Post-correction target manifest SHA-256: `45baf6fce79b861f4dbde1543a92e1ead0161d88ab2547a11439e23eba63666c`

## Boundary and provenance

This was a content-blind, independent, non-author review of final P253/0100.
The activation receipt contained exactly `0`, the activated README matched its
frozen hash, and every frozen target selector matched before body opening.
P253/0095 was not opened or used. No parent, P6, particle, electron or neutrino
conclusion was adjudicated.

The manifest uses attempt-relative paths for `0100` artifacts and
repository-root paths for the module and test. The final exact-v5 receipt has
eleven PASS lines, empty stderr and exit zero; the focused-v5 receipt has twelve
passing tests, empty stderr and exit zero. These receipts were reused without
rerunning unchanged evidence.

The substantive pass requested one bounded evidence-ledger correction. The
old continuation receipt called the Cao magnetization superpotential
unconditionally nonzero even though the construction proves its formula but
not nonvanishing of `grad Phi cross grad P` on the actual branch. The corrected
receipt now says `potentially nonzero`. Its SHA-256 is
`e44fcde0af9cb45217045cd0a94ccfcb0ca5c2f99f5c834f7cbebe123154e978`;
the correction receipt is
`5473d007777030b6c4f6c3641808dd3a1fa5da2d8988a24d48717065c1636fc2`;
the refreshed completion receipt is
`49711c3a0c6746f10a98c3af8eebaf9d839690a8a4be90c006033c53e51b4619`.
The correction-only check verified these hashes and the affected prose. No
equation, API, test, verifier, or captured output changed.

## Source and authority audit

The accepted inputs are consumed only for a smooth transported scalar, the
charged axisymmetric no-swirl Cao conventions, and the force law at their
declared scopes. The load-bearing Ertel, force, monodromy, cylindrical and
punctured-domain identities are rederived in `0100`; no author-stage motivation
is treated as authority. The derivation does not import a Lorentz current,
compact character, quantized action, force-selected lock, moving defect, or
persistent nonaxisymmetric carrier.

## A0 — unforced Ertel invariant and current

For a classical incompressible Euler solution and a transported true scalar,

    D_t omega=(grad u)omega,
    D_t grad chi=-(grad u)^T grad chi.

Consequently

    D_t(omega dot grad chi)
      =((grad u)omega) dot grad chi
       -omega dot ((grad u)^T grad chi)=0.

With `div u=0`, this is
`partial_t q_E+div(q_E u)=0`. A true scalar has polar gradient, while
vorticity is axial, so `q_E` is a pseudoscalar and `q_E u` is axial. This is a
Galilean continuity-current tuple, not a Lorentz four-current.

Since `div omega=0`, `q_E=div(chi omega)`. Thus its integral is exactly the
oriented boundary flux of `chi omega`; it vanishes for smooth globally
single-valued whole-space data when that flux tends to zero, including compact
support. The local identity only needs the displayed classical derivatives;
the integrated claim additionally needs trace/integrability and vanishing
boundary or infinity flux. Topology alone does not defeat Stokes while
`chi omega` remains a global smooth field.

**A0 verdict: established as stated.**

## A1 — forced Ertel source and flux

For acceleration convention
`D_t u=-grad p/rho_m+f`, curl gives

    D_t omega=(omega dot grad)u+curl f.

The same stretching cancellation yields

    D_t q_E=(curl f) dot grad chi
           =div(f cross grad chi).

Therefore the conserved local flux is exactly
`q_E u-f cross grad chi`; its minus sign is correct. Integrated conservation
again requires the corresponding total boundary flux to vanish. The force
does not change the assumption that `chi` is transported.

**A1 verdict: established as stated.**

## B0 — closed-line real constant lock

If `q_E=lambda chi`, then on a vorticity line parameterized by
`dx/dtau=omega`,

    d chi/dtau=lambda chi.

For a closed line carrying nonzero single-valued `chi`, return to the same point
requires `exp(integral lambda dtau)=1`. For real `lambda`, its loop integral is
zero. Hence a real nonzero constant is impossible on that line. The mechanism
does not apply to zero tags, open lines, complex monodromy, defects, or a
variable real coefficient with zero loop integral.

**B0 verdict: refuted only for a nonzero real constant lock on a nonzero-tagged
closed vorticity line.**

## B1 — exact forced lock condition

In unforced Euler, differentiating `q_E-lambda chi=0` gives
`chi D_t lambda=0`, so `D_t lambda=0` wherever `chi` is nonzero. Under forcing,
the correct and distinct compatibility equation is

    chi D_t lambda=(curl f) dot grad chi.

It remains meaningful without division on the zero set: there the source must
vanish if the lock is to persist. The closed-line monodromy is a simultaneous
kinematic restriction. The parity ledger also requires a pseudoscalar
`lambda` when `chi` is a true scalar.

**B1 verdict: established as an exact compatibility equation.**

## B2 — force-selected variable lock

The compatibility equation does not select a force, an initial coefficient,
or a globally regular solution satisfying its zero-set and closed-line rows.
Nor does it supply a material-leaf action or electron-linked endpoint
interaction. Those are explicitly left open rather than inferred from B1.

**B2 verdict: blocked at a same-carrier force/coefficient solution and its
action/interaction map.**

## C0 — regular Cao density and convective current

In the frozen cylindrical convention,
`omega=r zeta e_theta=zeta partial_theta`. Since
`grad chi` has azimuthal component `(1/r)partial_theta chi`,

    omega dot grad chi=zeta partial_theta chi.

Every smooth axisymmetric tag, including `chi=F(I)` on the regular tagged
band, therefore gives `q_E=0` pointwise. The convective part `q_E u` also
vanishes pointwise. This is a geometric statement about that tag/carrier, not
a universal current no-go.

**C0 verdict: the zero is established exactly; the regular Ertel-density
supplier route is refuted on this axisymmetric no-swirl geometry.**

## C1 — charged Cao magnetization superpotential

With

    f=(g chi_c/rho_m)(-grad Phi-H grad P),
    grad chi_c=chi_c'(P)grad P,

the nonconvective forced flux is

    -f cross grad chi_c
      =(g/rho_m)chi_c chi_c' grad Phi cross grad P
      =curl[(g Phi/(2rho_m))grad(chi_c^2)].

The sign, density factor and factor `1/2` are correct. For constant `g` and
`rho_m`, its divergence vanishes identically. Two poloidal gradients give a
toroidal axial current. It can be locally nonzero, but the target does not
prove that `grad Phi cross grad P` is nonzero on the actual charged Cao branch;
the final wording now correctly says *potentially* nonzero. Its integrated
content is a boundary row and it is not a net charge or weak current.

**C1 verdict: established as an exact local superpotential formula, with
actual branch nonvanishing unproved and unclaimed.**

## D1 — punctured-domain equilibrium phase

On `M=R^3` minus the symmetry axis,

    dtheta=(-y dx+x dy)/(x^2+y^2)

is a smooth global closed non-exact one-form. Its Euclidean contraction with
`omega=zeta partial_theta` is `q_theta=zeta`. Because the compact toroidal
vorticity support stays away from the axis,

    integral_M q_theta dx
      =2 pi integral zeta(r,z) r dr dz
      =2 pi kappa.

This does not use a global real lift of `theta`, so it does not contradict the
exact-divergence result in A0. For the fixed geometric azimuth,
`D_t theta=u_theta/r`; hence it is material at the axisymmetric no-swirl
equilibrium, and along an axisymmetric no-swirl evolution only while support
separation and inner/outer flux hypotheses persist. The target correctly stops
short of a general perturbation theorem and distinguishes this continuous
circulation-derived quantity from electromagnetic charge and internal action.

**D1 verdict: established at the stated punctured-domain equilibrium scope.**

## D2 — moving defect and nonaxisymmetric continuation

No independently advected circle phase, moving defect line/domain,
moving-inner-boundary flux theorem, persistent nonaxisymmetric carrier,
finite-energy defect action, or same-carrier endpoint interaction is
constructed. The fixed-angle calculation in D1 cannot provide those objects.

**D2 verdict: blocked at the moving-domain material phase and persistent
same-carrier construction.**

## Evidence attribution and combined verdict

The exact verifier genuinely checks the arbitrary-matrix stretching
cancellation, product/divergence identities, forced sign, parity table,
algebraic lock residuals, cylindrical contraction, superpotential curl and
fixed-angle derivative. Some global rows are necessarily analytic rather than
executable: the boundary theorem, nonzero-tag closed-line hypothesis, global
punctured-domain patching and absence of a moving defect. The API's
`ertel_charge_from_flux` and charge helpers are typed algebraic interfaces, not
independent proofs of those global facts. The focused suite is regression
evidence for those interfaces. Repository validation is workflow evidence.

After the single bounded correction, P253/0100 is **established at its exact
route-scoped supplier boundary**. It supplies the unforced and forced local
Ertel currents, the exact closed-line obstruction and forced compatibility
condition, the regular Cao zero, a potentially nonzero Cao superpotential
flux, and the punctured equilibrium integral `2 pi kappa`. It does not supply
force selection, a general advected defect, action quantization, a Lorentz
chiral current, a particle, or P6.

The strongest next construction is either (i) a force-selected variable lock
satisfying the zero-set, closed-line, parity and action rows on one persistent
carrier, or (ii) an independently advected circle phase with a moving defect
domain, complete inner/outer flux proof, finite action and endpoint interaction.
