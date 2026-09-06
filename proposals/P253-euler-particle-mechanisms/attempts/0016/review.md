# P253/0016 independent joined 0006/0012 audit

## Review boundary and independence

This is one non-author review of the joined same-family two-ring claim, not a
separate re-adjudication of attempts 0006 and 0012. The centrally registered
activation receipt is zero. The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`. The reviewer authored and
implemented neither source attempt, either exact verifier, nor their receipts.

The frozen 0016 README has SHA-256
`9222b3f38130ce4140edf15569b5a4413446ce9374087aec25cb4e972a219ead`.
The load-bearing reviewed artifacts are:

- 0006 construction:
  `c0ec3a813dd1a0df395552584b7e4b2858f4d486d1033684b77d684e1c79e7c5`;
- 0006 source receipt:
  `a8e745a436f144620295cba7fff1222ca8ff12bdda71b3f0ae7002f7dd83d071`;
- 0006 exact verifier:
  `ad9c74d15dab16c9d41d707069485ffe3783536d81e1b904ce1ee56eee1feb11`;
- 0006 final recorded output:
  `47c4ab82e65f9ceb44449732c52c48e66170355f5b94c2ff0456e2332fd38d60`;
- 0012 construction:
  `89d7f61eeefb38ece761d5823374c6065c637a54bcc0c94a4648ba5ca3276da3`;
- 0012 source receipt:
  `6ed6776f2354f0e0c18edd8fedfcb1117e0eeae484b1d1af57b70b7ba21c9764`;
- 0012 exact verifier:
  `bd2a9023da6c88d7ea6d7e97139ad088db83324c508e447f59ea8230042ef22d`;
- 0012 recorded output:
  `91093dc06de8e483bf7544f2b47c989e7135110f6ab95e5337ab0d9d2590ee82`.

The cached primary bodies match the source receipts: Cao et al.
`arXiv:2206.10165v2` PDF
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`
and text
`aa808fcf435a307c40815d4fdf81d4711c4a2750275e9810fbfc9b501f6085fd`;
Buttà--Marchioro `arXiv:1904.04785v2` PDF
`bdbc75defd6d5fb230e6784f0bd9566bbac69d2f4e831045d0a56478a86f7fb4`
and text
`f31552d96296eaaff360893835e2c993f111f9fad24ff534e3825c5b2ae77146`.
The latter is the corrected-v2 body required by the freeze.

## One full nonlinear Euler datum

The two labels are not independently evolved carriers. In the half-plane
`Pi={(r,z):r>0}` they obey

    partial_t xi_j+(K[xi_1]+K[xi_2]).grad xi_j=0,
    bold_omega=r(xi_1+xi_2)e_theta.

Thus their sum is one axisymmetric no-swirl Euler vorticity and both labels are
transported by its one common volume-preserving flow. Smooth nonnegative
compactly supported Cao data, and their small smooth weighted-volume-preserving
deformations away from the axis, remain in the natural classical
axisymmetric-no-swirl class used by the two sources. Injectivity preserves the
material labeling and initial disjointness until the support estimates below;
it does not turn either label into an autonomous solution.

The conventions are Cao's relative vorticity `xi=omega_theta/r`, measure
`dnu=r dr dz`, and energy/impulse after division by `2*pi`,

    E(xi)=(1/2) integral xi G_1 xi dnu,
    P(xi)=(1/2) integral r^2 xi dnu.

Buttà--Marchioro's intensity is consequently
`integral omega_theta dr dz=integral xi dnu`. A constant physical mass density
multiplies all physical energies, impulses, deficits, and the modulus by the
same positive factor and does not enter the vorticity evolution or alter any
sign or strict inequality. No variable-density term is silently omitted.

## Exact deficit and moment calculus

For one scaled or unscaled Cao carrier, let

    F(xi)=E(xi)-c P(xi),
    S=sup F,
    Delta_j(t)=S-F(xi_j(t)).

Each common-flow label remains on the exact rearrangement leaf, so
`Delta_j(t)>=0`. Individual impulse is not assumed conserved. Direct weighted
transport and axial kernel invariance give

    dot P_1=-double_integral xi_1 partial_zx G_1 xi_2 dnu_x dnu_y,
    dot P_2=-dot P_1.

Writing `E(xi_1+xi_2)=E_1+E_2+E_12` and using conservation of total energy and
total impulse yields, with no center ansatz,

    Delta_1(t)+Delta_2(t)
      =delta_0+E_12(t)-E_12(0).                         (1)

The signs in (1) are fixed: replacing `-E_12(0)` by `+E_12(0)` contradicts the
identity already at `t=0` for nonzero interaction.

If the axial support gap is `ell>0`, Cao's positive Green bound at exponent
`3/2` gives

    E_12 <= (4 C_G/ell^3) P_1 P_2
          <= C_G P_T^2/ell^3,                          (2)

where `P_T=P_1+P_2` is the actual conserved total impulse. This keeps all
exchanged individual moments and the correct factor four.

For the direct support bootstrap, equimeasurability and disjoint labeling give

    ||bold_omega||_1 <= 4*pi*sqrt(kappa P_T),
    ||bold_omega||_infinity <= 2 R_in Lambda

while `r<=2R_in`. The optimized three-dimensional Biot--Savart split therefore
gives

    U_A=C_BS(4*pi*sqrt(kappa P_T))^(1/3)
             (2 R_in Lambda)^(2/3),
    C_BS=3/(2^(4/3) pi^(1/3)),

and

    T_A=min(R_in/U_A,g_0/(4U_A)).

On `[0,T_A]`, `r<=2R_in` and the axial gap is at least `g_0/2`. Equations
(1)--(2) then imply

    0<=Delta_j(t)<=
      delta_0+8 C_G P_T^2/g_0^3-E_12(0).                (3)

The negative initial cross term in (3) is correct and useful.

For fixed radial support, Cao's maximizing-sequence compactness proves the
qualitative modulus

    m_R(eta)=inf{S-F(xi): xi in R(xi_*),
                 supp xi subset {r<=R}, d_R(xi,Sigma)>=eta}>0.   (4)

Indeed a zero-deficit sequence converges strongly in `L2(dnu)` modulo axial
translation; equimeasurability plus finite support upgrades this to `L1`, and
`P(|difference|)<=R^2 ||difference||_1/2`. This is a compactness implication,
not a coercive Hessian, rate, or uniform-in-core modulus.

## Exact Cao-to-Buttà--Marchioro scaling

For a Cao solution of core parameter `delta`, circulation `kappa`, and speed
`c_delta=W log(1/delta)`, the exact Euler scaling

    u^{A,B}(x,t)=A u_delta(Bx,ABt)

gives

    omega^{A,B}=AB omega_delta(Bx,ABt),
    xi^{A,B}=AB^2 xi_delta(Br,Bz,ABt),
    kappa^{A,B}=A kappa/B,
    E^{A,B}=A^2 E/B^3,
    P^{A,B}=A P/B^3,
    c^{A,B}=A c_delta.                                  (5)

Consequently `E^{A,B}-c^{A,B}P^{A,B}=A^2(E-c_delta P)/B^3`,
so the rearrangement maximizer and its axial-translation orbit scale with the
functional rather than being assumed to transfer.

Let the unscaled support lie in a `C_s delta` disk centered at radius
`r_delta`, with `r_delta -> kappa/(4*pi*W)`, and set

    B_delta=r_delta/r_0,
    rho_delta=2 C_s delta/B_delta,
    A_delta/B_delta=a/(kappa |log rho_delta|).           (6)

Then the scaled support lies inside the radius-`rho_delta/2` disk centered at
`r_0`, its intensity is exactly `a/|log rho_delta|`, and

    ||bar omega||_infinity
      <= 4 a C_omega C_s^2/
         (kappa rho_delta^2 |log rho_delta|).            (7)

The scaled speed tends to `a/(4*pi*r_0)`. Mutating the `B^-3` energy/impulse
power in (5), or dropping the logarithm in (6)--(7), breaks the direct change
of variables and the source hypotheses; the recorded symbolic verifier checks
the load-bearing functional and intensity scalings exactly.

## Primary theorem applicability, support, and time

Buttà--Marchioro Theorem 1.1 allows `N` definite-sign rings at fixed centers in
the open half-plane. It requires initial support in mutually disjoint
`rho`-disks, intensity `a_j/|log rho|`, and one uniform
`M/(rho^2 |log rho|)` bound. It does not require distinct radii. The scaled
Cao sequence and the small support-preserving same-leaf deformations satisfy
these hypotheses with one `rho`-independent `M`.

For two equal radii `r_0` and equal positive intensities `a`, the source
reference trajectories have the same velocity `a/(4*pi*r_0)` and retain their
initial axial separation `d`. For any fixed

    0<R<min(r_0,d/2),

the theorem supplies `rho_R>0` and `T_R>0` such that, for every constructed
sequence member with `rho_delta<=rho_R`, the complete labeled supports, not
merely most of their vorticity, obey

    supp xi_{j,rho}(t) subset Disk((z_j+a t/(4*pi*r_0),r_0),R)

for `0<=t<=T_R`, and those coarse disks remain disjoint. Hence
`r<=r_0+R` and the full axial support gap is at least `d-2R` throughout that
interval. The theorem separately proves a finer mass-concentration statement;
that weaker statement is not substituted for its support conclusion here.

The quantifiers are `for each admissible fixed R, there exist rho_R and a
positive T_R, independent of rho`. The source does not make `T_R` explicit or
arbitrary, and its output radius is `R`, not a renewed `rho` core. Therefore it
cannot be restarted to infer all-time control. The source trajectories are
asymptotic reference centers; the audit licenses zero limiting relative
velocity, not an exact finite-core centroid ODE or mechanical force law.

For each fixed scaled carrier, combining the full-support gap with (1)--(4)
gives the exact conditional splice

    delta_{0,rho}
      +C_G P_{T,rho}^2/(d-2R)^3-E_{12,rho}(0)
      < m_{rho,r_0+R}(eta)                               (8)

implies

    d_{r_0+R}(xi_{j,rho}(t),bar Sigma_rho)<eta,
    j=1,2, 0<=t<=T_R.                                   (9)

The time is uniform in concentration for the source family, but the Cao
modulus and the allowed initial neighborhood in (8) are carrier-dependent;
no uniform-in-`rho` orbital tube is proved.

## Minimum topology correction

The phrase “relative-open neighborhood” in 0012 needs one concrete topology
repair. The strict support/gap class is not open in Cao's
`L1+L2+P(|.|)` orbital topology on the full rearrangement leaf. A
weighted-volume-preserving rearrangement can move an arbitrarily small piece
of a bounded core outside its tube (or into the separating region) while its
`L1`, `L2`, and weighted-impulse distance tends to zero. Thus the one-parameter
diffeomorphism example proves nonempty nearby deformations, but by itself does
not prove openness in that full topology.

The supported claim becomes exact by naming the relative topology already
used by the construction. For each label, take the space of smooth
`dnu`-preserving diffeomorphisms supported in a prescribed initial tube and
small in `C1`, with the carrier initially having a strict margin to the
`rho`-disk boundary. Give their image pairs the quotient/relative topology
from this support-controlled diffeomorphism space. Energy, impulse, cross
energy, support margin, and deficit are continuous there, so all strict bounds
and (8) define an open set.

“Not generating a translation” should also be replaced by the exact condition
that the chosen image is outside `Sigma`: a nontranslation vector field tangent
to all carrier level sets can leave the vorticity unchanged. Choose one small
weighted-divergence-free deformation with a nonzero normal component to a
regular carrier level set. Its image is not in `Sigma`, so uniqueness gives a
strictly positive deficit. Centering a sufficiently small `C1` neighborhood at
that deformed image preserves positive deficit, support margin, separation,
and (8). This produces a genuine nonempty relative-open set of nontranslation,
nonzero-deficit data without future-history preparation.

This topology specification is the minimum correction. It does not change
the Euler state, any equation, the Cao metric used in conclusion (9), or the
Buttà--Marchioro theorem transfer.

The frozen 0016 README also calls the journal item “Article 9.” The source
receipt and DOI identify Journal of Mathematical Fluid Mechanics 22 (2020),
article 19. The DOI and audited primary body are correct, so this is a
non-scientific bibliography correction only.

## Oracle assessment

The existing four-check 0006 and five-check 0012 receipts are accepted as
recorded and were not rerun. They expose the deficit algebra, exchanged-moment
product bound, support-gap coefficient, optimized Biot--Savart constant,
Euler functional scaling, and intensity scaling. The initial 0006 structural
SymPy false negative is correctly classified as verifier implementation and
its repair compares a simplified difference with zero.

The scripts do not prove source applicability, support openness, or the time
quantifiers. The primary bodies are the stronger oracle for those points. A
sign mutation in (1), a `B`-power or logarithm mutation in (5)--(7), and a
replacement of “positive source interval” by arbitrary-time or of “complete
support” by mere concentration are all detected respectively by the algebraic
identity/change of variables and the literal primary theorem. No numerical or
small-ratio prescription binds.

The verdict YAML safe-loads. Final whitespace validation uses the scoped
tracked `git diff --check HEAD` for the 0011 correction and the equivalent
`git diff --no-index --check /dev/null` checks for each new 0016 artifact; all
produce no whitespace diagnostics. No unchanged scientific oracle was rerun.

## Single joined verdict and strongest theorem

The one joined 0006/0012 route is **established after the minimum topology
correction above**.

Precisely: for every sufficiently small member of the Euler-scaled Cao
polynomial ring sequence, one can place two identical equal-radius,
equal-intensity axial translates far enough apart and choose a nonempty open
set, relative to the named support-controlled `C1` same-leaf diffeomorphism
class, of genuinely deformed nonzero-deficit initial pairs. Each pair is one
classical axisymmetric no-swirl Euler datum. If its actual deficit, total
impulse, initial cross energy, and fixed-carrier Cao modulus satisfy (8), then
on the Buttà--Marchioro interval `[0,T_R]` both complete labeled supports stay
in disjoint radius-`R` tubes and both labels satisfy the orbital estimate (9).
The source interval is positive and independent of concentration; the initial
orbital neighborhood and modulus need not be uniform in concentration.

The remaining dependency for the intended persistent two-carrier mechanism is
a recurrent thin-support estimate, invariant pair tube, or genuine pair-orbit
compactness/stability theorem that closes the exit beyond `T_R`. The result
does not license all-time persistence, a quantitative Hessian gap, generic 3D
or swirl stability, an exact finite-core center law, a pair potential or
mechanical force, collision/exchange dynamics, spin, statistics, charge,
relativity, an electron/neutrino identification, parent completion, or a
global no-go.

## Bounded correction check

The requested single correction pass inspected only the revised topology
passages in `attempts/0012/construction.md`, their synchronized route summary
in `attempts/0012/verdicts.yaml`, and the new
`attempts/0012/review-correction-0016.md` receipt. It did not reopen the exact
deficit, impulse, cross-energy, scaling, support, or time calculations and did
not rerun the unchanged four-check or five-check oracles.

The repair implements the review's minimum correction. The admissible
parameter space is now the group of smooth `dnu`-preserving diffeomorphisms
supported in prescribed initial tubes, with its `C1` topology and strict
support margin. The text chooses an actual deformation whose normal trace on a
regular carrier level set differs from that of an axial translation, so its
image is outside `Sigma`; it no longer infers deformation merely from calling
a generator “nontranslation.” Uniqueness of the maximizing orbit gives a
strictly positive deficit at that image. Continuity then supplies a sufficiently
small open neighborhood preserving positive deficit, support, separation, and
the strict modulus inequality. Passing to the image's group-orbit quotient is
consistent because the quotient map by the carrier stabilizer is open.

The corrected construction explicitly disclaims openness in the full Cao
`L1+L2+P(|.|)` topology, while retaining that metric for the orbital-control
conclusion. Route B2 now uses the same corrected support-controlled topology.
The receipt also records article 19 and correctly leaves the frozen 0016 README
unchanged. The DOI and audited source body remain the same.

Corrected artifacts inspected:

- `attempts/0012/construction.md`:
  `1ca6727f1f797a536e224903036a3246cd7b6603b1fc420207222bb2b2464ec7`;
- `attempts/0012/verdicts.yaml`:
  `0a42dfdecfd6c369c64fd75dd5b41524480bfb66c0ae0ddea6f6a223d034b83b`;
- `attempts/0012/review-correction-0016.md`:
  `3ef474230e67dc67b8ccf046e20dc595e8bd9a08482fa040b42ef8caf142f25c`.

The final joined verdict is **established after bounded correction** at the
exact scope already stated above. The remaining dependency and every exclusion
are unchanged.
