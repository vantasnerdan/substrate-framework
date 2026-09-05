# Independent review: actual quadrature phase control

Reviewer: `/root/smooth_core_review`, 2026-09-05. Author of 0210: `/root`.
The registered 0214 boundary was activated after central schema confirmation
at 269 accepted claims. Base authority: v0.181.0. Write surface: 0214 only.

## Result and independence

**Established as stated at the declared fixed-cell, fixed-time control
scope. No scientific correction is requested.** The construction adds an
arbitrary finite signed phase using actual positive-density Euler/Lin
preparations, with zero complete initial energy at K=0. Positive whole-field
inversion cancels its genuinely nonzero first spatial energy term. Its
second energy coefficient and remote physical observation errors tend to
zero with the stated, compatible high-harmonic/macroscopic-scale ordering.

I did not author or implement 0210 or its 0205 passive/off-flow source.
I authored earlier optical evidence including 0155, which 0205 separately
uses for its optical supplier. That supplier is an unchanged, separately
scoped input here: this review does not count my reading of that part of
0205 as independent validation of 0155. The new control proof uses the
actual X-independent background, passive sector and off-flow propagation;
it does not need my earlier optical coefficients to establish its phase
and energy assertions. No new clock, EPS or autonomous-parent condition
has been added to this review.

## Independent phase and constrained-energy calculation

The actual coarea measure is invariant under T=omega(c) partial_theta.
For both quadratures the displayed h and g satisfy Th=g/2 pointwise.
The exact histories have eta=exp(-tT)(h+t g) and canonical cotangent
pi=rho exp(-tT)g. Consequently their phase is

    Omega_ij=rho <h_i g_j-g_i h_j>,

because the two terms proportional to t cancel. Direct substitution,
before angular integration, gives

    h_1 g_2-g_1 h_2=G(c)^2/(2 N omega(c)).

Thus B_N is strictly positive, with the factor in (1), not half or twice
that factor. The complete energy matrix is zero because
g_i g_j-g_i Th_j-Th_i g_j=0. Conservation at later times also follows
directly: its remaining term is proportional to the integral of
T(g_i g_j), which vanishes in the invariant measure. Column interchange
reverses phase and preserves the zero energy matrix. Neither operation
changes the sign of mass density or ensemble probability.

For the finite-K claim I independently eliminated the constrained rate.
Write P=P_K, A=A_K, eta=P(h e_X), and projected cotangent P(g e_X).
The solenoidal initial rate is

    eta_t=P(g e_X-A eta).

The divergence identity for two solenoidal fields gives
(I-P)A eta=(I-P)Du eta. Therefore its actual Euler velocity is precisely
P(g e_X-Du eta), as used in the proof. Expanding
rho[||eta_t||^2-||A eta||^2+<eta,Hess(p0)eta>]/2 gives (2), including
the NEGATIVE longitudinal-pressure norm. The longitudinal constraint has
not been dropped or replaced by a favorable quadratic form.

At K=0 both h and g point along e_X, while the first projection derivative
is transverse. Since Hess(p0)e_X=0, all first-order projection and pressure
cross terms in the energy vanish by this vector structure. The surviving
Hermitian matrix is exactly

    H_ij^(1)=i rho <(K.u)(h_i g_j-g_i h_j)>.

It is generally nonzero even though its real diagonal entries vanish.
Taking only a real diagonal would miss it. Under the whole-field pair
R and -R the axial input n=det(R)R e_X is the same, but this polar
transport coefficient changes sign. Positive averaging cancels H^(1),
while both phase contributions retain their sign. This is a whole-law
statement, not cancellation of an individual realization's current.

## Pressure orders, physical observations and the scale window

The nonzero transverse Fourier multiplier has P_XX=1-K_X^2/|n|^2+O(K^3).
The absence of a first phase correction and the order-minus-two second
correction follow from this actual projection, not from omitting pressure.
After phase normalization,

    ||g||=O(N^(1/2)),       ||h||=O(N^(-1/2)),
    ||P' h||=O(N^(-3/2)),  ||P'' h||=O(N^(-5/2)).

Transport raises the carrier order by one. In particular the quadratic
pressure term involving A P' h has size O(N^-1), not order one.
The other second-order terms g A P'' h, P'g A P'h, and the explicit
iK.u h terms are also O(N^-1) or smaller. The same counting gives
O(N^-2) for the second phase correction. Smooth nonstationary angle
integration controls the low Fourier tails to every fixed negative
Sobolev order; the fixed cell has no shrinking nonzero Fourier gap.

The mean-harmonic convention is the one explicitly supplied by
0196/fixed-cell-diagonal-limit.md: differentiate K=k kappa with the
physical mean projector P_kappa fixed, uniformly in unit kappa. The
nonzero-mode derivatives obey its multiplier bounds. No smooth extension
of K K^T/|K|^2 at the origin is needed. The initial projected configuration
and projected cotangent have zero mean. The physical velocity completion
can have a small first-K mean from -P_kappa<Du eta_K>; this is retained
in the proof's complete physical-current estimate, not asserted zero.
These clarifications preserve the actual source boundary and all the
claimed diagonal estimates.

The 0205 off-flow argument applies to the full Euler/Lin system here:
its scalar principal transport has the actual flow graph, with order-zero
matrix/pressure corrections. The fixed source bands and transported tag
are separated from that graph by their invariant psi gap. The restricted
propagator and its finite K/time derivatives are therefore smoothing.
Integration by parts in theta gives arbitrarily high negative powers of
N, including after the square-root-N normalization. Complete mean-current
rows can instead be estimated by propagating their smooth adjoint test
functions; the finite-rank mean projector is spatially smoothing as well.
Thus this assertion is not a local-pressure approximation or an
identification of a canonical row with tagged mechanical spin.

All profile, field, time and nonzero tag-denominator constants are fixed
first. Sobolev growth of the preparations is polynomial, and the actual
linear Euler/Lin evolution bound depends on the fixed background, not
on their amplitudes. Enlarge the finite exponent L to cover the complete
finite list of cubic action and observation remainder constants. Then

    K_N=c_0 N^(-L-1),
    N^L K_N=O(N^-1),
    N^-q/K_N=O(N^(L+1-q))

prove the two required normalized errors tend to zero for q>L+2.
The zero-order tagged error is exactly zero. This explicitly resolves
the competing small-K and first-gradient-error requirements. It is a
fixed-time diagonal construction, not a uniform acoustic-time estimate
or a claim about one preparation with bounded norm as N tends to infinity.

## Finite constraints and physical meaning

The full initial cross phase/energy coefficients against a previously
fixed finite collection are linear in the new G. Complex rows can be
split into their real and imaginary parts before counting constraints.
For M real rows, M+1 disjoint fixed band profiles supply a nonzero
homogeneous kernel vector without a generic-rank assumption. Its unit
coefficient normalization gives a strictly positive uniform lower bound
on integral mu G^2/omega by disjoint support. N-dependent kernel vectors
therefore do not defeat B_N~N^-1 or the observation estimates. Subsequent
controls include earlier cross rows among the finite constraints.

The raw phase b contributes b/3 to the isotropic vector phase because
E[n n^T]=I/3. The reconstruction Phi=3E[n theta] is an observation map;
it does not multiply the complete-fluid action by three. Thus raw
b=3 delta_j, with any family probability included, has the stated meaning.

Zero added energy and finite added phase are compatible because these
are actual solution columns with a time-dependent embedding, not an
invariant oscillator plane assigned a new generator. The proof explicitly
retains the embedding connection and conserved energy separately. Its
application to a separately supplied literal spin coefficient is valid
at that interface; this control does not create or alter that supplier's
physical clock at K=0. That distinction is part of the positive theorem,
not an additional debt item.

## Evidence and decision

The strongest oracle is the analytic constrained-energy, projection and
off-flow estimate above. I read the complete 0210 proof, verifier, first
output and receipt, the complete relevant 0205 proof bodies and receipt,
the passive_packet implementation and the direct 0196 Bloch estimate.
The captured 16/16 first-run receipt is reused: its predicates test the
actual API-generated quadratures, both Lin columns, zero energy matrix,
phase sign/factor, full projection orders, nonzero Hermitian first energy,
inversion and compatible scaling. They are exact anchors, not a substitute
for the analytic operator estimates. No numerical floor or new numerical
design is involved, and no unchanged full validation was rerun.

- Verification: analytic derivation supported; 16 exact anchors reused.
- Review: one independent substantive pass; no load-bearing finding.
- Compatibility: actual C016 cell and 0205 hypotheses retained.
- Epistemic role: independently reviewed constructive evidence attachment.
- Parent relationship: useful phase/energy interface, not campaign completion.
- Correction check: not needed.

## Content-addressed source boundary

The actual file hashes match the author's captured receipt:

| Source | SHA256 |
| --- | --- |
| [0210 proof](../0210/quadrature-phase-control.md) | `452c68fc4c86bb4111d8087adc0017ec35a664a9b5eaeaa976bd4b1d1bb2578a` |
| [0210 verifier](../0210/verify.py) | `bfc5aa2df1441f326e2d4dcb23e2e1a5d60a40ab63929a439794c397061621c0` |
| [0210 first output](../0210/first.stdout) | `f6839eece3a38ee2d609662548417500e62ff3a9b67c369c605ffa5dbbb02019` |
| [0205 passive/off-flow proof](../0205/physical-energy-returns.md) | `1c066076af365e9817d57e00e4ecdd36107c47a15d5d445fff756b6c6575bc6b` |
| `src/substrate_framework/euler_passive_control.py` | `79c7191821a7aa9ea04f201dfb103d086b82b1610c28e5a512789b6a9e2e1dd4` |

Signed: `/root/smooth_core_review`.
