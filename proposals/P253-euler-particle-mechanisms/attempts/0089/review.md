# P253/0089 independent joined review of corrected 0083 and 0087

## Integrity and correction check

This review executed against the centrally activated contract at README
SHA-256
`500c428c5d5e119b5656ebfcef5a5b3ac3956f195cee2cc8b9d898782778ff69`.
`particle-balance-review` authored or implemented neither target. P253/0088
remained unopened and excluded. No production numerical run or unchanged
oracle rerun was performed.

The single bounded correction passed. Its receipt
`0083/0089-bounded-correction-receipt.md` has SHA-256
`8894fdca60b3b39daad8c0759720aed27145923d5beb0583648a89892aefe256`.
The corrected 0083 README, derivation, result, source audit, and validation
hashes are respectively
`9030a1dc7763414df79968854ac0fc678004ccaba0e253b4070f543931388a3c`,
`0dc3a837a201a97fd0bd89d9316992884e4b82f4d4fed8027d30f8473f5a8003`,
`635216a567b2bdf92a12eba4a97cd46244f79fc800a55692a8f9204d45b07534`,
`326aff97dbe136b20aef4f05b15efa8f65088f81066ce9f1303c147e02bc2a03`,
and
`bb727a5c6cc53303669767e55213c5653992216066a4a15f0ef833dafc09b216`.
The corrected 0087 README, source audit, validation, verifier, module, and test
hashes are respectively
`2582259baff42b64333affb85c12cb2a217d26717f50796ba76b14936824b77f`,
`5cdaf18dc60c880d353fd2d7444ccc8a697b9768a01a3df5b6e0ef6246042e4a`,
`bcfcb009b44d649a598c3fe2bc54006e9323f8b811f1cbeb9b07c0bf49a7d6a4`,
`740be5d6480738e051778374b68b0790a61c3dc0cf2c0f97badf470a2f11c3a9`,
`8b3c664e19fd64f1b2067fd83b9613b2a27c0b77d2dab8be5cd3b72e5dff5aec`,
and
`c05480cd5432b08f33f6babe3f5e934eb4a1d4e77b6a55b4ac186243013bcf51`.

## Source and supplier boundary

Cao et al., arXiv:1910.07493, supplies the thin fixed-parameter vortex-ring
source construction and its local asymptotics. It does not itself state the
target's fixed-mean-radius connected external-epsilon path, uniform
same-family spectral envelope, or exact two-frequency crossing. Those are
direct constructions using the accepted 0080/0084 `B_R` inverse and the
independently reviewed 0066/0073 and 0074/0078 graph/Riesz machinery. The
0084 result is consumed only at its fixed-circulation, exact-mean-radius
scope; no sharp Schur theorem is imported.

The Bessel recurrences and Debye inequalities support 0087's special-function
ledger. They do not provide a physical current, shell trace, KKS
normalization, or radiation rate. Corrected 0085 is consumed only through the
independent 0086 review at its final hashes.

## Unit A — external-epsilon Cao path

The scale coefficient is consistent with the local core equation and
circulation row:

    C_s=(Lambda_p/kappa)^((p-1)/2) R^(-(p+1)/2).

Writing `r=x_epsilon,1+xi`, the leading radial cell gives a second moment of
order `kappa epsilon^2`; the relative-`O(epsilon)` odd correction contributes
at the same order after multiplication by the center; and the remaining row
is `O(kappa epsilon^3 |log epsilon|)`. Thus the exact mean constraint gives
`R-x_epsilon,1=O(epsilon^2)`. Lemma A.2 only supplies the required relative
boundary control, so the supported value statement is

    s_area=C_s epsilon+O(epsilon^2),
    delta=(C_s/R)epsilon+O(epsilon^2),

without differentiating either remainder.

The correction supplies the missing global-in-the-thin-interval bridge. The
parameter-only dilation and centering identify all source maps with one fixed
compact `C^(2,alpha)` space. One uniform inverse, residual modulus,
contraction radius, positive-interface margin, and `C^(2,beta)` compactness
bound give overlapping multiplicative relative-epsilon charts. Local
uniqueness identifies overlaps and prevents positive-epsilon termination by
norm, support, or domain drift inside the declared thin interval.

The crossing envelope is now noncircular. `bar rho_N` is defined first over
the preliminary `1/L` bands and the actual integer fibers `LP,LQ`. Sequential
profile compactness, normalized-column uniqueness, common-domain graph
convergence, the fixed-contour resolvent identity, and rank-one simplicity
give `bar rho_N -> 0`. Only then is
`h_N=sqrt(max(bar rho_N,1/N))` defined. The threshold

    N>4 P A_cov/k_*

makes the lower endpoint positive and places both endpoints inside the same
preliminary band and connected thin path. The scale and spectral errors are
`o(h_N)`, so the accepted 0079/0082 IVT argument applies uniformly along the
actual carrier segment.

**Unit A verdict: established as stated after the bounded correction.** The
fixed-`(kappa,R)` path and its exact same-family crossing supplier are proved
on the declared thin interval. This is local continuation plus finite
buffered coverage, not global branch uniqueness.

## Unit B — response reduction

For a compact dynamically accessible displacement `h=[xi,Omega]`, integration
by parts and the full Hodge velocity give the exact off-diagonal response in
the form

    integral xi dot F_12.

Finite physical rows are correctly handled by a bounded primal oblique map
`P_row`. The dual witness is `P_row^* F_12`, since

    <P_row xi_0,F_12>=<xi_0,P_row^* F_12>.

No self-adjoint projection is assumed. Hence nonzero
`curl(P_row^* F_12)` is a sufficient local exposing predicate. The target
does not evaluate that predicate, normalize the physical KKS covector, or
derive a gate-time coefficient.

**Unit B verdict: the constrained dual response identity is established.**
Physical nonvanishing and KKS/action normalization remain open constructions,
not defects in the identity.

## Unit C — circulation clock, Doppler maximum, and Bessel ledger

The cross-sectional circulation convention gives `Gamma=kappa` with no extra
azimuthal `2 pi`; the associated core clock is
`Omega_a=kappa/(2 pi a^2)`. Maximizing the full subluminal Doppler shell gives

    k_perp,max=|omega|/sqrt(c_EM^2-c^2),

with the stated optimizing direction. Converting cylindrical components
before transverse projection produces the exact `n-1,n,n+1` vector-Bessel
orders. The Debye exponent
`acosh(1/q)-sqrt(1-q^2)` is positive precisely for `0<q<1`.

**Unit C verdict: established as stated.** These are physical-clock and exact
shell/special-function identities, not a current or flux theorem.

## Unit D — fixed radial-index route

Holding the radial/Sturm label `J` fixed while `n=Np` grows makes the relevant
Bessel ratio grow linearly, so the protection inequality has only a finite
upper ceiling. It does not produce a nonempty protected integer pair. Here
`J` is a mode index, not physical action.

**Unit D verdict: established finite-ceiling obstruction.** It neither
quantizes action nor rules out other finite or scale-breaking routes.

## Unit E — fixed harmonic-step/high-index route

For fixed nonzero integer toroidal step `ell`, the pair
`n_1=J ell`, `n_2=(J+1)ell` has the finite-`J` leading predictor

    q_0(J)=kappa ell L_phi /
      (2 pi^2 R k_* sqrt(c_EM^2-c_g(J)^2)).

Agreement with the exact finite carrier requires an actual
`J>=J_asym` and remainders `|R_i|<=C_asym/J`, all below the independent
`H_loss`, Cao-speed, `J_rad`, and electromagnetic subluminality ceilings.
Because `c_g(J)` grows logarithmically, there is no fixed-`c_EM` subluminal
`J -> infinity` Cao family. The public helper keeps its historical API name
but now truthfully documents only a finite-`J` leading predictor at a supplied
carrier speed; verifier and test labels agree.

**Unit E verdict: the leading finite-window algebra and conditional agreement
criterion are established.** No nonempty protected pair or asymptotic family
has been constructed.

## Unit F — flux and gate boundary

A nonzero transverse Fourier-current trace on the real Maxwell shell is an
outgoing-radiation obstruction to an `L2` time-harmonic field. It is not
itself an `L2` mode, a flux magnitude, or a decay rate. Turning it into a
leakage coefficient requires the actual normalized current, a
limiting-absorption/time-domain flux construction, physical KKS mode action,
and a gate-time comparison.

**Unit F verdict: established as a conditional boundary only.** No radiation,
gate, particle, P4/P5, electron, or neutrino conclusion follows.

## Strongest supported result and next construction

Corrected 0083 supplies the missing connected fixed-`(kappa,R)` thin-carrier
path and uniform actual-fiber envelope, so the independently reviewed column
crossing transfers to an exact same-family Cao crossing. Corrected 0087
establishes the physical circulation clock, sharp Doppler and vector-Bessel
ledger, a fixed-index finite ceiling, and a finite-`J` conditional agreement
window. The response functional is exactly reduced, but its constrained curl,
physical KKS normalization, actual current trace, and leakage/gate conversion
remain unevaluated.

The next positive construction is to compute
`curl(P_row^* F_12)` for the exact crossing modes, normalize their physical
KKS action, and evaluate the resulting transverse current on the Maxwell
shell. Only then can a limiting-absorption or time-domain calculation decide
whether a nonzero leakage coefficient and usable gate window exist.
