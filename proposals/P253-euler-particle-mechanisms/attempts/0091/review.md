# P253/0091 independent review of corrected 0090

## Integrity and correction closure

This review executed against the centrally activated README at SHA-256
`9dbd1057bd928c75ec3252ec880be860f15a6af448f30f3057853403711998ec`.
`particle-balance-review` authored or implemented none of 0090. P253/0088
and all P253/0089 evidence remained unopened and excluded. No production
numerics or unchanged oracle/API rerun was performed.

The single bounded correction passed. Its receipt
`0090/0091-bounded-correction-receipt.md` has SHA-256
`0a82e8be80ed32caa2d7b8338e0f1740fd2807dbd4f4c288c945e4c13b5a632f`.
The corrected derivation, source audit, result, and validation hashes are
respectively
`c2c39e2ef0585647dfe3d3238c9c2999cf6bbe970389df8528706697985b0360`,
`43140e6eecfe2411fd7b31f550eaba6416ac040d6babdbebe1077c74305a789e`,
`4a5ff2e67ea6d3713dede6fd269eb4402ac8d5e9c1973c338376546578f87b3c`,
and
`a5cd1523bc6c653d5902577a68159de8196645125f173517de1927323436a2ea`.
The same-transaction ledger repair passed at completion-receipt SHA-256
`ccb6a22ed484cc1389b4ad7e768d6f69a59ff4a400ebf66dd8ea2103575c17e6`
and manifest SHA-256
`71f49aa3c3b8d4f1fd3b28a975c397be54732b5d95c392649561084039da9678`.
The historical quartet is labeled pre-review/pre-correction and the final
quartet and manifest agree.

## Source and analytic boundary

The reviewed 0068/0071 Maxwell action supplies Poynting balance, gauge and
constraint conventions, and Maxwell-only causal propagation. Corrected
0085/0086 supplies the comoving shell and dark-current necessity. Neither
supplier proves the outgoing normalization, the finite-window limit, an
actual carrier current, or a reciprocal resonance. Those are correctly kept
as direct 0090 calculations or open dependencies.

The comoving symbol completes exactly as

    D=c_EM^2 |k_perp|^2
      +(c_EM^2-c_g^2)(k_z-omega c_g/(c_EM^2-c_g^2))^2
      -omega^2 c_EM^2/(c_EM^2-c_g^2).

Modulation and a positive anisotropic dilation reduce it to a standard
nonthreshold Helmholtz symbol. On a fixed subluminal margin they preserve the
weighted spaces, so the stated `L2_1 -> L2_-1` limiting-absorption bridge is
applicable. The exact algebraic oracle checks constants and signs only; it is
not counted as proof of this analytic map or of the time limit.

## Unit A — neutral source and outgoing normalization

For the convention `exp(-i omega t)` and
`Omega=omega+c_g k_z`, retarded evolution gives

    A_T=mu_EM c_EM^2 (D-i0 sign Omega)^(-1) P_T J,
    E_T=i Omega A_T.

Combining complex conjugation with
`(D+i0 s)^(-1)=PV(1/D)-i pi s delta(D)` makes the dissipative part of
`J dot conjugate(E)` negative. Poynting's theorem therefore gives

    P_out=pi/(2 epsilon_EM)
      integral |Omega| delta(D) |P_T J_hat|^2 dk >=0.

The factor `1/2` is precisely the real-phasor time average; adding a second
identical shell would double count. Net neutrality is earned either as a
literal zero-mean hypothesis or from integrable charge and profile flux, for
which continuity is evaluable at `k=0`. The separate weighted-current trace
hypothesis alone is not used for that point value.

**Unit A verdict: established as stated on the declared integrable/weighted,
nonzero-frequency, strictly subluminal prescribed-current class.**

## Unit B — shell geometry, trace, and Gaussian source

The shell has the unique positive radial root

    r(n)=|omega|/(c_EM-sign(omega)c_g n_z),

and

    grad D=2(c_EM^2 k-c_g Omega e_z),
    |partial_r D|=2 c_EM |omega|.

Its full gradient is bounded away from zero for every fixed subluminal margin.
Spherical-coordinate coarea gives the exact weight

    pi |omega|^2 /
      [4 epsilon_EM (c_EM-sign(omega)c_g n_z)^3].

The condition `<x>J in L2` gives `J_hat in H1`; the compact shell is away
from `k=0`, so its `L2` trace and the transverse projector are well defined.
For `J_hat=i j0(k cross a)exp(-sigma^2|k|^2/2)`, direct angular integration
gives

    P_G=2 pi^2 j0^2 |a|^2 |omega|^4
      exp(-sigma^2 omega^2/c_EM^2)/(3 epsilon_EM c_EM^5).

**Unit B verdict: established as stated.** The public helper and exact oracle
are faithful algebraic support for these identities, not evidence for a
general LAP or physical carrier current.

## Unit C — conserved switching and late free field

For zero-mean `rho`, a declared Bogovskii/Hodge completion
`div K=-rho` makes

    rho_T=a exp(-i omega t) rho,
    J_T=exp(-i omega t)(a J+a' K)

exactly conserved, including translation of the profile. Its transverse
Fourier numerator is

    A_T(xi) P_T(J-i xi K),
    xi=c_EM|k|-Omega.

For a smooth compact-time envelope and zero incoming radiation, the full real
current gives the exact late transverse energy

    E_late=(2 epsilon_EM)^(-1) integral |F_T(k)|^2 dk.

Both conjugate phasors and their finite-window interference occur inside this
single `F_T`; no one-sided phasor substitution is made. The zero-mean charge
condition and the weighted `P_T K` energy condition are correctly distinct.

**Unit C verdict: established as stated for smooth conserved envelopes in the
declared source spaces.** This is a prescribed-current Maxwell theorem, not a
carrier-specific radiation result.

## Unit D — plateau limit, sharp endpoints, and moving volume

The correction makes the plateau theorem explicit. With

    H(xi)=(8 epsilon_EM)^(-1)
      integral delta(xi-[c_EM|k|-Omega]) |P_T J_hat|^2 dk,

one has `P_out=2 pi H(0)`. The exact rectangle autocorrelation identity yields

    |E_rect(T)-T P_out| <= 2 integral |tau| |C(tau)| d tau.

Fixed ramp profiles include the entire `a'K` endpoint pair. Their finite
energies and integrable rectangle-endpoint correlations uniformly control
the `J-K`, `K-K`, and ramp cross terms. The conjugate branch is separated by
at least `|omega|`; its self term and its cross correlation are uniformly
bounded under the stated weighted and first-moment hypotheses. Consequently

    E_late(T)=T P_out+O(1)

is established with an explicit `T`-independent norm bound under exactly
those hypotheses. Without them, only the exact spectral norm survives.

A characteristic switch remains a distributional weak limit whose endpoint
energy must be computed for its chosen smoothing. Reynolds transport proves
the exact finite moving-volume identity with relative flux
`(S-c_g u_EM e_z) dot n`. Equality of an infinite-world-tube limit with the
late free-field energy requires an additional decay and endpoint-limit
theorem and is not claimed.

**Unit D verdict: the smooth fixed-ramp `T P_out+O(1)` theorem and finite
moving-volume balance are established after correction.** Sharp-switch
endpoint evaluation and infinite-world-tube equivalence remain open.

## Unit E — reciprocal factor of two

If a mode with action `A_mode |z|^2`, energy
`nu_phys A_mode |z|^2`, and current `z j+conjugate(zj)` belongs to one
reciprocal joint action and obeys
`dot z=(-i nu_phys-gamma)z+...`, energy balance gives

    gamma=P_out[j]/(2 nu_phys A_mode),

while the fractional energy/action decay rate is twice `gamma`. This is a
dimensionally correct conditional identity. The prescribed-current problem
does not construct the amplitude equation, pole, real shift, Markov closure,
or physical action normalization.

**Unit E verdict: established as an exact conditional conversion, not as a
damping theorem.**

## Unit F — interpretation boundary

0090 supplies no actual Cao/analyzer current, nonzero physical shell trace,
KKS action, mixer response, gate history, or reciprocal outgoing Feshbach
map. Switching is externally prescribed. The result therefore supplies no
autonomous reset, Born rule, statistics/exchange character, universal action,
particle, P4/P5, electron, neutrino, stability, or bare-Euler finite-speed
conclusion.

**Unit F verdict: established boundary.** The strongest positive statement is
the exact prescribed-current radiation and smooth finite-window theorem.

## Strongest result and next construction

Corrected 0090 establishes an exactly normalized nonnegative outgoing power
for a localized conserved subluminal comoving Maxwell current, its exact
shell/coarea representation and Gaussian check, an exactly conserved smooth
switch, the full-real late-field spectral energy, and—under explicit
correlation and endpoint norms—the uniform law
`E_late(T)=T P_out+O(1)`. The factor-two amplitude/energy conversion is exact
only conditional on a reciprocal normalized mode equation.

The next construction is to supply an actual constrained crossing-mode
current with physical KKS normalization, evaluate its transverse shell trace,
and build the reciprocal limiting-absorption/Feshbach closure. Sharp-switch
endpoint energy and the infinite moving-tube limit remain separate analytic
routes rather than prerequisites for the established smooth-source theorem.
