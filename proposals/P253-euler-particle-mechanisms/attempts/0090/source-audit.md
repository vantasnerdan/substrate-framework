# Source and authority audit

## Reviewed framework inputs

| Input | Frozen hashes | Exact use |
|---|---|---|
| P253/0068 with independent P253/0071 review | `0068/derivation.md` `41bfff94e0c72d0c028564b2ee7c04d231fc3ea3eca75d1d7c556facb86d7816`; `0068/result.yaml` `ec029045305fee0bfabf53520fad2f4820b7ec775367d630e64848c0bdfe0d5e`; `0071/review.md` `50e1cb63d6c115a92c3fe6720eaeae0b707c5b81871db663a7a804a9a3127155`; `0071/verdicts.yaml` `fd59cc3209460f85aa97cacf66ea5a88d6d5cd9ce8a1630b8e37cf1f31d0f1c0` | Supplies the exact conserved U(1) current, Maxwell action, Poynting work balance, gauge/constraint boundary, and Maxwell-only causal propagation. It supplies no harmonic radiation rate or gate theorem. |
| Corrected P253/0085 with independent P253/0086 review | `0085/derivation.md` `3bca39e05d66891d5392390f371a84a7c748893f2fe2c5dec915523a9c67a01e`; `0085/result.yaml` `7fb6c0af9a7530506471c1f69a3fe271c09bd4792ec95c32cc8ad85eb264d299`; `0086/review.md` `4c5c535e60a134ba589ac066aa906249952d47d843c96e7f6a3764686f85ea67`; `0086/verdicts.yaml` `afd7e563d782d56f8d9f383a8611caf41cf4f2f4bbeb1410dc089262348f9c3d` | Supplies the comoving Doppler shell, its smooth radial-root/coarea structure, the transverse dark-current necessity, and the full joint `i R` essential-spectrum boundary. It supplies no outgoing normalization, finite-time flux, resonance, or mode current. |

## Derived here

Poynting's theorem, the transverse Maxwell wave equation, the retarded
distribution identity, coarea, the Fourier trace theorem, and the elementary
forced-oscillator representation are applied directly under the convention
frozen in the README.  No empirical comparator or fitted constant is used.
The sign and constants are rederived in `derivation.md`; the importable helper
and verifier encode their algebraic consequences rather than serving as a
source.

The corrected plateau proof is internal as well.  It pushes the exact
resonant rectangle energy to the scalar mismatch variable, uses the exact
rectangle autocorrelation, and bounds the remainder by the first time moment
of that correlation.  Fixed endpoint profiles include ramp corrections and
the complete conserved `a_prime*K` term.  Their cross correlations, endpoint
energy, and the uniformly off-shell conjugate temporal branch are separate
hypotheses.  Completing the anisotropic comoving symbol to a shifted and
rescaled standard Helmholtz symbol types the weighted LAP bridge.

P253/0087 and its later independent P253/0089 review remain outside this
attempt's input boundary.  The reviewed Bessel shell estimate can consume the
present power functional later without being used to prove it.
P253/0088 remains active and is excluded; no response, KKS normalization, or
gate history is imported from it.

## Exact boundary

The prescribed-current outgoing and late-field identities belong to the
declared Euler--Maxwell extension with explicit `epsilon_EM`, `mu_EM`, and
`c_EM`.  They neither derive those constants from Euler nor select a charge or
action.  The reciprocal amplitude conversion remains conditional on an actual
same-action current map and analytic one-block closure.  No electron,
neutrino, Born/reset, exchange/statistics, P4, P5, or parent result follows.
The moving-control-volume identity is exact at finite volume; equality of an
infinite-world-tube limit with late free-field energy is not imported.
