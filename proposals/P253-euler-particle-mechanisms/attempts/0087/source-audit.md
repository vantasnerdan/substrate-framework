# Source and authority audit

## Reviewed inputs

| Input | Frozen hashes | Exact use |
|---|---|---|
| P253/0052 with independent joined P253/0057 review | `0052/derivation.md` `6ea25cd99b91686d0d7efada0bc756656b82b4ebe44f310848511ab63da23e76`; `0052/result.yaml` `1d0781660745266a447bf8ecf2f2e67c7804bfb1572060ad05ce56ec7824cc7a`; `0057/review.md` `c3a01215f13a4ea1bce9ac4b9c86d88f491458c412381a12b86b2e87853d4293`; `0057/verdicts.yaml` `ba850f0fe814a165f5620e5da4f7e1a76109846ec9557981971a92c4de758cd6` | Supplies the physical clock `Omega_a=Gamma/(2*pi*a^2)`, `omega=Omega_a sigma`, and the exact density/Fourier normalization boundary. It does not supply the present radiation theorem. |
| P253/0066 with independent P253/0073 review | Frozen by final P253/0073 | Supplies only the fixed-positive-`k` generalized Sturm--Liouville law `sigma_J=k L_Phi/(pi J)+O(k/J^2)` and positive sectoral energy. It does not supply a joint high-`J`/thin-ring transfer. |
| P253/0074 with independent P253/0078 review | `0074/derivation.md` `ab7695b026e3008a05869a9887cb114c06c6d736e5c6fc0ea95fd8459c50dace`; `0074/result.yaml` `b2ca7ddeacf3ea88e0f1e560ce03c756afa535e66788e4e712c6884c6430a09f`; `0078/review.md` `676fd487d310690cfe89246a5f0a8ba71c676c39c02aa69b4dd2a393a4b2c95e`; `0078/verdicts.yaml` `33725bbee5e149ca152146583614a8b8cc02bcad6ec49ed9e24ea77a06c2fe3f` | Supplies fixed-`J`, fixed-massive-`k` graph/Riesz transfer and positive-Krein modes. It does not supply equal frequency, high-`J` uniformity, or Maxwell radiation control. |
| P253/0079 with independent P253/0082 review | `0079/derivation.md` `ad80c2804798c9144b33bb78d1b9438af0e6f23928edc970825179eb597ad038`; `0079/result.yaml` `8f9747fcba0dfa7698de066d63c8641b501880b9c6274169b2db2e8f5a355cde`; `0082/review.md` `75e55672e359fdf91b91017d6098e266eb605f4ac7e289132153e22a4f7cf334`; `0082/verdicts.yaml` `69eee8688f986e7ee0463d92d9a115f6e397a74f2bbc13370d8d0ba414e175a4` | Supplies the rational-ray crossing atoms, the adjacent fixed-`ell`/high-`J` candidate, and the two-sided gate ledger. The exact same-carrier crossing, physical KKS response, and gate time remain open. |

## Corrected and independently reviewed input

P253/0085 derivation SHA-256
`3bca39e05d66891d5392390f371a84a7c748893f2fe2c5dec915523a9c67a01e`
and result SHA-256
`7fb6c0af9a7530506471c1f69a3fe271c09bd4792ec95c32cc8ad85eb264d299`,
with corrected source-audit SHA-256
`fa0e570fea2554354e50d4340b1884f769cadb7de62785eb0e2dc68189648f37`,
are now consumed through independent P253/0086 review
`4c5c535e60a134ba589ac066aa906249952d47d843c96e7f6a3764686f85ea67`
and verdict
`afd7e563d782d56f8d9f383a8611caf41cf4f2f4bbeb1410dc089262348f9c3d`.
They supply the Doppler shell, coarea row, necessary transverse dark-current
condition, and the full constrained-generator essential inclusion.  No 0085
stability or resonance conclusion is imported.

## Derived in P253/0087

The exact `Gamma=kappa` convention, sharp transverse shell maximum, covariant
three-order vector Bessel ledger, optimized fixed-margin contour bound,
fixed-`J` finite radiation ceiling, and fixed-`ell` finite-`J` leading
predictor at the supplied carrier speed and speed ceiling are derived here.
The helper module checks these
algebraic identities only; it does not supply the joint-transfer remainder.

No value from active P253/0083 is used as authority.  No external empirical
comparator, production numerical result, electron property, neutrino property,
Planck constant, or charge value enters the calculation.
