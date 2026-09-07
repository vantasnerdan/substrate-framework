# P253/0087 author completion receipt

The final replayed README SHA-256 is
`8e4fd5886a6bd7ff28e51786e12639ac875df55caf3530b4980da6bc9dcea7f2`.
The activation chronology and both bounded README precision repairs are
recorded in `postactivation-precision-receipt.md`, SHA-256
`1c1fcbdc12e323c9841bc96e8e793491c1e0c09cad5c965ebc0e20f7987db9c6`.
Every schema replay printed `WORKFLOW VALID`, left stderr empty, and exited
exactly `0`.

## Strongest author-stage result

The exact convention bridge is `Gamma=kappa`, so

    Omega_a=kappa/(2*pi*a^2),
    omega_physical=Omega_a*sigma.

On the comoving Maxwell shell the sharp transverse maximum is
`|omega|/sqrt(c_EM^2-c_g^2)`.  A cylindrical vector `n` character has exactly
the three Bessel orders `n-1,n,n+1`, including after the covariant transverse
projection.  An optimized contour shift gives a rigorous exponential Debye
bound under a fixed argument/order margin.

For fixed-`J` massive modes the physical clock grows as `n^2`; the resulting
argument/order ratio grows linearly and gives only a finite `N_rad` ceiling.
For an actual adjacent fixed-`ell`, high-`J` pair, the finite-`J` leading
predictor is

    q_rad=kappa*ell*L_Phi/
      [2*pi^2*R*k_*sqrt(c_EM^2-c_g^2)].

A buffered `q_0(J)<=1-2*eta_rad` yields the exact memberwise margin only when
the open same-carrier transfer proves `Q_i=q_0+R_i`,
`|R_i|<=C_asym/J`, and the actual integer satisfies
`J>=J_asym>=C_asym/eta_rad`.  There is no fixed-`c_EM` subluminal
`J`-to-infinity limit because the Cao translation speed grows logarithmically.
Physical KKS/current normalization, quantitative response and gate time, and
shell-trace-to-power conversion remain the next achievements.

## Executable evidence

The first focused run exposed only three SymPy representation/assumption
comparisons and is described in `first-focused-failure.md`.  The final
`focused-v3.*` receipt reports seven passing tests, empty stderr, and exit
`0`; stdout SHA-256 is
`42e694012b7e27bb5c20c10b875ab87f9e2ff07c01903380ab24bab29de80308`.
The exact analytic helper reports six passes, empty stderr, and exit `0`;
stdout SHA-256 is
`87237155f39d4cdc5164517f1fc74224199b9a4a53319e8131da2d9dd8174896`.
The repository validator reports `WORKFLOW VALID`, empty stderr, and exit
`0`.

A repository-wide run reported 2784 passing tests, but it was launched before
the 0087 module and tests existed.  It is retained only as a pre-0087
checkpoint result.  It is not claimed as evidence for this attempt.

## Final core hashes before this receipt

- `derivation.md`: `247079dd059c574ee3b96b98305d7504403106e5e5f84ea08acbc2c0d8a17f82`
- `source-audit.md`: `e29f6812f432efc995f0e3c6dd205905b1e5a75722f5994d9bd5b16dc6f2206f`
- `result.yaml`: `356ebe4a5ed0b66dabb26bb6c4294896f502aa08808d586ad2d8abaa0457aabf`
- `validation.md`: `484e5ee5885ad9dd2f5f684aed4bc8c28e602edbcfc2363b6a2d6f6179647780`
- `verify_radiation_barrier.py`: `dd3a31ba54fff032ae8255531b5650ccb0050f03e06a0bfffe0e6de76d12434e`
- `src/substrate_framework/euler_high_harmonic_radiation.py`: `815ffb8334bba2d75b936cdf73f2e39ba5c72e85ecdb9c2a133870c0fc7d263b`
- `tests/test_euler_high_harmonic_radiation.py`: `39ff19b3eb3581941b3a95f449344b4a0ef3ff46af4067d8ca267876bedd67ae`

No production numerics were run.  No particle, quantum, stability, exact dark
mode, resonance-width, or parent-campaign completion claim is made.
