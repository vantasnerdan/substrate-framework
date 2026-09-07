# P253/0085 author completion receipt

The corrected activated README remains SHA-256
`874522770f9272d73e538435d016ce7a6f1dae7a60217c3ff77f33f7acc07847`.
The postactivation spectral-scope receipt is SHA-256
`f7af4be104e354534c01cc62f12d335d5c54f1fa9233fd4e0843f45281fcc34c`.
The exact activation replay exited `0` with `WORKFLOW VALID` and empty stderr.

## Strongest author-stage result

The exact subluminal field relative energy completes to a positive quadratic.
Minimization on the instantaneous Gauss/`div B` tangent gives the positive
anisotropic charge form

    (1-beta^2)/(2 epsilon)
      integral |qhat|^2/[|k_perp|^2+(1-beta^2)k_z^2],

which is equivalent to homogeneous `H^-1` for a fixed subluminal margin and
is `O(g^2)` on `delta chi`.

The full joint relative Hessian nevertheless has infinite positive and
negative index.  Pure compact Maxwell radiation gives an infinite positive
subspace after the finite total-momentum/gauge rows are removed.  On a common
untagged regular Cao cell, the dynamically accessible fluid form is a
strictly negative multiplication operator plus compact Green operator; a
finite-codimension subspace is uniformly negative.  Those two strict signs
persist for sufficiently small charge using the independently accepted
P253/0080/0084 charged branch.  This refutes ordinary one-sign joint
coercivity and finite Morse index, not spectral stability.

The full translating-frame Maxwell spectrum fills `i R`.  Scaled compact-curl
cutoff plane waves, including a `k_L->0` sequence at zero frequency, give
normalized constrained weak-null graph sequences and prove that inclusion in
the joint Fredholm essential spectrum at a localized charged carrier.  An
imaginary Cao mode is therefore embedded.  Its exact necessary
dark-current condition is vanishing, to sufficient order, of the transverse
trace of `g(chi delta_u+u delta_chi)` on every component of
`|k|^2=(omega+c k_z)^2/c_EM^2`.  The derived star-shaped shell and coarea
weight give a nonnegative functional for the next source-specific
limiting-absorption/Feshbach calculation.  A nonzero trace excludes an `L2`
Maxwell field at the prescribed real frequency but does not itself construct
a resonance or determine its width.

The exact regular-band tag identities show that the chosen `chi=F(I)` removes
same-orbit tag-only stabilizer directions where `zeta zeta'` is nonzero.  The
electric energy still controls only a `g^2 H^-1` norm and degenerates as
`g->0`.

## Executable provenance

The early test failures were representation-only comparisons between factored
and expanded SymPy expressions and are preserved in `first-focused.*`,
`focused-tests.*`, and `shell-first.*`, with the repair documented in
`first-focused-failure.md`.  Final `focused-v3.*` reports seven tests passed,
exit `0`; final `exact-v2.*` reports six independently derived assertions,
exit `0`; both stderr files are empty.  No production numerics were run.

## Final hashes before this receipt

- `derivation.md`: `da2bd59eb3e26a17fa81955d889bc3e0b2b1f0137e8c6a00b3061e7c173a5082`
- `source-audit.md`: `82f7f130154a7048ef4fb986ff10d4737eb9572e7945056e9e754f15a12c29b5`
- `result.yaml`: `e86304d7682b4fb532b8c3324077f1c7259e1b5a9f7e9841d54f4330a18c295b`
- `validation.md`: `21721bcd9c1dcfb2200840d56bec615def159dbce1280ad5be0eeaf8d4b89579`
- `verify_charged_hessian.py`: `b055628193247da57958d2e3d5f12beea8cf8443fe920e96eb12d9a47e3f01ea`
- `src/substrate_framework/euler_charged_hessian.py`:
  `c17fc51ec7f98051b082432b9dc9f2d807685c1d91c704ae799f45eef037ac30`
- `tests/test_euler_charged_hessian.py`:
  `5acb87a2a9140a8b42830a2d3ec625e74b86a9f2a4ecb5e77b9d95b0de9a0fd6`
- final focused stdout:
  `0e0177d95805cd7e9ccb1859b562b7fe1d8ffaed73559d85239e8afb6ed7a331`
- final exact stdout:
  `da4adb00fe224bd2004e8a18f34013d1e92402954c82567f49dfe270bafb81cb`

## Continuing achievements

Evaluate the on-shell current for the two actual P253/0079 modes with physical
KKS normalization, then construct the one-block weighted Maxwell
limiting-absorption/Feshbach map.  In parallel, evaluate the finite low block
of the steady-slaved axisymmetric Hessian and formulate a modulated theorem
that retains outgoing radiation.  The bare-Euler carrier and universal-action
selection routes remain active, as do every neutrino-specific obligation.
