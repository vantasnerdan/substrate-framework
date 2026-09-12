# Validation

## Strongest exact oracle

The bounded exact oracle is

    /home/dan/substrate-framework/.venv/bin/python \
      proposals/P253-euler-particle-mechanisms/attempts/0079/verify_schwinger_hopf_blocks.py

The corrected run is preserved in `symbolic-v4.command.txt`,
`symbolic-v4.stdout`, `symbolic-v4.stderr` and `symbolic-v4.exit`; exit is
exactly `0`.  It verifies by exact SymPy identities:

1. `J^2=-1`, `[V_C,J]=0`, `{V_A,J}=0`, and `V_C+V_A=V` for a generic real
   four-by-four matrix;
2. realification of a complex-linear Hamiltonian generator commutes with `J`;
3. the diagonal/off-diagonal Pauli axes are non-collinear, while scalar
   multiples are rejected;
4. the corrected full-system gate ledger depends separately on outbound
   leakage and the return Volterra norm, and one-way leakage alone cannot
   alter the asserted two-sided absolute `P`-action formula;
5. the optional `C^1` `O(N^-2)` rational-ray shift cancels the leading
   `O(N^-1)` curvature mismatch against the `Theta(N)` derivative; and
6. the stronger `C^0` construction places the two physical wave numbers at
   `k_*+/-sqrt(bar e_N)` and makes `e_N/h_N->0`, the IVT sign scale.

The earlier `symbolic.*` receipt is retained as the preregistration-stage
one-way formula. `symbolic-v2.*` records the first two-sided repair,
`symbolic-v3.*` the optional differentiable rational-ray refinement, and v4
the final `C^0` bracketing check. None of the earlier receipts is used as the
final oracle.

## Artifact consistency

The final repository-interpreter command in `validation-v2.command.txt` exits
`0` and
checks the frozen README hash and activation receipt, YAML parsing, terminal
route-scoped verdict vocabulary, active parent state, the `C^0` rational-ray
and return-kernel result fields, and the corrected symbolic receipt.

The first development invocation is retained in `validation-first.*`: it
exited `1` because the validator looked for the notation `O(N^-2)` while the
result statement used the equivalent `O(1/N^2)`.  The predicate was corrected
to the artifact's literal notation and the successful run was captured
without changing any scientific result.

## Analytic checks not delegated to the oracle

- For the generalized Sturm pencil, differentiating the normalized simple
  eigenvalue gives (4), and substituting it into `d(k^2/lambda)/d(k^2)` gives
  the positive expression (5).
- Differentiating `sigma_1(k)=sigma_2(r_12(k)k)` proves (10c). Expanding that
  identity through cubic order gives (10d).
- Uniform graph-Riesz error `e_N->0` and `h_N=sqrt(bar e_N)` give
  `e_N=o(h_N)`. The simple column root has opposite signs at
  `k_*+/-h_N`; conditional on a continuous exact-carrier path covering the
  bracket, IVT gives (10i). Uniform error alone does not supply that path.
- The polarization alternative compares `A_J k+O(k^3)` with
  `C_b r^2 k^2 log(1/(rka))+O(r^2k^2)`. For fixed large rational `r` the
  first dominates positively as `k->0`, while at `k=x_0/r` the fixed positive
  bending term dominates. A sign-changing boundary root therefore has fixed
  opposite-sign endpoints; no derivative is inferred from the value-level
  remainder. Positivity of the bending Krein sign is not inferred from this
  frequency calculation.
- Under the explicitly conditional uniform physical KKS-dual normalization,
  three unit pure Fourier factors and one physical volume integration leave
  `(R a_c^2)^(-1/2)`. Since `delta ell_N=O(1)`, covariant angular derivatives
  add no `N^s` loss; the corresponding gate-time lower bound is (18f).
- Duhamel substitution gives (27). If `kappa_back<1`, the Volterra inequality
  gives (28a). The isolated Bogoliubov propagator and its inverse are both
  bounded by `S_A`; combining that lower/upper estimate with return closeness
  gives the absolute squared-action bound (28c). The presence of `PVQ` in
  `kappa_back` is essential.

## Claim boundary

This validation does not evaluate the Sturm nonidentity coefficient, prove a
continuous exact-carrier path covering geometric `delta`, prove the continuous
dual functionals `G_12` or `G_3` nonzero, derive an expanded `ad/ad^*`
representative for them, derive uniform physical KKS-dual normalization, build a gate
history, or control its nonlinear Euler remainder. It therefore does not
validate an equal-frequency doublet, an autonomous `SU(2)` analyzer, P2/P4,
or any quantum or particle claim. No production numerics were used.
