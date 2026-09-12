# Validation and oracle scope

## Activation boundary

The frozen README has SHA-256
`f1d54e9aa4400ca663b732e46f7795f561f0b336d9511e0b5623d992d7e3da08`.
Central activation used the repository schema command whose receipt hash is
`0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8`;
its stdout hash is
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`
and `activation-schema.exit` contains exactly `0`.

## Exact executable oracle

The repository interpreter executed

    python proposals/P253-euler-particle-mechanisms/attempts/0062/verify_block_topology.py

once at its frozen predicate boundary. The exact command, stdout, empty
stderr, and exit `0` are materialized as `block-topology.command.txt`,
`block-topology.stdout`, `block-topology.stderr`, and
`block-topology.exit`. The run is exact SymPy algebra, not production
numerics; no long process or `w4:p6` allocation was needed.

The thirteen assertions derive and expose:

1. all three coordinate components of
   `A eta=-[W,eta]-[B eta,omega]`, including the Hodge signs;
2. the three-component nonzero-harmonic DA map;
3. its exact positive-core inverse including the `grad zeta` correction;
4. tangency of the `F(zeta)` centralizer field;
5. cylindrical divergence-free closure of that field;
6. the physical KKS density `r^2 zeta` before the `2 pi rho_0` factor;
7. exact divergence cancellation for the `J`-weighted Hamiltonian inner
   generator;
8. the `H^s`, `H^(s+1)`, and induced-vorticity `H^(s-1)` cutoff exponents;
9. opposite Hessian determinants at the two pendulum critical phases.

The oracle is sign sensitive. Reversing the Hodge term in any component makes
the computed Lie bracket disagree with the displayed block. Omitting the
`grad zeta` term prevents the DA inverse from recovering the theta component.
Removing the cylindrical factor `r` changes both the centralizer divergence
and KKS density. Reversing one pendulum determinant destroys their exact
opposition.

## Proof-level checks outside finite algebra

The following conclusions are exact analytic deductions narrated in
`derivation.md` rather than finite symbolic substitutes:

- Jacobi and `[W_0,omega_0]=0` give `A_0 C_0=C_0 M_0`.
- Global Hodge decomposition places the exterior circulation field in the
  zero toroidal harmonic and requires full transmission in every block.
- A smooth push-forward `W=g_*Y` conjugates the complete flows and hence
  preserves orbit and separatrix incidence.
- Direct rescaling of a nontrivial width-`h` cutoff gives the Sobolev lower
  exponents checked algebraically by the oracle.

## Maximum verdict

The oracle and derivation validate the harmonic differential expression,
whole-space representation, positive-core DA inverse, carrier-map
intertwiner, KKS measure, centralizer topology, exact-volume generator, and
self-similar cutoff scaling.

They do not validate:

- the reverse boundary Hardy estimate or closed global DA graph;
- bounded whole-space core/exterior Hodge transmission in that graph;
- the two-index commutator estimates, `HJ2`, or `GR`;
- a limiting-absorption boundary value or distorted adjoint trace;
- existence of the spectral seed `q`, the exact first critical pair, or an
  evaluated physical `V_*`;
- a convergent rotating Cao branch, stability, or a particle/quantum claim.

This separation prevents a passing finite algebra tally from being used as a
continuum resolvent or nonlinear existence oracle.

## 0063 bounded correction validation

The existing thirteen-predicate oracle was not rerun: none of its predicates
contains the coordinate-basis realization of the Biot--Savart integral, the
free-boundary collar family, or a claimed spectral seed. The correction is
instead exposed by exact analytic checks in `derivation.md`:

- conversion of cylindrical coordinate coefficients to and from an
  orthonormal Cartesian frame supplies the two `D` factors in (7a);
- the h-by-h support measure and one-derivative cost reproduce each exponent
  in (18b), (18e), (18g), and (18h);
- fixed `|n|>=2` character orthogonality, positive-core injectivity of `C_0`,
  and support separation from exterior stabilizers give the quotient lower
  bound (18c);
- `C_0 xi=curl(xi cross omega_0)` gives
  `B_R3 C_0 xi=P_L(xi cross omega_0)` and explicitly resolves the zero-moment
  low-frequency row without invoking an unconditional inhomogeneous Hodge
  gain.

The correction therefore changes the continuum route verdict: the ordinary
ambient and generator graph lower bounds are refuted at `s=4`, `p>=6`, while
the positive-core inverse, weighted orbit topology, and source-specific
sandwiched trace remain separate live constructions. It introduces no
`q_*`, `V_*`, branch, or stability claim. Exact pre/post hashes and unchanged
oracle hashes are recorded in `0063-bounded-correction-receipt.md`.
