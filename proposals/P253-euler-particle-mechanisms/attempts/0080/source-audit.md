# Source and dependency audit

## Primary source

Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2, cached at
`/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.

Used at exact scope:

- (2.1)--(2.6): the axisymmetric operator, symmetric Green kernel, measure
  `dnu=r dr dz`, local logarithm, far decay, and compact-source bounds;
- (3.1)--(3.4): the exact steady equation, circulation, concentration, axis,
  and decay rows;
- Lemma 3.8: the limiting two-translation kernel;
- Lemma 3.9: the normalized contradiction estimate, supplemented here by an
  explicit compact Fredholm range argument;
- Corollary 3.12 and the local Pohozaev calculation: value-level checks for
  the leading rows;
- Appendix C: the symmetry-reduced limiting Green/Fredholm mechanism.

Not imported:

- auxiliary system (3.36) is not an exact steady row;
- Proposition 3.13 does not give the missing C1 differentiated remainder;
- Theorem 1.6 uniqueness is not derivative injectivity;
- source orbital stability does not transfer to the charged extension.

## Framework inputs

P253/0068 and independent P253/0071 supply the transported current, joint
action, local constrained flow, and affine Coulomb domain. P253/0077 and
independent P253/0081 supply the exact traveling Maxwell elimination, modified
Grad--Shafranov/Bernoulli rows, tag stabilizer, and `tau=g^2` parity while
leaving branch existence open. The P253/0063 unweighted displacement-graph
counterexample does not obstruct this fixed-source steady Green space. No
dynamical Kelvin zero resolvent or P253/0074 spectral inverse is used.

## Derived here

The supported `C^(2,alpha)` Green map, explicit right-cell/adjoint split,
complement estimate, exact Lane--Emden scaling density, rescaled physical
moment matrix, positive constant-linear affine matching matrix, exact augmented
circulation--impulse bordered isomorphism, distinct circulation--mean-radius
bordered isomorphism, and subluminal-window charged IFT are new author results.
The extended-growth Lane--Emden decomposition uses the source Green
normalization to retain the scaling `m=0` cell and affine/decaying `m=1` row;
Cao Lemma 3.8 removes every additional bounded/decaying harmonic, while
evenness removes axial translation.
The two bordered realizations share the same Fredholm/blow-up proof and have
different finite rows; they are not asserted to be identical operators. The
impulse and first-moment rows use the frozen base targets
`kappa_*`, `I_*`, and `R_*=M_1(zeta_0)/kappa_*`; no functional
`R_mean(zeta)` is differentiated as though it were a constant. The
contradiction uses the scaled source norm divided by `H_epsilon^p`, which is
uniformly equivalent to the pulled-back `C^(2,alpha)` source norm; the raw
pointwise density scale is `epsilon_core^(-2)H_epsilon^p`. At
fixed circulation and exact mean radius the external-`epsilon_core` IFT gives
a local continuous uncharged path and the value relation
`s_epsilon/R=C(kappa,R)epsilon_core(1+o(1))`. A differentiated scale remainder
or direct two-sided coverage estimate is still required before that path earns
the geometric-delta bracket used by P253/0079. The projected-complement HSE is
the `Q A Q` isomorphism proved in (17). The full fixed-`(mu,c)`
`A_epsilon` inverse and the sharp two-by-two Schur sign remain open on the
separate full-kernel/LS-tangent and differentiated-remainder constructions.
The reusable module checks only the leading finite-dimensional algebra.

P253/0079 additionally shows why its uncharged `N->infinity` construction is
not automatically a charged construction: `delta_N=Theta(N^(-1))` makes the
Cao speed grow like `log N`, while the Maxwell block has the finite subluminal
ceiling (37). A charged crossing must prove a nonempty interval above all
graph, IVT, and response thresholds and below that ceiling.

The charged branch is explicitly an Euler--Maxwell extension consuming the
new `U(1)` field and coefficients `g`, `epsilon_EM`, and `mu_EM`. Nothing here
derives or adopts those inputs as bare-Euler physics, and the continuous action
scaling established in P253/0055 remains unaffected.
