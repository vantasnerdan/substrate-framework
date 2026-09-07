# Source and dependency audit

## Cao--Lai--Qin--Zhan--Zou

Primary source: Cao et al., arXiv:2206.10165v2, cached at
`/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.

Used at these exact scopes:

- Proposition 1.4: existence of the exact thin polynomial vortex-ring family
  with fixed circulation and prescribed translating speed.
- Equations (3.1)--(3.4): the exact whole-space streamfunction equation,
  circulation row, concentration class, axis boundary, and decay.
- Lemma 3.8: the limiting Lane--Emden linearized kernel consists of the two
  translations.
- Lemma 3.9 and Appendix C: the contradiction/Fredholm mechanism for a
  normalized concentrated linear problem and the odd translation cell.
- Theorem 1.6: uniqueness inside the stated concentrated family, used only to
  identify the uncharged carrier and never as derivative invertibility.

Not imported:

- equation (3.36) is an auxiliary asymptotic parameter system.  Corollary
  3.12 and Proposition 3.13 compare exact parameters to it with nonzero
  errors.  It is not the exact augmented charged steady map and does not prove
  the Schur determinant (31);
- nonlinear orbital stability does not transfer to the Euler--Maxwell
  extension or its tag directions;
- the source proves no Maxwell, charged-carrier, KKS, or joint-orbit theorem.

## Reviewed framework inputs

- P253/0068 and independent P253/0071 establish the declared transported-tag
  U(1) extension, action, local constraint-preserving flow, static Coulomb
  sign, and retained-history map.  They explicitly leave a stationary charged
  carrier open.
- P253/0075 and independent P253/0076 establish local propagation of a fixed
  degree-minus-two inertial asymptotic coefficient.  The fixed-frame texture's
  stationary residual is nonzero.
- P253/0063 supplies the unweighted collar counterexample.  It prevents an
  unproved ordinary displacement/graph inverse from entering the persistence
  claim; it does not obstruct the steady elliptic map above.
- P253/0074 is active author work and is excluded as a premise.  In particular
  no dynamical zero-frequency inverse or Cao Riesz cluster is used.

## New derivation boundary

Equations (5)--(25) are derived directly from the reviewed joint action and
axisymmetric calculus.  The exact charged map is (21)--(23).  The source and
reviewed inputs support its uncharged base, elliptic Maxwell block, and limiting
kernel structure.  They do not construct the exact complement Fredholm/
surjectivity theorem HSE on the charged map's weighted spaces, and they do not
support the finite-core speed--impulse Schur determinant (31).  Consequently
the exact charged branch implication is conditional on both named analytic
achievements; it is not promoted as an existence theorem in this attempt.
