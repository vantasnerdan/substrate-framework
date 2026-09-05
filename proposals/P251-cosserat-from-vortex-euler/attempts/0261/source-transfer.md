# Primary-source transfer for the degenerate inverse

This receipt records the exact theorem scope consumed by
`inverse-analysis.md`.  None of these sources is cited as a substitute for the
one-sided trace/right-inverse lemma isolated there.

## Local hypoellipticity

[Hörmander, *Acta Mathematica* 119 (1967), 147--171,
DOI 10.1007/BF02392081](https://doi.org/10.1007/BF02392081) proves local
hypoellipticity for second-order operators expressed as sums of squares plus a
drift when the required Lie algebra spans.  For the actual linearization in
0261, the diffusion fields vanish at `T=0`, but their first commutators with
the transverse drift are the two coordinate derivatives in equation (21) of
`inverse-analysis.md`.  The theorem therefore licenses local smoothness for a
distributional solution on a smooth two-sided extension of the collar.

It does not prove solvability, a boundary trace, Fredholm index, compatibility
with the physical quadratic-form domain, uniformity under a moving boundary,
or a global right inverse.

## Intrinsic local estimates with drift

[Bramanti--Zhu, *Analysis & PDE* 6 (2013), 1793--1855,
DOI 10.2140/apde.2013.6.1793](https://doi.org/10.2140/apde.2013.6.1793)
proves local Schauder estimates for `X_i X_j u` and `X_0u` when the real smooth
fields satisfy the paper's Hörmander hypothesis and the second-order
coefficient matrix is uniformly positive in the vector-field directions.  It
also proves the corresponding local `L^p` estimates under its VMO hypotheses.
The 0261 operator has identity coefficient matrix in the fields
`X_1=T d/dx`, `X_2=T d/dz` and the drift (18), so the algebraic
operator class and local bracket hypothesis match on a smooth extension.

The paper's estimates are local in an open domain.  They do not provide the
one-sided characteristic-boundary Poisson operator, the physical moderate
branch, the boundary-graph trace, or the global mode gluing required by (31).
Those exclusions are essential: applying an interior estimate directly up to
the unextended free boundary would change the theorem's hypotheses.

## Nonlinear inverse theorem

[Hamilton, *Bulletin of the AMS* 7 (1982), 65--222,
DOI 10.1090/S0273-0979-1982-15004-2](https://doi.org/10.1090/S0273-0979-1982-15004-2)
proves the Nash--Moser inverse theorem in tame Fréchet spaces.  Its relevant
license is conditional: a smooth tame nonlinear map whose nearby derivatives
possess a smooth tame family of inverses or right inverses admits the
corresponding local inverse or solution map.

Hamilton does not manufacture the tame right inverse.  In 0261, the family
required by that hypothesis is exactly the unproved estimate (34), including
the free-boundary trace and uniform dependence on the Hanzawa graph and
profile.  Invoking Nash--Moser before proving that family would be circular.

## Net license

The primary literature validates the local Hörmander interpretation and gives
the correct intrinsic regularity technology once a one-sided moderate solution
has been constructed.  It also supplies the final abstract nonlinear theorem
after a smooth tame right-inverse family exists.  No located theorem closes
the intervening characteristic-boundary trace/Fredholm problem for the exact
quadratic degeneracy and physical domain in 0261.
