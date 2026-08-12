# P229 promotion frontier for issue #59

This record separates the two narrow review candidates from attractive but
currently unsupported P229 conclusions. It is a handoff boundary, not a
negative campaign completion verdict.

## C-CMV-002 remains proposed

The current `mpmath` route differentiates the same zeta construction used by
the SymPy route, so it is regression coverage rather than an independent
two-loop master-formula derivation. The combinatoric prefactor, subtraction
scheme, and tachyon branch still need a genuinely independent derivation. The
one-loop-running substitution in `rg_improved_potential` also has opposite-sign
infinite one-sided limits at `x=b/Lambda^2=1` whenever a nonzero two-loop term
is selected. It cannot support a global two-loop potential or corrected
minimum. The next decisive action is an independent diagrammatic/master-
integral derivation with scheme and branch reconciliation before any new
minimum is claimed.

## C-CMV-004 remains proposed

The existing lattice work compares only square, triangular, and one rectangular
ansatz and does not establish a global optimum. More fundamentally, it samples
the homogeneous one-loop potential over an inhomogeneous classical field; it
does not compute the fluctuation determinant of that inhomogeneous background.
No amount of additional sampling repairs that missing requantization step. The
next decisive action is to construct the actual gauge-fixed fluctuation
operator for at least one periodic background, demonstrate convergence and an
independent soluble or numerical oracle, and only then compare preregistered
lattice candidates.

## Excluded narrative

P229's source audit, printed tables, literature summaries, ring-resummation
narrative, Preparata variational interpretation, and physical QCD/substrate
vacuum language remain provenance or research frontier. They are not separate
accepted claims and are not imported into C-CMV-001 or C-CMV-003.
