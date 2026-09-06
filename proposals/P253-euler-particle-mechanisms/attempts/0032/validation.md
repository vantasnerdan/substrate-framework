# 0032 analytic trace/metric validation receipt

## Boundary

Only activated attempt `0032` artifacts were written. No central proposal,
source, API, test, governance, memory, or other attempt was edited. No
numerical activation or production computation was performed.

## Exact oracle

`verify_trace_metric.py` checks five exact predicates: the full-pressure
Kelvin identity and transversality, polarization-area conservation, the
pressure-normal neutral amplitude, separation of hyperbolic/Jordan/`-I`
returns, and the positive-metric obstruction/classification. The first
execution's command, stdout, and exit are preserved alongside the verifier.

The append-only pullback verifier `verify_pullback_coefficient.py` independently
checks the Baldi metric entries, `det(D Phi)`, closure of the pulled cross
operator on the `kappa`-orthogonal basis, and the resonant `D(Omega)` upper
entry. Its first stdout and exit are preserved in `second-run.stdout.txt` and
`second-run.exit`.

The exposing failure inherited from 0025 is retained: treating `c=k cross A`
as a passive Cauchy vector omits `(vorticity dot k)A`. The corrected identity
is rederived here; no trace conclusion is based on the false shortcut.

## Analytic result and remaining construction

Pressure-localizability and dense irrational action-angle trajectories solve
the pressure-normal invariant sector with an explicit bounded positive metric.
For the derivative-resonant angle sector, conserved area proves only
`det(M)=1`. Frequencies do not determine the physical gradient and frame
connection, so the exact source trace remains the ordered integral of the
returned coefficient. Near-Jordan asymptotics do not decide its discriminant.

The append-only continuation now evaluates the pullback explicitly. With
`h=g^{-1}kappa`, `d=kappa^T h`, and `alpha=det(D Phi)/d`, the scalar coefficient
on `c_z=x(m,n,0)+y(0,0,1)` is
`[[mu alpha E_ee, a/m+mu alpha E_ef],[mu alpha E_fe,mu alpha E_ff]]`, with the
four `E` entries given as explicit contractions in `derivation.md`. The
small-shell trace correction is reduced to the finite scalar `T_1` Duhamel
integral and a controlled `O(I_0)` Taylor remainder. Its sign remains the
named exact analytic remainder; no generic Pexp label is being substituted for
that coefficient.

The route is therefore blocked with one named analytic dependency: evaluate
or sign-bound `Delta=(tr M)^2-4` for the source-defined ordered integral,
including the physical metric connection. If and only if that remainder is
irreducible after further exact reduction, a separately activated
small-ratio-numerics design may be proposed. No numerical verdict is claimed
here.

## Append-only continuation receipt (trace coefficient)

The explicit returned-frame matrices are
`C0=[[0,0],[-1/sqrt(2),0]]` and
`C1=[[3*cos(sigma)/2,sin(sigma)/sqrt(2)],[sin(sigma)/(2*sqrt(2)),-cos(sigma)]]`.
The exact ordered integral over `0..2*pi` is `diag(pi,-pi)`, so its trace is
zero. Hence `T1=0` and the justified bound is `tr(M)=2+O(I)`, not a claimed
nonzero `sqrt(I)` splitting. The first nonzero coefficient needs cubic source
Taylor data (`C2`) and is explicitly left as the next analytic dependency.

## Sol-High source-equation correction receipt

The preceding dependency statement is superseded. Baldi (4.29)--(4.33)
recursively fixes `W2`; no free cubic datum remains. The earlier `C1` omitted
rotation of the cylindrical orthonormal frame, so its matrix integral is kept
only as provenance and withdrawn as a physical coefficient.

`verify_c2_source_recurrence.py` derives `W2`, the inverse-action coefficients,
the cubic physical trajectory, and `C0,C1,C2` from Baldi (3.39), (3.48),
(4.29)--(4.33) and Gavrilov (4.2)--(4.4). Its exact ordered integrals are
`9*pi^2` from the single `C2` term and `2*pi^2` from the double `C1` term. On
the flat cutoff plateau the period correction is `O(I^2)` and cannot change
the order-`I` trace. Thus `tr(M)=2+22*pi^2*I+O(I^(3/2))` and
`Delta=88*pi^2*I+O(I^(3/2))`.

The development transcript contained two interrupted symbolic expansions and
one recurrence-transcription assertion failure before the final source-solved
verifier. They were implementation/OCR diagnostics, not scientific attempts.
The first successful frozen execution is preserved in `c2-first-success.*`.

## Bounded 0038 oracle correction receipt

The previous execution remains preserved, but its printed `C2[1,0]` omitted a
spatial cubic in the physical velocity gradient.  Because that term is
multiplied by `H'(c)/sqrt(H(c))=O(epsilon^-1)`, the order-correct truncation is
`mp(z,4)`.  The repaired verifier compares against its former spatial-degree
two truncation and asserts the exact restored increment

    sqrt(2)*(3*cos(q)^4/2-13*cos(q)^2/8-9/16).

The first corrected execution is captured verbatim in
`c2-corrected-spatial-jet.command.txt`,
`c2-corrected-spatial-jet.stdout.txt`, and
`c2-corrected-spatial-jet.exit`.  It exits `0` with all eight exact checks
passing.  The correction changes only `C2[1,0]`; the trace oracle again returns
`single_C2_trace=9*pi^2`, `double_C1_trace=2*pi^2`, and
`I_discriminant=88*pi^2`.  No numerical computation or unchanged downstream
oracle was run.

Because the first corrected command resolved to the ambient Anaconda
interpreter, the oracle was executed once more with the explicit repository
interpreter.  The `c2-repository-interpreter.*` quartet captures command,
stdout, empty stderr, and exit `0`; stdout is byte-identical to the first
corrected run.  `c2-spatial-jet-correction-receipt.md` pins the requested
pre-verifier, pre-stdout, post-verifier, and both corrected-run hash sets.
