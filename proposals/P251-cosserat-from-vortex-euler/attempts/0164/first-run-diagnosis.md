# First execution: expression-normal-form assertion repair

The first run passed all fourteen Euler, parity, nodal and curvature
checks, then exited1 at the physical mass comparison. Its derived entry
was `I*k**2/2-j*k**2/4+rho`; the asserted entry was
`rho+(I/2-j/4)*k**2`. These are exactly equal but have different SymPy
expression trees. The first output is preserved unchanged.

Repair: test exact simplified matrix DIFFERENCE equal to zero, using
the existing local `zero` helper, for both jet matrix comparisons.
No equation, target coefficient, parameter, sign, tolerance or scientific
predicate changes. The first failure is an assertion-representation
defect, not evidence against the acoustic current or its mass formula.
The repaired execution is captured separately in `repaired.stdout`.

That execution passed all seventeen general Euler/marker/matrix checks,
then was interrupted (exit130) inside the configuration API's second
fully generic symplectic inversion, where multivariate rational expression
swell was consuming the run. Its output and Python stack are preserved.
The implementation repair keeps the GENERAL current bracket as a direct
exact symbolic contraction, and calls the already-tested infrastructure
on one exact rational exposing instance while retaining the noncommuting
current parameter symbolically. This is an oracle-representation repair,
not a numerical approximation or a special-case proof of the general
construction. The new output is `final.stdout`.
