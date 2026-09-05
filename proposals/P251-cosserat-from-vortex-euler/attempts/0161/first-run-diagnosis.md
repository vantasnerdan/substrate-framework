# First scientific execution: fixture assembly repair

The first execution passed the first fifteen predicates, then the computed
velocity covariance failed. Its actual wavevectors printed as
`(lambda,0), (-lambda/2,sqrt(3)*lambda**2/2),
(-lambda/2,-sqrt(3)*lambda**2/2)`, giving covariance
`diag(3*Psi**2*lambda**4/4,3*Psi**2*lambda**2/4)`.

Classification: implementation error. Python method-call precedence made
the external lambda multiply an already scaled second component of the
column join. The repair constructs the two Cartesian entries explicitly
and scales the complete vector once. No theorem coefficient, field,
criterion or numerical tolerance changes. The failed execution is retained
in first-run.stdout.

The second execution passed nineteen predicates, then a structural SymPy
equality compared expanded `rho*Ux*Vx+rho*Uy*Vy` with factored
`rho*(Ux*Vx+Uy*Vy)`. The exact residual is zero. The final predicate now
tests that expanded residual; second-run.stdout preserves the native
failure. This is an oracle representation repair, not a normalization
change or tolerance relaxation.
