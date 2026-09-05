# First execution: representation-only equality repair

The first execution passed the six orbit and physical-moment checks, then
failed the scalar ratio comparison. `cancel` represented its denominator
as 4 Omega^2-4 with the opposite numerator sign, while the expected
expression used 4(1-Omega^2). Structural Python equality is not rational
function equality. The minimum repair replaces that one predicate

    spin_ratio4 == rho*q*(1+w**2)/(4*(1-w**2))

with the cancelled difference equal to zero. No input, equation, profile,
coefficient, sign or tolerance changes. first.stdout is retained, and the
corrected first successful execution is captured separately.
