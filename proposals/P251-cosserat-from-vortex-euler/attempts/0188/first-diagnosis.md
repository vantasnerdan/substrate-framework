# First execution diagnosis

The first execution passed the full stationary Euler and curl checks, then
failed at the determinant assertion. The determinant is
`3*A*B*(B-A)*(B+A)/16`, whereas the comparison target was left expanded as
`3*A*B*(B**2-A**2)/16`. Python structural equality is not polynomial equality.
The repair factors their difference and checks that it vanishes. Neither the
field, determinant, construction nor comparison value changed. The original
`first.stdout` is retained.
