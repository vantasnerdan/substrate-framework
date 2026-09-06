# Bounded material-domain correction from independent review 0036

Independent review `0036` found that equations (18)--(19), although locally
and variationally exact, were written as absolute whole-space integrals. The
accepted `0027` reference has `u_theta=L_infinity/r` throughout the infinite
axial exterior, so those absolute integrals diverge. The pre-correction
`canonical-construction.md` SHA-256 was
`5a72e4e5b7bf7a83d8d5add8b0f2e87f2d6bf961f4c9a57764f40137b259ce12`.

The correction states (18)--(19) on finite material domains or for compact
variations and adds the whole-space fixed-background expression

    R_rel=2pi rho_m int dt int a da db
      [(r_t^2+z_t^2)/2-xi^2/(2r^2)+xi^2/(2a^2)]

with the relative volume-constraint term and a declared finite-excess domain.
The added term is fixed label data, so the Euler equations, cyclic momentum,
Hamiltonian, and centrifugal sign are unchanged. The first-order action is
also scoped to fixed-end relative or compact variations. No API, test, or
source predicate changed, so no oracle was rerun. Independent `0036` will
check this single repair against the post-correction hash.
