# Fixed physical cutoff: exact phase and curvature jets

This is a named constructive correction, first exposed while deriving0166.
It changes no frozen0155 files or C-CST-011 evidence. The latter's separate
high-order packet/moment hierarchy is not being reviewed here.

## The exact distinction

For a fixed physical radial label, an equality for its phase coefficient
at carrier p_star does not imply equality for that coefficient's carrier
derivatives. A single interior bump can fix q_chi(p_star)=q0 but generally
leaves nonzero q_chi'(p_star),q_chi''(p_star). They enter the leading
natural-scale curvature, not merely an O(delta²) remainder.

For reference chi=x and a Laguerre mode(n,m), put

    f(x)=exp(-x/2)L_n^(m-1)(x), E=2n+m, q0=m/E,
    R(x)=(x/2-E-q0)f(x), D0=(3/2)x partial_x.

Since x=r²/ell(p)² is proportional to p^(3/2) at fixed physical r,
D0 is its exact relative carrier derivative. Common scalar prefactors
cancel from the ratio. The phase equality through the two needed jets
is therefore equivalent to the THREE conditions

    F_j(chi)=integral chi x^(m-1) D0^j R(x) dx=0,
    j=0,1,2.                                             (1)

For untruncated chi=x all three vanish. The j=0 identity is the exact
Laguerre Laplace calculation. The other two also follow by integration
by parts against the homogeneous measure x^m dx, with vanishing
Gaussian endpoint terms. Fixing only F_0 does not fix F_1,F_2.

## Three finite bumps repair all three equations

Let chi0=x times a sufficiently distant smooth cutoff, and choose three
small disjoint interior bumps q_b on regions with chi0>0. The matrix

    M_jb=F_j(q_b)

can be made invertible. Indeed the three analytic densities in(1) are
independent: after removing the common Gaussian and positive power,
their polynomial degrees are n+1,n+2,n+3, with nonzero leading
coefficients. Independent continuous functions admit an invertible
point-evaluation minor; sufficiently small smooth bumps preserve it.

Set a=M^(-1)(F_0(chi0),F_1(chi0),F_2(chi0))^T and

    chi=chi0-sum_b a_b q_b.                               (2)

All three residuals in(1) are exponentially small as the cutoff grows,
while M is a fixed finite matrix. Thus one sufficiently large but
FINITE cutoff makes(2) smooth and nonnegative, with a nonzero angle
denominator. Rescaling chi by one common positive constant can enforce
chi<=1 without changing(1). Equation(2) is an explicit tag-geometry
construction, not an empirical fit or a free frequency parameter.

For0166,n=2,m=2, this gives exactly q0=1/3, q'=q''=0 at p_star and
the stated Omega² delta curvature. Without(2), the exact coefficient
must include the actual cutoff jets; positivity still follows from
a sufficiently small tail, but the ideal displayed coefficient is
not exact.

## Named correction to the0155 coefficient, not its positive mechanism

0155's `cutoff-coefficient-repair.md` imposes only the j=0 condition.
It establishes the exact value5/21 there but does not by itself make
the finite physical cutoff's curvature coefficient exactly5/7. The
strongest positive finite-cell mechanism survives: all three tail
errors can be made smaller than its strict positive coefficient.
For its exact5/7 assertion, apply(1)--(2) with n=8,m=5, giving
q0=5/21 and its two zero carrier derivatives. The same degree argument
proves the required three-bump matrix exists, and its eight pressure/
reference marking rows remain available on open intervals of chi>0.
This supplies the missing fixed-tag construction append-only, without
retroactively treating a one-value match as a three-jet identity.
