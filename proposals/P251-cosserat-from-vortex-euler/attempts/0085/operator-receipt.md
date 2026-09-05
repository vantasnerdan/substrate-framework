# Analytic operator construction before finite-jet computation

Let omega be real analytic and nonzero on a small ball. The compact
velocity equations are D xi=0, where

    D xi=(div xi, div(xi cross omega)).

Seek a three-component differential operator in right-normal form

    S_j=sum_{|alpha|<=N} partial^alpha composed with multiplication by s_jalpha.

In D S, commuting a known coefficient through a derivative differentiates
only the known coefficient. The unknown s_jalpha remains on the RIGHT.
Thus D S=0 is a finite homogeneous linear system over the field of
meromorphic analytic germs, not a differential equation for these
unknown coefficients. At N=6 it has 3 binom(9,3)=252 unknowns and at most
2 binom(10,3)=240 equations. Both constant-output rows vanish because
each row of D is an outer divergence, so its actual matrix has238 rows.
A nonzero null vector exists over that field. On a smaller open ball
avoiding its finitely many poles, it defines a nonzero analytic S.

For every compact smooth bump psi in this ball, xi=S psi and
v=xi cross omega are compact and divergence-free exactly. A nonzero S
cannot have v identically zero as an operator: that would give S=omega T
for a scalar differential operator T, and div S=(omega dot grad)T=0.
The nonzero principal symbols of differential operators over the germ
field multiply without zero divisors, so T=0, a contradiction.

This establishes nontrivial compact induced velocities, but neither
their angular moment nor a nonzero KKS pair. These require a separate
exact test. For K=e cross x and f=omega cross K,

    integral K dot v = integral psi S^* f,
    S^* f=sum_j,alpha (-1)^|alpha| s_jalpha partial^alpha f_j.

There are NO derivatives of the unknown s_jalpha in this adjoint row.
Consequently a finite-jet matrix can directly test whether the angular
row is outside the constraint row space. If so, a meromorphic null
vector has S^* f nonzero on a smaller ball, and a sign-definite compact
bump gives a nonzero exact mechanical spin. This would be an actual
construction, not an imposed modulus.

The first finite-jet probe is frozen as follows. N=6, unit curl eigenvalue,
eight prescribed rational unit wavevectors (the axes, three3-4-5
directions, and two1-2-2 directions), deterministic transverse rational
amplitudes, and their exact cosine/sine helicity partners. Derivatives at
the origin are rational. The initial oracle uses arithmetic modulo101;
it is a rigorous nonzero-minor test, not floating-point numerics. It earns
an angular-rank conclusion only if the unaugmented rank reaches a proved
characteristic-zero upper bound, or a separate exact rational rank proves
the bound. A lower modular rank by itself is not a no-spin verdict.
All three physical spin rows are retained and reported separately.

This local Beltrami prototype is a construction probe, not itself an EPS
knotted-tube background. Transferring an open finite-jet rank property to
an actual EPS good patch still requires the appropriate allowed local
approximation/support argument. No empirical comparator, discretized
PDE eigenvalue, or fitted constant enters this exact finite-jet oracle.
