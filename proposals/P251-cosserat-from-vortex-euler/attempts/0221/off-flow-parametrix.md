# Why the remote pressure response has arbitrarily high decay order

This supplies the analytic mechanism behind the finite-wave proof's
off-flow observation bound, not merely a compact-remainder citation.
The background and fixed regular source bands are those already frozen
there. The observed stationary tag has a positive invariant psi gap
from those bands. All times lie in one fixed compact interval.

The full linear Euler generator on solenoidal velocities is
`L_K=-u.grad+A_K`, where A_K is a matrix pseudodifferential operator of
order zero. This follows directly by commuting the full Leray projector
past u.grad: its first-order symbol is the scalar transport symbol and
the commutator is order zero. Multiplication by iK.u, Du, and the actual
pressure symbol are included in A_K. The physical mean projector has
the fixed ray convention; its finite-dimensional terms are smoothing.
Appending the Lin equation gives the same scalar transport principal
symbol on the displacement/velocity block, with order-zero coupling.

For orientation, Shvydkoy's primary paper states this Euler form in
equations(4)--(5) and constructs a transported order-zero parametrix
with an order-minus-one residual in the proof of Proposition4.2,
equations(65)--(71):
[The essential spectrum of advective equations](https://arxiv.org/pdf/math-ph/0412019).
That paper's compact-remainder statement alone would not imply the
arbitrary-power estimate needed here. The following finite-order
continuation of its symbol construction supplies that stronger step.

Let U_t be scalar pullback by the actual steady flow. Conjugating the
system by U_t yields a time-dependent order-zero matrix PDO B_K(t).
On a bounded time interval, all of its symbol seminorms needed at any
fixed finite order are bounded. Construct a parametrix Q_M(t) of its
evolution with symbol `q_0+q_-1+...+q_-(M-1)`. The leading symbol solves
the finite-dimensional matrix transport ODE with identity initial data.
At each subsequent order the symbol composition formula gives a linear
matrix ODE with the same leading coefficient and a known source from
the preceding symbols. Variation of constants therefore constructs it
with zero initial correction. No inverse spectral gap is used.

The resulting `U_t Q_M(t)` has residual `U_t R_-M(t)`, where R_-M has
order at most -M and is bounded uniformly in t between the corresponding
Sobolev spaces. Duhamel with the actual Euler/Lin evolution gives a
remainder mapping H^(s-M) to H^s. Differentiating a fixed finite number
of times in t or K requires only a larger finite symbol order M and
the actual finite-interval Sobolev energy bounds. These operations do
not change the principal flow graph.

Multiply on the left by the tag cutoff and on the right by a cutoff
covering the fixed source bands. Since the flow preserves psi, the two
cutoffs are separated from that graph for every selected t. The kernel
of the transported PDO term is then smooth: away from its diagonal,
integration by parts in the Fourier covector gives every finite kernel
derivative bound. The Duhamel remainder gains M derivatives, with M
arbitrary. Consequently the complete cut-off propagator maps H^-q to
any prescribed finite H^s, at a fixed constant depending on q,s,T and
the source/tag gap. This is the precise smoothing statement being used.

The initial source `G(c) exp(iNtheta)` has nonvanishing phase gradient
on the fixed regular band. Its negative Sobolev norms decay as any
chosen inverse power of N, after including the actual derivatives and
normalization. The full pressure projection of the initial Kelvin
force may be noncompact, but it is itself a PDO: composing it into the
same parametrix preserves the flow graph and only changes the finite
order of the source map. The negative-sector order-zero tagged response
is therefore small, not identically zero.

The full-fluid mean/current uses propagated smooth adjoint tests rather
than a local tag cutoff. Finite-time adjoint Sobolev bounds keep those
tests smooth to every selected finite order; source-angle integration
by parts gives the same inverse-power estimate. The bounded ray-wise
mean projector is included, uniformly over the compact sphere of rays.
No differentiability of that projector across all Cartesian K at zero
is asserted. All K jets here follow the radial convention already
frozen for the common laboratory wave vector.

These arguments establish the stated finite list of remote physical
observations with arbitrary chosen inverse-power accuracy, including
the generally nonzero order-zero pressure tail. They do not normalize
an absent clock or construct a missing acoustic/optical cross-output.
