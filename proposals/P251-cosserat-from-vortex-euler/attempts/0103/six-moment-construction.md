# Six independent compact angular response moments

Let A(x) be0085's right-normal constraint matrix for
div xi=div(xi cross omega)=0. Its universal characteristic-zero rank is at
most235, proved there from three independent highest-degree adjoint nulls.
The same rational eight-wave prototype attains that rank over F_101; none
of its denominators is divisible by101.

For an angular axis e define two actual functionals, including density,

    G_e(xi)=rho integral (e cross r).xi,
    L_e(xi)=rho integral (omega cross(e cross r)).xi.

Their right-normal adjoint coefficient rows are respectively

    R_G,e,jalpha=rho (-1)^|alpha| partial^alpha(e cross r)_j,
    R_L,e,jalpha=rho (-1)^|alpha| partial^alpha[omega cross(e cross r)]_j.

At the probe center r=0, the first row is supported only on |alpha|=1:
its value is minus the corresponding component of e cross e_alpha. The
second is exactly0085's physical-spin row. Multiplication by rho>0 does
not affect any rank.

The new exact arithmetic yields

    rank A=235,
    rank[A;R_L]=238,
    rank[A;R_G]=238,
    rank[A;R_L;R_G]=241.

The last nonzero minor modulo101 is a nonzero rational minor. Conversely,
the universal bound235 plus six added rows is at most241 in characteristic
zero. Thus the SIX rows are independent modulo the constraint row space
over the rationals, not merely over a numerical tolerance. This proves the
new independence at the prototype and on an open analytic-jet neighborhood.
The existing all-L-only upper bound cannot prove this new statement by
itself: replacing G by L loses three dimensions, as the mutation verifies.

Use the same235 pivot rows/columns of A as in0085. The remaining17 free
columns define analytic null coefficient sections s^l(x). Choose six whose
matrix J(x) of the six angular rows is invertible at the center. Invert J
on a smaller ball and right-multiply the null sections to obtain six
analytic coefficient sections with angular rows exactly the coordinate
basis. All constraints still vanish identically; differentiating unknown
coefficients is not omitted, because the operators remain in right-normal
form. For any compact smooth bump psi of integral1, form

    xi_l=sum_(j,alpha) e_j partial^alpha(s^l_jalpha psi).

Integration by parts gives exactly (L,G)(xi_l)=e_l. Both xi_l and its
actual induced velocity xi_l cross omega are compact and divergence free.
Their centroids and symmetric first moments vanish separately. Changing
the origin changes neither angular functional because both zeroth moments
vanish. The normalization may use the full adjoint rows on the support
ball; the simplified center formula is only the arithmetic witness.

The transfer to an actual EPS field is also constructive. Add epsilon times
the same translated, correctly scaled prototype to the existing EPS seed
at the same curl eigenvalue. In the chosen minor, A's first-divergence rows
and the G rows are independent of omega; the other rows are linear in its
jet. Its highest epsilon coefficient is the nonzero prototype minor.
Only finitely many epsilon values are excluded. An arbitrarily small good
epsilon preserves the EPS tube and the six-row rank simultaneously.

Rank openness supplies six disjoint smaller interior balls, each admitting
all six prescribed moments. In ball l choose only its l-th normalized
response. Distinct response supports have zero KKS and zero FULL H cross,
because their induced velocities are compact as well as their generators.
Every diagonal response energy is a finite computed number; its sign is
not assigned. These are actual six-moment adjustment directions on one
smooth Euler field, with their complete finite energy costs retained.

The finite geometry and strict minor bounds fit0098's good-patch event.
They can now adjust G_Q and G_S independently of L_Q and L_S while leaving
the core observation and raw cage pair unchanged on disjoint supports.
The nonlinear identity G_Q=B^2/P and its positivity hierarchy still require
0102's finite root construction; this rank theorem does not assume them.

Evidence: `joint_moments.py`, `first.stdout`, first execution status0 and
five passing checks. Its immutable prototype import prints the older
unrepaired row-count caution as historical output; the universal bound
above is the already established analytic repair, not inferred from that
old sentence. No new numerical eigenvalue or repeated rational-rank replay
was required to establish the six-row determinant argument.

Route verdict: established as stated. Next parent action: use these actual
responses in the full canonical phase action's physical-spin normalization,
then retain the complete second-gradient field map and reaction reduction.
