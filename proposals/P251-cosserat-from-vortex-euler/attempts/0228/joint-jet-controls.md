# Triangular construction of complete joint physical-action jets

This uses the ACTUAL same-cell Kelvin control families of reviewed0221.
Their physical output errors and normalization costs stay part of the
construction. Let z be any fixed finite list of real initial phase
amplitudes, including acoustic and optical configuration/rate columns.
K is the common laboratory wave vector. The desired finite phase and
energy jets are the difference between two specified full physical
action forms, not coefficients inserted into an Euler equation.

## Reality types and the available actual building blocks

For real Bloch preparations, the phase is skew Hermitian and the energy
Hermitian, with `F(-K)=conj(F(K))`. Write through degree2

    Omega=Omega0+i S1(K)+Omega2(K),
    H=H0+i A1(K)+H2(K),

where Omega0,Omega2,A1 are real skew matrices and H0,H2,S1 real symmetric.
The polynomials indicated by subscripts have their respective degrees.
This distinction retains allowed odd polar/axial coupling terms.

Reviewed0221 supplies arbitrarily many actual mutually cross-orthogonal
control blocks at K=0. Strengthen its finite linear cross constraints
to include every first and second K derivative against previous and
baseline columns. All these rows remain linear in a new source profile:
the true solenoidal generator and full Leray/curl operators depend
linearly on that profile. There are finitely many real rows. Select
M+1 fixed disjoint subprofiles and a unit homogeneous-kernel vector,
then whiten the actual nonzero forms. The reviewed uniform norm and
signature margins survive this stronger finite list. Select these
profiles ONCE for the whole K jet, not as an uncontrolled K-dependent
kernel. Subsequent polynomial amplitude maps preserve the cross zeros.

The two building blocks are therefore:

    P: two actual columns with Omega0=J, H0=0,
    E_sigma: one actual column with Omega0=0, H0=sigma, sigma=+1 or -1.

J=[[0,1],[-1,0]]. Either physical phase orientation is available by real
column interchange. Every block retains its actual higher K jets.
In particular the scalar energy block has real EVEN H(K) and purely
imaginary ODD Omega(K), by the two reality identities above. A matrix
phase block can have a nonzero imaginary-skew H1; it is not discarded.

Assign a block to z using a matrix R(K) with real even and imaginary
odd coefficients. Its actual forms are `R^* Omega R` and `R^* H R`.
It represents actual complex Bloch amplitudes with their real conjugate
partners, not negative density. The following sequence makes the jet
matching triangular.

## Order zero

Decompose the desired real skew phase correction into wedges of real
input vectors. A P block with constant rows (v^T,w^T) contributes
`v w^T-w v^T` to Omega0 and zero to H0. Decompose the desired real
symmetric energy correction into signed rank-one terms; an E_sigma
block with constant row v^T contributes sigma v v^T and zero phase0.
All finite source/control cross jets have already been set to zero.
Compute and retain the actual first and second jets these controls add.

## Order one: energy first, then phase

For an energy target `i K_a(v w^T-w v^T)`, take two independent actual
scalar energy blocks of signs + and -, with input rows

    R_+(K)=v^T+i K_a w^T,       R_-(K)=v^T.

Their constant energies cancel EXACTLY. Their first energy coefficient
is the desired wedge, since each scalar block's own H1 is zero. Their
constant phases are zero. Their induced phase1 is a known imaginary
symmetric rank-one row, `i[sigma_+(K)+sigma_-(K)]v v^T`, where sigma_±
are the ACTUAL scalar phase1 coefficients. Retain it in the next step.
Real wedges span every skew input matrix; treat each spatial direction.

After all energy1 corrections, decompose the remaining symmetric
phase1 polynomial into rank-one input rows. For a target i f(K) v v^T,
where f is a real linear polynomial, use one P block with rows

    R(K)=(v^T, i l(K) v^T).

If its actual intrinsic entry Omega_11(K)=i d(K)+O(K^3), choose
`l(K)=[f(K)-d(K)]/2`. Then its phase1 is EXACTLY i f(K) v v^T.
Its phase0 and energy0 vanish. Its energy1 also vanishes: H0 is the
zero matrix, while the diagonal entry of its imaginary-skew H1 is zero.
Thus this step does not undo the energy1 matching. It may change degree2
forms, which are retained for the final step. The complete linear
polynomial d, not just the targeted K_a coefficient, is subtracted.

This proves the nontrivial triangularity: energy1 controls can alter
phase1, but the correcting phase1 controls do not alter energy1.
No inverse generic controllability matrix or discarded current is used.

## Order two: all new amplitude rows vanish at zero

For a phase target `K_a K_b(v w^T-w v^T)`, use a P block with rows
`(K_a v^T,K_b w^T)`. To satisfy Bloch reality take both rows multiplied
by i; the common factor leaves its pulled-back forms unchanged. Its
phase2 is the desired wedge, its energy2 is zero, and it has no lower
forms. Its own higher jets start contributing only at degree3.

For an energy target `sigma (l.K)^2 v v^T`, use E_sigma with row
`i(l.K)v^T`. Its energy2 is the target, its phase2 is zero, and all
lower forms vanish. These span every real symmetric matrix-valued
quadratic polynomial: rank-one vectors span symmetric input matrices,
and `K_a K_b=[(K_a+K_b)^2-(K_a-K_b)^2]/4` supplies the mixed spatial
monomials with actual positive/negative energy blocks. Thus phase2 and
energy2 can be corrected independently after all their actual previous
contributions are computed. The full desired joint jet is now matched.

## Whole-field law, pressure and physical-error window

Each added field is an actual same-cell Kelvin preparation evolved by
full Euler/Lin. Its finite profile constraints include every original
acoustic/optical cross row through order2. Orthogonality is not inferred
from disjoint velocities, whose pressure tails are noncompact.

For an isotropic target, rotate the entire construction together with
its physical input representation: U is polar and Phi axial, with the
same representation repeated on their rate coordinates. The transformed
control input map uses the inverse physical representation, and its K
argument is R^-1 K. Positive whole-O(3) averaging therefore preserves
the covariant prescribed jet. It also retains any allowed odd mixed
polar/axial form. An axial-only R/-R cancellation is not imposed on a
different representation. Signed amplitudes and indefinite Jacobi
energies remain distinct from the ensemble's positive probabilities.

All operations are finite algebra, definite2x2 normalization and a
finite homogeneous profile kernel. For source jet coefficients bounded
by some fixed inverse power of a narrow-band width h, their amplitude
and K-derivative costs are bounded by finite powers of h^-1 and N. The
number of controls and constraints is finite and independent of h,N.
One may use fixed coordinate wedges and sums/differences of coordinate
rank-one matrices, rather than differentiating an eigenbasis at an
eigenvalue crossing. The chosen h-indexed profiles require no h-derivative
in the physical K jets.

The exact jet match removes any need to argue that an uncanceled first
form coefficient is smaller than K. For physical outputs, however, the
negative-control K=0 pressure tail is still generally nonzero. Apply
0221's full off-flow bound to each of the finite controls, including
the new polynomial amplitude and preparation norms. With N=h^-M fixed
first, bound the total cubic remainder and required norms by h^-D.
Use a fixed source moment order r>2D+2 and integration order q large
enough that every polynomial prefactor times N^-q is o(h^(2D+2)). Set
K=h^(D+1). The source/remote observation errors are then o(K^2), and
the actual cubic remainder divided by K^2 vanishes. As in0221, r and q
change fixed constants, not the h exponents of the fixed derivative list.

This last diagonal requires the stated polynomial source-cost bounds;
an abstract dense-range approximation without an accuracy-versus-norm
bound does not automatically meet it. The construction supplies the
full action-normalization interface on the reviewed cell. Missing
physical acoustic angle, optical current evolution, stationary EPS
identification or density remains actual supplier work, not a free
coefficient assignment licensed by this normalization theorem.
