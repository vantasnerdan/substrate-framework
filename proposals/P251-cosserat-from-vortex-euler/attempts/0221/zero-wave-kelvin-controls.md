# Exact zero-wave phase and energy controls on the actual Kelvin orbit

Use the SAME C016 field and the complete0218 coadjoint forms, with the
actual stationary observed core fixed first. The source bands here are
fixed regular wrapped streamlines outside that core and separated from
it by a positive psi gap. The full negative-sector velocity has pressure
tails; zero remote observation is not assumed for that sector.

## Actual opposite signatures and their uniform margins

Let T=omega(c)partial_theta on a fixed regular band with omega bounded
above and below by positive constants. Choose smooth real G(c) there,
and two real streams

    s1=G cos(Ntheta), s2=G sin(Ntheta), N>=1.

Impose the finite first-shell moments P1s1=P1s2=0 when using the negative
sector. They make the planar mean of its full Kelvin velocity zero, as
required by0218. Let brackets below use the ACTUAL invariant coarea
measure mu(c)dc dtheta/(2pi), including the chosen cell normalization.
Define

    E_N=rho N^2/2 integral mu omega^2 G^2,
    P_N=rho N/2 integral mu omega G^2.

The h=-s positive sector has exactly

    H_+=E_N I, Omega_+=-P_N J,

where H is the Hessian matrix, so scalar energy is z^T H z/2.
The h=0 negative sector has exactly

    H_-=-E_N I+R_N, Omega_-=+P_N J,
    (R_N)_ij=rho <Tsi,H^-1 Tsj>.                         (1)

R_N is positive semidefinite and uniformly O(||G||^2), because
H^-1/2 T is an order-zero bounded operator on this fixed smooth cell.
Thus -H_- is strictly positive definite for all sufficiently large N,
with eigenvalues E_N[1+O(N^-2)]. These are exact physical signatures,
not eigenvalues of a numerical truncation.

Choose orientation-preserving real matrices C_+,C_- from the ACTUAL
definite forms so C_+^T H_+ C_+=I and C_-^T H_- C_-=-I. For example an
inverse-transpose Cholesky factor is explicit for a2x2 positive form.
Their transformed phase scalars are

    kappa_+=-P_N/E_N,
    kappa_-=+P_N/sqrt(det(-H_-)).                       (2)

Here kappa_-=N^-1[<omega G^2>/<omega^2 G^2>][1+O(N^-2)]
and kappa_+ is the same expression with a minus sign and no remainder.
The ratio <omega G^2>/<omega^2 G^2> lies strictly between reciprocal
frequency bounds of the chosen band. The actual C016 wrapped period is
nonconstant (in particular it diverges at the separatrix), so two FIXED
regular frequency bands with disjoint reciprocal intervals exist.
Assigning the positive and negative sectors to these bands gives

    |kappa_++kappa_-|>=c/N>0                           (3)

for all sufficiently large N, with either sign available by swapping
their band assignments. This does not require the negative field itself
to be a monochromatic Euler mode. Its full solution and embedding stay
in the actual Euler action.

## Exact finite cross constraints preserve those margins

Against any previously fixed finite list of preparation columns, initial
phase and energy cross rows are LINEAR in each new G. Include both
quadratures, real/imaginary parts when needed, first-shell moments, and
any previously installed controls among the rows. If there are M real
constraints, take M+1 fixed smooth profiles on disjoint subbands of the
same frequency interval and select a unit vector in their homogeneous
kernel. No generic rank or nonzero determinant is assumed for that kernel.

Disjoint supports give uniform positive lower and upper bounds for
||G||^2 and the weighted integrals in E_N,P_N. All fixed profile derivative
norms also stay bounded. Thus (1)--(3) hold even when the kernel vector
depends on N or on previous controls. Construct the positive pair first;
include its four real columns' bilinear rows before constructing the
negative pair. Subsequent whitening preserves all zero cross rows.
The actual Leray pressure tails are in the energy constraints, not
discarded by an assertion of disjoint velocity support.

## Arbitrary phase with zero complete energy

After making the two pairs mutually phase/energy orthogonal, add the two
whitened pairs with the SAME real amplitude a. Their complete forms are

    H_total=a^2(I-I)=0,
    Omega_total=a^2(kappa_++kappa_-)J.

For a prescribed finite nonzero b, choose the band assignment matching
its sign and a^2=|b|/|kappa_++kappa_-|. The resulting phase is EXACTLY
bJ and its full initial energy matrix is EXACTLY zero. For b=0 use zero
fields. This construction uses actual opposite-sign Kelvin energies,
not negative probability or a subtraction of kinetic fluid mass.

The normalization retains its physical cost. Since E_N~N^2 and a~sqrt(N),
the normalized stream has size N^-1/2, its transverse displacement and
axial Euler velocity have size N^1/2, and its axial displacement has size
N^-1/2. All fixed Sobolev norms grow polynomially. The large transverse
displacement is different from0210's axial-only control and must be
included in the finite-K pressure estimate; its jets are not copied.

## Independent arbitrary energy with zero phase

A single actual generator assigned to the existing phase amplitudes by
a real vector v has energy H0 vv^T and identically zero self phase.
Use a normalized h=-s generator for H0>0 and a sufficiently high-N h=0
generator for H0<0. Impose all baseline and previous-control cross rows
as above. Spectral decomposition of a prescribed real symmetric finite
energy correction then supplies its signed rank-one terms. These are
actual physical Kelvin energies on the same field. They change the
complete conserved energy without changing the restricted phase.

This completes the zero-wave algebraic normalization interface, including
actual pressure and cross constraints. Its REMOTE physical-current and
finite spatial action errors are the next part of0221, not assertions
earned by the zero-wave form alone. At K=0 the negative fields generally
already have small nonzero off-tag pressure response. The full Euler/Lin
propagator and its derivatives, the initial/current rows and actual
norms must bound that error relative to the chosen macro wave number.
