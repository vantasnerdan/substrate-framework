# Actual invisible phase control with retained full Euler energy

Use the SAME fixed C016 cell, and fixed smooth regular wrapped bands
separated from the actual observed core by a positive psi gap. The full
Euler/Lin and off-flow pressure statements of0205 are used at exactly
that scope. This supplies a new phase control, not an optical clock.

## Exact phase and energy on the actual passive sector

Let T=omega(c)partial_theta, with omega bounded above and below by
positive constants on these fixed bands. Select a real smooth G(c),
compactly supported in them and nonzero. For an integer N>=1 put

    g1=G cos(Ntheta),       g2=G sin(Ntheta),
    h1=G sin(Ntheta)/(2N omega),
    h2=-G cos(Ntheta)/(2N omega).

All fields multiply the physical e_X vector. They have zero mean under
the ACTUAL coarea measure mu(c)dc dtheta/(2pi). The exact histories are
w_i=e^(-tT)g_i e_X, eta_i=e^(-tT)(h_i+t g_i)e_X, not two assigned
oscillators. Their initial full cotangents are rho g_i e_X. Since
Th_i=g_i/2, their complete initial Jacobi energy MATRIX vanishes:

    H_ij=rho[<g_i,g_j>-<g_i,Th_j>-<Th_i,g_j>]=0.

Their full initial phase is instead

    Omega12=rho<h1 g2-g1 h2>
            =rho integral mu(c)G(c)^2/(2N omega(c))dc=:B_N>0.    (1)

The actual full linear Euler action conserves these restricted phase
and energy pairings along its solution columns. The fields themselves
need not span an invariant two-dimensional subspace. Their material
configuration and velocity do not vanish, and their embedding connection
is retained when restricting an action to physical observations.

To add a prescribed nonzero scalar phase b, multiply both columns by
sqrt(|b|/B_N); for b<0 interchange the two columns. This makes the added
phase exactly b J, with zero added K=0 energy. For b=0 use zero fields.
No negative density or probability is introduced. Because the bands and
G are fixed, B_N=B_1/N. Actual normalized velocity norm is O(sqrt(N)),
configuration norm O(N^-1/2), and Sobolev norms grow polynomially.
These are retained preparation costs, not a zero-cost normalization.

## Full finite-K preparation and the genuine linear current term

Use eta_i,K=P_K(h_i e_X), projected canonical cotangent rho P_K(g_i e_X),
and the exact Euler velocity w_i,K=P_K(g_i e_X-Du eta_i,K). As in0205,
the constrained longitudinal cotangent is restored when needed; it pairs
to zero with solenoidal eta but is not deleted from the Euler equation.

At X Fourier index0, for nonzero transverse n,

    (P_K)XX=1-KX²/|n|²+O(|K|³/|n|³),
    (partial_K P_K)e_X is transverse and of order -1,
    partial_K² P_K is of order -2.

Thus there is no first-K phase correction, and the second correction
to(1) is O(|b| |K|²/N²). The mean harmonic remains zero for these
initial fields. Smooth fixed bands have a nonstationary angle phase:
Ttheta=omega>0. Pseudodifferential integration by parts gives
||G exp(iNtheta)||_(H^-j)<=C_j N^-j for every fixed j. This is an
analytic high-harmonic bound with fixed profile constants, not a sampled
small eigenvalue or a numerical floor.

The energy deserves a separate calculation. With A_K=u.grad+iK.u,
the exact initial Hamiltonian is

    H=rho/2[||P_K g||²-2 Re<P_K g,A_K eta>
                -||(Id-P_K)A_K eta||²+<eta,Hess(p0)eta>].       (2)

It follows by using eta_t=P_K(g-A_K eta) in the actual Jacobi energy.
At K=0 equation(2) vanishes as proved above. Its first-K Hermitian
matrix is the actual transport/Noether term

    H_ij^(1)=i rho integral (K.u)(h_i g_j-g_i h_j).              (3)

This term is generally NONZERO and is retained, not called a small
pressure error. The nonzero wrapped-streamline drift can enter it.
For example, the phase-weighted background mean includes both the axial
speed c and the wrapped transverse drift. No material-spin interpretation
is attached to that vector.

To use the control for an axial physical optical input, transform its
entire field, tag and phase data by the same whole-field O(3) law as the
optical supplier. Pair R with -R, with positive half-weights. The measured
axial input axis det(R)R e_X is identical on this inversion pair. At fixed
LAB K the body wavevector changes sign, while the source body u stays
fixed; equivalently the laboratory polar vector coefficient changes sign.
Thus(3) cancels between the two realizations. Both
phases bJ agree, and neither energy is subtracted. Equivalently its
isotropic rank-three tensor is zero. This is cancellation in the full
declared whole-field law, not a per-realization zero-current claim.

After the normalization in(1), all second-K terms in(2) are O(|b|/N).
Indeed g has norm O(sqrt(N)), h has norm O(N^-1/2), first projection
derivatives lower the oscillatory order by one and second derivatives
by two. A_0 h=g/2; commuting transport through a projection derivative
retains its negative order. The explicit iK.u h terms and all first-
derivative products are O(N^-1) at second order. Hess(p0)e_X=0, so its
projected terms are no larger. The restored constraint term obeys the
same bound. These estimates include the low Fourier tails, which are
superalgebraically small on the fixed smooth band. Hence the inversion-
paired control has

    Omega=bJ+O(|b| |K|²/N²)+O_N(|K|³),
    H=O(|b| |K|²/N)+O_N(|K|³).                       (4)

Its leading phase can be finite without a persistent second-jet energy
cost. Each O_N remainder is finite only after N is chosen, and enters
the actual macro-scale diagonal as in C016.

The vector-angle convention remains Phi=3E[n theta], with full phase
and energy averaged WITHOUT that factor. Consequently a raw phase b
per family contributes b/3 to its isotropic vector phase. To change a
specified physical whole-law phase density by delta_j, use raw b=3delta_j
(and retain any declared family probability). This normalization does
not rescale the actual physical tag or its mechanical spin.

## Finite cross constraints and actual observed smallness

Initial cross phase/energy coefficients against any previously fixed
finite list of acoustic or optical preparation columns through degree
two are linear functionals of G for each of these two quadratures.
Collect both sets of finite tensor rows BEFORE whole-law averaging.
If there are M scalar rows, choose M+1 fixed smooth profiles on disjoint
regular bands and a nonzero vector in their homogeneous kernel. No
generic row independence is assumed. Unit coefficient normalization and
disjoint support bound the positive integral B_1 above and below. This
continues to hold if the kernel coefficients depend on N. Phase
normalization then grows at most sqrt(N) and preserves the exact cross
constraints. Fresh bands permit successive finite controls.

At K=0 these exact axial histories vanish on the transported core tag.
Every literal angle, spin, centroid, displacement G and symmetric shape
row there is unchanged. For first and second K derivatives, the full
Euler/Lin off-flow kernel from0205 is smooth between these separated
invariant regions. Repeated angular integration by parts gives an
O(N^-q) bound for any chosen q after absorbing the sqrt(N) normalization.
The same holds for a finite number of time derivatives, with constants
relative to the actual fixed nonzero tag reference and optical phase.
Complete Euler mean-current rows are also smooth linear observations
and obey this bound. No canonical momentum is substituted for those
physical currents.

The scale ordering can meet the first-gradient observation error and
the normalized cubic action remainder simultaneously. For the fixed
source profiles and finite Sobolev hierarchy, the latter constant is
bounded by C0 N^L for some finite L: all differentiated preparations have
polynomial growth. Choose an integration-by-parts order q>L+2 and then
K_N=c0 N^(-L-1), with fixed sufficiently small c0>0. The cubic remainder
divided by K_N² tends to zero. The first-K observation error divided by
K_N² is bounded by C N^(L+1-q), also tending to zero; second-K errors
and the energy/phase coefficients in(4) tend to zero directly. The
zeroth tagged error is exactly zero. This explicit window avoids a
hidden incompatible requirement that K be both too small and too large.
Constants include every previously fixed source/control norm and actual
tag denominator; no uniform simultaneous shrinkage of that tag is claimed.

Consequently this actual return can change inherited optical mass to a
separately known positive literal spin coefficient, while leaving the
physical clock and tag observations unchanged in the controlled limit.
It cannot construct that clock or guarantee that an independently
shrinking tag retains nonzero action density. The already established
phase-null energy controls may then adjust the conserved initial energy
to the same physical oscillator form, retaining all finite cross rows
and their own full current/gradient estimates. This is a constructive
phase/energy interface on one Euler field, not the completed parent law.
