# Exact parameter-dependent normalization of both physical forms

This is the failure-derived continuation registered in README. Its input
blocks are the actual local Euler preparations constructed in
`local-normalizer.md`. The finite algebra below also applies independently
to any two forms with those blocks. The periodic Euler application retains
that file's actual-background and invariant source/tag gap hypotheses.

## 1. A regular simultaneous null frame

For n desired columns, take mutually orthogonal blocks for **both** forms:
energy vectors e_i^+, e_i^- with H(e_i^+,e_j^+)=delta_ij,
H(e_i^-,e_j^-)=-delta_ij and all their Omega entries zero; phase vectors
p_i,q_i with H zero and Omega(p_i,q_j)=delta_ij. Every other phase entry
vanishes. Here H is Hermitian, Omega is skew-Hermitian, and the real case
uses symmetric/skew-symmetric forms. The real physical construction in the
preceding file supplies precisely these blocks by additional finite cross
constraints. Complexification does not require new physical generators.

Use matrices of columns

    V_i=e_i^++e_i^-+p_i,
    W_i=(e_i^+-e_i^-)/2,
    Q_i=q_i.                                             (1)

Then V is null for both forms, the combined span of W and Q is null for
both forms, and

    V* H W=I,   V* H Q=0,
    V* Omega W=0,   V* Omega Q=I.                         (2)

For any Hermitian target h and skew-Hermitian target o, set

    Y=V+(W h+Q o)/2.                                     (3)

The exact identities are

    Y* H Y=h,       Y* Omega Y=o.                        (4)

Indeed the quadratic correction vanishes on the common null dual span.
The cross terms are (h+h*)/2 and (o-o*)/2. Equation (4) is an exact
identity, with no smallness condition on its targets. More generally, the
derivative of the joint Gram map at V has the explicit right inverse

    R(h,o)=(W h+Q o)/2.                                  (5)

An imaginary diagonal of o is included; restricting the complex form to a
real skew matrix would have omitted that Bloch degree of freedom.

## 2. The actual affine problem and a uniform local inverse

Let C(t) be the matrix of actual auxiliary Euler sources, with full forms
H_c(t), Omega_c(t), and let B(t) be the already selected observed physical
sources. Here t is the signed wave magnitude on one fixed unit ray; the
same argument works uniformly on a compact set of rays whenever the input
forms have uniform bounds there. Set

    H_b=B* H B,       X_H=B* H C,
    Omega_b=B* Omega B, X_O=B* Omega C.

All these are actual forms, including the full Leray projector. The total
source B+C L Y has forms

    H_b+L(X_H Y+Y* X_H*)+L^2 Y* H_c Y,
    Omega_b+L(X_O Y-Y* X_O*)+L^2 Y* Omega_c Y.             (6)

In particular disjoint supports have not been used to delete kinetic
pressure cross terms. Choose a fixed L sufficiently large that, throughout
the parameter neighborhood,

    ||X_H||/L+||X_O||/L <= eta,
    ||H_target-H_b||/L^2+||Omega_target-Omega_b||/L^2
        <= eta.                                         (7)

After division by L^2, (6) is a perturbation of the joint Gram map at V.
Write Y=V+R(z), restricting z to a pair consisting of a Hermitian and a
skew-Hermitian matrix. At t=0 and zero affine perturbation, this map is
exactly z by (4). If

    ||H_c(t)-H_c(0)||+||Omega_c(t)-Omega_c(0)|| <= eta,    (8)

its derivative on a fixed small ball in z differs from the identity by
at most C eta. Its value at z=0 is O(eta). Choose eta so C eta<1/2 and
use the contraction z -> z-[normalized forms-targets]. This gives a
unique small solution of (6) with both forms exactly equal to their
targets. Differentiating this finite-dimensional equation yields smooth
parameter dependence and explicit bounds in terms of the finite input
derivatives and (1-C eta)^-1. No eigenvalue tending to zero is inverted:
the right inverse is the explicit constant map (5).

Reality follows from uniqueness. When the source families and targets
satisfy the physical conjugation relation at -t, the same relation holds
for z and for the total sources. Thus one does not choose unrelated
normalizers for a Fourier mode and its conjugate.

## 3. Uniformity in the auxiliary carrier

The physical construction at t=0 uses finitely many profiles and constraints.
After actual whitening, phase blocks have generators of order one and
energy blocks have generators of order N_c^-1/2. The actual full forms at
t=0 are the fixed matrices in section 1. Before their positive/negative
cancellation, phase-block energies are O(N_c), but a derivative with
respect to Bloch wave number lowers this symbol order by one. Thus their
first parameter derivatives are O(1); higher fixed derivatives are bounded
as well. Energy-block derivatives are smaller. Finite profile constraints
have uniformly bounded smooth representatives: choose a unit kernel vector
in a fixed finite-dimensional family of disjoint bumps. There is no need
to differentiate that choice with respect to N_c.

The zero-force-mean constraints from `local-normalizer.md` remove the
initial projector singularity for the auxiliary columns. Consequently
(8) holds for |t|<=t_0 with t_0 independent of sufficiently large N_c.
The baseline columns may have a ray-wise mean projector. Along a fixed
ray its bounded angular matrix is retained in X and H_b; differentiating
the radial parameter does not differentiate an artificial singular mean.

For each already fixed smooth baseline family, all X derivatives have
finite bounds independent of large N_c. One can see this directly by
moving the auxiliary oscillation onto the smooth baseline test in the
full kinetic pairing. Local helicity and KKS terms also have these
bounds, whether retained as local integrals or eliminated by exact
support separation. Polynomially normalized control amplitudes are
included before using arbitrary-order integration by parts. Accordingly
L in (7) can depend on the baseline costs but is independent of large
N_c. The construction therefore has no requirement of the form
|t|<N_c^-D or |t|<N_s^-D. If baseline parameter derivatives grow with N_s,
their finite cost enters L and the differentiated solution; it does not
change the small neighborhood (8).

## 4. Observation preservation and the actual license

For a fixed physical supplier inventory, t, and L, choose N_c after those
quantities. Invariant source/tag separation gives arbitrary-order smoothing
of every retained remote finite-window observation, including its time and
parameter derivatives. Since the normalized auxiliary family and its
coefficient derivatives have only finite polynomial costs, N_c can make
its measured contribution smaller than any prescribed positive tolerance.
The nonzero Leray tails are included in this estimate. This ordering is
compatible with the baseline measured-jet ordering under investigation in
0250: the baseline carrier is selected first, then the auxiliary carrier.

The resulting full action forms match their targets exactly at the selected
wave numbers, throughout the neighborhood in which (6) is solved. There
is no action cubic Taylor remainder to compare with the observation
accuracy. This establishes the simultaneous normalization algebra and its
transfer under the local-block and off-flow hypotheses; it is stronger
than only matching two Taylor jets.

The same statement does not, by itself, identify the physical momentum and
angular currents with the prescribed target. That conclusion still needs
the actual 0246 Ward/current transfer for these parameter-dependent sources.
Likewise an isolated fixed ring is not a periodic positive-density compact
Euler assembly. Those are explicit remaining parent joins, not supplied
matrices or assumptions hidden in (4).

`route_verdict: established as stated for the exact simultaneous Gram
construction and its stated uniform physical-block transfer hypotheses`

`evidence_scope: exact finite algebra and analytic perturbation; full
current transfer and compact stationary carrier remain separate`
