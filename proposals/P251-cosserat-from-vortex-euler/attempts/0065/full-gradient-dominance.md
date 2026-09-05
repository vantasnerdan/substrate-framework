# Positive curvature after the complete stationary reduction

## 1. Freeze the geometry and base operator before the gradient carrier

Use the normalized stationary actual-EPS-containing ensemble of 0057 and
the physical compact pair of 0059. First choose the entire finite coherence
geometry: physical core jets, raw base cage, eleven disjoint affine response
balls, and THREE additional disjoint cage balls away from every raw base
support and response ball. The additional balls lie where one local
vorticity component omega_z has a fixed positive lower bound. Choose their
cutoffs and three orthonormal reference bond directions t1,t2,t3, with
one fixed nonzero length d within the coherence scale.

Select the good-patch event including these bounds before any carrier
comparison. It has fixed positive intensity nu by the source-support
construction. Choose a finite BASE carrier k0 giving the complete positive
angle/reaction block, positive affine shear after reaction, and coercive
stationary momentum operator P. All base geometry, response functions,
amplitude conventions, operator jets and norm bounds are then fixed.

Only now choose the independent gradient carrier kappa, with
`kappa/lambda>0`, eventually large enough for the derived bound below.
Neither the good-patch event nor the normalization is reselected as kappa
changes. This hierarchy distinguishes actual carrier dominance from a
coefficient fitted to the desired dispersion.

## 2. The geometric gradient lift and its exact KKS orthogonality

In additional ball j take a SINGLE circular negative-helicity carrier

    zeta_j=-curl(phi_j p_j)/kappa,
    p_j=(cos(kappa z),sin(kappa z),0).

Project it with the SAME eleven dual responses eta^A as the base:
`Z_j=Pi_11 zeta_j=zeta_j-sum_A eta^A F_A(zeta_j)`.
The fixed response fields depend on the actual sampled background but not
on kappa. Their coefficients are uniformly bounded; integration by parts
against their fixed smooth moment kernels even gives inverse-carrier decay.
The original physical core jets are unchanged.

Let e be the physical local relative-angle axis. Its scalar physical
angle is `q_a(X)=e.q(X)`. The prescribed extra material displacement is

    sum_j [q_a(X+d t_j/2)-q_a(X-d t_j/2)] Z_j.

This is an independently testable neighboring-frame kinematic attachment.
No elastic energy is assigned to the difference. Evaluate its energy and
reaction using the full stationary Euler action. It vanishes identically
for uniform q and changes no zeroth-order physical angle, inertia or
locking coefficient. Marks and bonds are held fixed under variations.

The raw additional cages are disjoint from raw Q,S and from each other.
At each response ball all corrected fields are scalar multiples of the
ONE vector field eta^A, while distinct response balls are disjoint. Thus
their KKS crosses are zero POINTWISE, not merely by an averaged cancellation:

    Omega(Z_j,Q)=Omega(Z_j,S)=Omega(Z_j,Z_l)=0.

The eleven affine pairings also vanish exactly. This remains true when a
base column retains a prescribed affine moment by adding arbitrary fixed
coefficients of eta^A; zero base rotation moments are not needed for this
orthogonality argument.

For the selected patchwise Bloch amplitudes these identities hold at every
macro wave vector. The full momentum Hessian P(k) is unchanged because its
columns have not changed. Its entire KKS coupling D(k) is unchanged by the
displayed orthogonality. Therefore the complete reduced kinetic operator

    J(k)=D(k)^* P(k)^-1 D(k)

is unchanged, INCLUDING all momentum-gradient, mixed and mass-gradient
terms. This is an exact property of the chosen lift, not an instruction
to drop such terms. Affine symplectic columns are unaffected through the
second slow order; their first omitted macro jet is quadratic in k and
its pairing with this gradient lift is cubic.

## 3. Uniform mixed-jet bounds despite the large carrier

Let epsilon n be a macro wave vector, |n|=1, and use the stationary
representation of 0057. The extra configuration column is

    Q_lift(epsilon)=Q_base(epsilon)+i epsilon Z_n(epsilon)+O(epsilon³),
    Z_n=sum_j d (t_j.n) Z_j (e.q).

Here the stationary profile of each patch carries its exact factor
`exp(-i epsilon n.y)`; the centered bond symbol is
`2i sin(epsilon d t_j.n/2)`. Momentum fields remain independent stationary
reaction amplitudes. No isolated-patch inverse is substituted.

Write the full mixed Hessian between any fixed base column B and an added
column Z as

    H_epsilon(B,Z)=rho <V_B(epsilon),F_Z(epsilon)>,
    V_B=P_epsilon F_B-curl_epsilon F_B/lambda.

This identity moves BOTH Leray and curl onto the fixed base column. It
holds in the stationary inner product, including all interpatch tails.
The base response V_B is C^1 with a finite norm bound, by 0057's exact
moment/double-divergence construction. Its eventual curl-plus-double-div
extension also suffices: the curl part has a polynomial symbol.

The added force has a uniform carrier-independent L² bound A_Z. Its first
macro derivative is multiplication by bounded local y, so has a bound
R A_Z. No spatial derivative of its fast phase is taken here. Consequently

    ||H_0(B,Z)|| <=M_B0=rho ||V_B(0)|| A_Z,
    ||partial_epsilon H_0(B,Z)||
        <=M_B1=rho [||V'_B(0)||+R ||V_B(0)||] A_Z.

For reaction columns these are operator norms over ALL square-integrable
stationary patch amplitudes. Coercivity and the uniform fixed-column
Sobolev bounds establish the finite norms. They do not assume independent
kinetic energies of different patches.

## 4. Positive full gradient energy and the exact Schur remainder

The complete stationary Hessian is

    H=rho E[|P0 F|²-F.curl F/lambda].

Keep the first term as its full nonnegative Gram, including every cross.
The local helicity matrix of the three projected added cages obeys,
uniformly on the good event,

    -lambda^-1 integral F_Zj.curl F_Zl
      =(abs(kappa)/abs(lambda)) A_j delta_jl+R_jl,
    A_j=integral(phi_j omega_z)²>=A_*>0,
    ||R||<=C_*.

The remainder contains all compact returns and shared-response crosses.
Moving curl to each fixed response bounds its carrier cross independently
of kappa. Different raw cage supports have no local helicity cross.

For a coherent vector q and any n, the triad identity gives
`sum_j (t_j.n)²=1`. The full Haar rotation of the field and geometry gives
`E_Palm[(e.q)²]=|q|²/3`. Hence

    H(Z_n,Z_n)>= (rho nu d²/3)
                    [abs(kappa) A_*/abs(lambda)-C_*] |q|².

Let base momentum, mixed and retained Taylor coefficients be P_i,N_i,H_i,
i=0,1,2, with no factorial in coefficient two; put R_i=(P^-1)_i. The new
mixed block has

    Delta N_0=0, Delta N_1=L1, Delta N_2=L2,
    ||L1||<=M_P0, ||L2||<=M_P1.

The retained H_2 gains H(Z_n,Z_n) plus a mixed-derivative remainder bounded
by `2 M_Q1`. Its first-order change is bounded by `2 M_Q0`.
The complete Schur gradient coefficient is
`K_2=H_2-sum_{a+b+c=2} N_a^* R_b N_c`.
The change in its reaction subtraction is bounded by the explicit fixed
polynomial

    R_Schur=2||N0||||R0|| M_P1
       +2||N1||||R0|| M_P0+||R0|| M_P0²
       +2||N0||||R1|| M_P0.

No commutation or diagonal assumption is used. P_2 and R_2 remain in the
base coefficient; they cancel only from this DIFFERENCE because P and N0
are unchanged. The canonical `hermitian_schur_jet` verifies this full
noncommuting coefficient algebra.

## 5. Pair the full action, then normalize it

For the two circulation signs, vary the entire reaction fields separately
before averaging. P,H,N are even, D changes sign, and the odd gyroscopic
term cancels after the inverse operator is taken. Then average the full
reflected/rotated actions. The odd spin-gradient potential term is removed
by reflection, not by pretending its unsigned microscopic square is zero.
Mixed strain/relative-gradient terms are parity odd as well. Every even
reaction correction above remains.

The lift leaves the zeroth coefficients alpha,j,rho and the full kinetic
jet unchanged. Thus, in the 0059 relative-angle representation, the physical
normal-form correction `-alpha j/rho`, and any additional fixed kinetic-jet
corrections, are unchanged. They are included in the fixed base coefficient.
More generally, if a fixed kinetic normalizing map has first jet T1 and
one retains an unaveraged first potential jet, its extra difference is
bounded by `R_map=2 ||T1|| (2 M_Q0+2||N0||||R0||M_P0)`.
After the selected parity reduction this optional bound can be set to zero;
keeping it is a valid conservative bound for another fixed representation.

Let R_base bound the operator norm of the COMPLETE base normal-form
gradient symbol uniformly over |n|=1. It includes P_2,D_2, all mixed jets,
the relative kinetic correction, and any physical current-field map. All
are finite by the earlier exact second-jet construction. Define

    R_total=R_base+(rho nu d²/3) C_*+2 M_Q1+R_Schur+R_map.

Then both transverse and longitudinal isotropic physical normal-form
curvatures satisfy the SAME strict lower bound

    min(C_t,C_l)>=rho nu d² A_* abs(kappa)/(3 abs(lambda))-R_total.

Selecting a finite signed carrier with
`rho nu d² A_* abs(kappa)/(3 abs(lambda))>R_total`
establishes C_t>0 and C_l>0. Every item on the right was frozen before
this carrier selection. The equation is a finite analytic sufficient
threshold, not a numerical comparison to a target elastic coefficient.

## Route result and transferable scope

Established as stated: the selected stationary actual-EPS relative-angle
assembly admits geometric gradient-only cages whose COMPLETE action,
reaction and kinetic normal-form curvatures are positive at a finite
analytic carrier. All coefficients remain full stationary operator
functionals. No isolated-cell factorization or all-wave-number PDE is used.

The proof is a gradient-dominance theorem, not a claim that the earlier
relative kinetic convention has an independently nonzero physical
translation/spin contrast. The parent's later rotation-retaining base
construction can use this same theorem: its base S may retain affine
rotation moments, while the NEW cages still have all eleven zero moments.
The disjoint/isotropic response argument, fixed kinetic operator and
uniform mixed-jet bound are unchanged. That use requires its own fixed
base physical/current normal form and remains a parent joining task.
