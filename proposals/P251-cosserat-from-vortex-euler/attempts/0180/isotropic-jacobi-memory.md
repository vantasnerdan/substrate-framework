# An all-time identity for the actual whole-law stress-memory row

This is a response identity on the bare Kelvin-D preparation, not a claim
that an isolated mean must close to solve the coupled parent. It sharpens
the finite-moment test by exposing precisely what cancels after isotropy.

For each symmetric unit strain pair m,l define the ACTUAL cell data

    S^{ml}=P(u_m e_l+u_l e_m),
    F^{ml}=P(p_m e_l+p_l e_m), p_m=partial_m p,
    K_J=P(A²+Hess p),
    chi_tt^{ml}+2PA chi_t^{ml}+K_J chi^{ml}=F^{ml},
    chi^{ml}(0)=0, chi_t^{ml}(0)=-S^{ml}.              (1)

These fields enter the common-direction transverse contraction0175; they
are not separately declared realizations with different macro inputs.
Let E=<|u|²>. The exact full-field isotropic displacement stress is

    R_D,iso=E/3+Z/5-Y/15,
    Z=sum_ml <u_m chi_t,l^{ml}+p_m chi_l^{ml}>,
    Y=sum_il <u_i chi_t,i^{ll}+p_i chi_i^{ll}>.          (2)

The trace cell is exactly sum_l chi^{ll}=-2t u. Indeed sum_l S^{ll}=2u,
sum_l F^{ll}=2P grad p=0, PAu=0 and K_J u=0 by stationary Euler.
Consequently Y=-2E and

    R_D,iso=7E/15+(1/10)sum_ml
                         [<S^{ml},chi_t^{ml}>+<F^{ml},chi^{ml}>]. (3)

For every nonzero Fourier wave, the solenoidal polarization gives
sum_ml ||P(u_m e_l+u_l e_m)||²=6|u|² modewise. Hence
sum_ml ||S^{ml}||²=6E and R_D,iso(0)=-2E/15.

## Exact forced energy, with its real sign

K_J is self-adjoint on the solenoidal cell space, while PA is skew
adjoint. For each actual solution(1), the conserved forced Jacobi energy is

    (1/2)||chi_t||²+(1/2)<chi,K_J chi>-<F,chi>
                                                   =(1/2)||S||². (4)

It is an inherited material quadratic form; it is NOT assumed positive.
Using(4) to complete the square in(3) gives the all-time identity

    R_D,iso(t)=-2E/15+(1/20)sum_ml {
                  ||chi_t^{ml}+S^{ml}||²
                      +<chi^{ml},K_J chi^{ml}>}.       (5)

This is a genuine whole-law output relation. It does not require each
orientation's stress to be constant, nor a stationary chi. It keeps all
pressure and circulation-return modes through their actual operator.
The material K_J includes A², and can be indefinite even when a selected
optical KKS plane has positive mass. That optical sign does not make(5)
a sum of nonnegative terms.

For the exact one-wave finite space, K_J vanishes and chi_t=-S for the
displacement column, so(5) gives the known constant response. For0175's
actual equal two-wave field, the fourth derivative of the brace sum at
zero is -8/5, since R_D''''(0)=-2/25. This is a concrete negative internal
Jacobi contribution, not a failure of the square-completion identity.

## Construction meaning

The bare isotropic row is constant exactly when the displayed aggregate
quadratic form vanishes along its actual forced orbit. That can involve
indefinite cancellations; setting every internal field stationary is an
unjustified extra condition. Corrected preparations change both the
constant in(4) and the initial/output slip in(3), which must be recalculated.
The physical tagged angle/spin also has its own moving row equations.

Thus(5) is a new exposing interface for0179 and a candidate finite physical
closure, not a proof of generic autonomy or a no-go over corrected data.
It explains why the actual first several zero time derivatives do not
settle closure, while preserving the strongest exact stationary-memory
and finite one-wave constructions already obtained here.

## Corrected initial cells: the temporal identity still has a clean form

Keep the same physical strain forcing and symmetric input tensor, but let
chi(0)=xi and chi_t(0)=v be any actual smooth solenoidal corrected data.
For each cell the constant in(4) becomes
E_ml=(||v||²+<xi,K_J xi>)/2-<F,xi>. The trace forcing remains zero,
and Y=<u,sum_l chi_t^{ll}> because each pressure integral against a
solenoidal chi is zero. Since K_J u=0 and PAu=0, Y_t=0 for ANY such
corrected trace data, not only the bare trace -2tu.

Completing the same square in(3) therefore gives

    R_D,iso(t)-R_D,iso(0)
       =(1/20)sum_ml [Q_ml(t)-Q_ml(0)],
    Q_ml=||chi_t^{ml}+S^{ml}||²+<chi^{ml},K_J chi^{ml}>. (6)

The initial stress coefficient changes with the actual data, but its
constancy is exactly the constancy of this aggregate quadratic observable.
This extends the whole-law interface to0179's corrected preparation without
assuming stationary chi or positive K_J. If a corrected preparation changes
the macroscopic observation or forcing too, those changed rows must first
be inserted; (6) concerns the stated physical stress and strain forcing.
