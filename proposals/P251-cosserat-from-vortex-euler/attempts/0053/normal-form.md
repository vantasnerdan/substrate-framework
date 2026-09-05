# Exact isotropic average and second-gradient normal form

Let a be a unit cell axis, C a positive real symmetric gradient matrix, and
R uniformly distributed in SO(3). Rotate the entire cell, not a alone.
For T=tr C and L=a.C.a, isotropy and index-pair symmetry give

    <(Ra)_i (Ra)_k (R C R^T)_jl>
      = A delta_ik delta_jl + B(delta_ij delta_kl+delta_il delta_kj).

Its two invariant contractions are 9A+6B=T and 3A+12B=L. Therefore
A=(2T-L)/15, B=(3L-T)/30. For density nu and G_ij=partial_j Phi_i,

    W_curv=nu/2 [A ||G||²+B((tr G)²+tr(G²))].

In W=c_tr(tr G)²+c_s||sym G||²+c_a||skew G||² the coefficients are
nu*B/2, nu*(A+B)/2, nu*(A-B)/2. Since 0<L<T,
c_s=nu*(3T+L)/60>0, c_a=nu*(T-L)/12>0,
3c_tr+c_s=nu*L/6>0. No false separate c_tr>0 premise is needed.
For periodic fields or retained boundary flux, the bulk transverse and
longitudinal coefficients are nu*A and nu*(A+2B), both positive.

The uniform cell action J_Psi*(a.Phi_dot)²/2
-K_Psi*(a.(Phi-beta))²/2 averages to j|Phi_dot|²/2
-alpha|curl U-2Phi|²/2, where beta=curl U/2,
j=nu*J_Psi/3 and alpha=nu*K_Psi/12. This coefficient map takes the
complete same-orbit Schur quantities, not the retired tension formula.
The residual cage inertia gives nu*J_beta*|curl U_dot|²/24 and remains
until the following kinetic normal form is taken.

## Retaining the kinetic-gradient terms

Assume the complete time-even reflection-paired isotropic reduced action has,
in a transverse curl helicity h=+/-1 (curl has Fourier eigenvalue h*k),

    M = [[rho+m_U k², b h k], [b h k, j+m_Phi k²]],
    K = [[A_U k², -2 alpha h k], [-2 alpha h k, 4 alpha+C k²]].

These are jet coefficients from the microscopic action. Neither positivity
nor origin of a microscopic coefficient follows merely from this display.
Let d=m_Phi-b²/rho. The near-identity map from normal fields to physical
fields is

    U_phys=(1-m_U k²/(2rho))*U - (b/rho)*h*k*Phi,
    Phi_phys=(1-d k²/(2j))*Phi.

Pull back BOTH matrices by this same map. Coefficientwise through k²,

    M_normal=diag(rho,j),
    K_normal=[[A_U k², -2alpha h k],
              [-2alpha h k, 4alpha+C_eff k²]],
    C_eff=C+4alpha*b/rho-4alpha*d/j.

The mixed inertia and gradient mass are not zeroed by assertion: the
potential correction is their observable effect. The transformation uses
curl and second derivatives, so it preserves uniform frame rotations and
translations. It changes nonuniform physical displacement/spin and must
be retained when mapping centroid or boundary-current observables. It is
invertible as a formal low-gradient series. Positivity and an all-k inverse
are not asserted outside that expansion. In the longitudinal spin sector
there is no curl mixing, giving C_L,eff=C_L-4alpha*m_Phi,L/j.

This is an exact normal-form statement for the original slow-affine
coarse-graining target, not a new claim that unrestricted Euler evolution
is exactly a finite-order differential model. The zero-density-of-structures
branch is evaluated in the unreduced action before dividing by j. The
microscopic cell/pressure/translation construction remains an independently
required input; this algebra does not supply it.
