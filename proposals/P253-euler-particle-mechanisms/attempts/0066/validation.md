# Analytic validation and claim-boundary checks

No production numerics were used in 0066. The exposing checks are exact
calculus, domain, and scaling checks against the physical whole-space Euler
operator.

## Fixed-carrier source and foliation checks

1. Cao Appendix A (A.8) gives a nonzero normal derivative on the actual
   scaled free boundary. Local elliptic regularity for the exact positive-part
   equation and `p>=6` gives `C^2` regularity at that fixed boundary. Normal
   Taylor expansion therefore yields `P=kappa_P d+O(d^2)` and
   `zeta=kappa_zeta d^p+O(d^(p+1))` with two-sided fixed-carrier comparison.
2. Proposition 3.2 is used only for `C^1` convergence. Applying Schauder to
   the difference of the exact rescaled equation and the radial limit upgrades
   this to `C^(2,alpha)` on a fixed ball for the sufficiently thin family.
   The radial center's nondegenerate Hessian, a positive gradient lower bound
   off a center ball, and Appendix A's single regular outer boundary then
   exclude secondary critical points and separatrices. Center linearization
   and the nonzero boundary speed give the positive finite endpoint limits of
   the period. This proves the fixed-carrier foliation used below; it is not
   attributed to Proposition 3.2 alone.

## Physical DA and KKS checks

In physical orthonormal cylindrical components, direct substitution of (15)
gives

    div xi=partial_r f+f/r+(i n/r)xi_theta+partial_z(i g)=0.

The definition `C_0 xi=-[xi,omega_0]` therefore produces a compact smooth DA
test seed. Pairing its real and imaginary polarizations in the physical KKS
form gives exactly

    Omega(q_test,conjugate(q_test))
      =4*pi*i*rho_0*integral r^2*zeta*f*g dr dz,

which is nonzero by the frozen choice of `f,g`. This validates a physical
symplectic seed but not an eigenmode.

## Exposing constrained Weyl sequence

On one regular cell choose volume action--angle coordinates. For `n!=0`,

    a_(m,n)=partial_beta-(m/n)partial_theta

is tangent to `I=constant`, kills the transport shear, and makes (26a)
exactly divergence-free. Since `zeta=zeta(I)`, the exact bracket is

    C_0 xi=i*n*zeta*xi-xi^I*zeta'(I)partial_theta.

Thus on this packet `C_0 xi=i*n*zeta*xi`: `C_0` is order zero in the large
radial frequency and preserves the principal polarization. For a packet of
point amplitude `A_N`, width `epsilon`, and radial frequency `N` with
`N*epsilon->infinity`, equivalence of coordinate and physical norms gives

    ||C_0 xi||_H3=Theta(A_N*N^3*sqrt(epsilon)).

Taking `A_N` equal to the reciprocal scale and applying the product rule to
`(omega(I)-omega(I_0))q_N` gives

    ||(T_n+i*m*omega(I_0))q_N||_H3
      <= C*(epsilon+N^-1).

There is no divergence correction. Axisymmetric physical rows vanish by
toroidal orthogonality. Any additional continuous same-character finite rows
are `o(1)` on the weakly null packet; subtraction of fixed smooth DA profiles
with an invertible row matrix changes both its norm and residual by `o(1)`.
Fixed-`n` Hodge compactness then gives `K_n q_N->0`. Hence the regular-cell
transport-band inclusion in (27) is established in the constrained DA
Fredholm essential spectrum. The converse direct-integral theorem and all
boundary/separatrix spectral contributions remain open, so (28) is only a
candidate gap.

## Conditional-transfer tripwires

The following predicates were checked textually in the claim-bearing files:

- equations (34)--(39) are explicitly targets or conditional consequences;
- no scalar coercivity estimate is promoted to the missing full vector
  Piola/curl/div/Leray graph estimate at `delta*n_delta=O(1)`;
- two limiting column modes are not called Cao modes before entire-column
  common-domain isolation and contour-relative convergence;
- the weighted topology is a tangent quotient, not a descended nonlinear
  orbit chart;
- the charged route requires every streamline compatibility mean, the
  nonlinear steady free-boundary map, and its translation Lyapunov--Schmidt
  step; no charged branch or stability signature is claimed.

These checks expose the exact positive results while keeping the converse
essential-spectrum theorem, full gap, Cao eigenmodes, nonlinear branch,
stability, and P2/LP2 outside the 0066 claim.
