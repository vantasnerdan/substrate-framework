# Actual energy and the velocity Jordan equation

## 1. All contrasts of the cubic-first-integral preparation

Put E=A²+B² and Delta=B²-A². Use the particular range solution phi_* in
`unequal-stationary-corrector.md`, and include all three indicated stationary
first-shell/passive freedoms:

    phi=d(phi_*+eta alpha+xi psi),
    r=d(-alpha/2+zeta psi), z=(r+Hphi,-phi_Z,phi_Y).

The complete initial material rate is q_D-aD+z. Its actual Hamiltonian
coefficient H_DD=rho k² h2/2 is

    h2=E(eta²+xi²)/5+E zeta²/10
       +2Delta eta xi/5+Delta eta zeta/5+E xi zeta/5
       -E eta/5-Delta xi/5-Delta zeta/5
       +E[25E²/(64Delta²)-371/960].                  (1)

This is a calculation of the complete Euler energy, not energy assigned
to an observed frequency. Orthogonality of Fourier shells makes the range
solution independent of the first-shell terms. Its actual velocity norm is

    ||z_*||²=5(A^6+B^6)/(16Delta²)
               +15(A²B^4+A^4B²)/(2Delta²)
             =5E(25E²-21Delta²)/(64Delta²).

Its whole-law weight is E[d²]=1/5. On the first shell,
<alpha²>=<psi²>=E/2, <alpha psi>=Delta/2; the two derivative rows have
the same covariance. Contraction of the actual q_D-aD gives the remaining
linear terms and the constant -7E/120. These facts give(1) directly.
The verifier independently constructs all Fourier velocity rows at A=1,B=2
before comparing the general formula.

The observed restoring coefficient is

    a=(E eta+Delta xi)/10-E/6.                       (2)

The Hessian of h2-a is positive definite when 0<|Delta|<E. Its leading
principal minors are 2E/5, 4(E²-Delta²)/25 and
2E(E²-Delta²)/125. The unique minimum occurs at

    eta=3/4, xi=-Delta/(4E), zeta=Delta/(2E).

Writing r=|Delta|/E gives

    min(h2-a)=E(375/r²-319-12r²)/960
             =11E/240+E(1-r²)(375+12r²)/(960r²)>0.   (3)

Thus this COMPLETE specified cubic-first-integral/three-control family
does not match the physical D energy to its observed acceleration, even
where that acceleration is restoring. This is not a no-go for other
first integrals, initial configuration returns, coupled optical energy,
or the parent objective. It is the exact defect a repair has to cancel.

## 2. The actual common-velocity equation

For common V perpendicular to kappa, the physical zero-K velocity column
is w_V=V+tT_V, T_V=-(Du)V. The chosen divergence-compatible translation
Bloch lift is T_V+ik q_V, q_V=kappa cross T_V/lambda. It is not the
literal Kelvin lift P_K(V cross omega). Since div q_V=-kappa.T_V,
the pressure projection derivative gives

    L1 V=-grad Delta^(-1)(kappa.T_V)-P[(kappa.u)V],
    b_V=L1 V-q_V=-P[(kappa.u)V+q_V].                 (4)

Here L_K V=T_V+ik L1 V+O(k²). This includes the full Leray derivative;
using the Kelvin-lift expression instead is the representation mismatch
recorded in `velocity-lift-diagnosis.md`.
Write w_V=V+tT_V+ik(tq_V+z_V)+O(k²). Then

    z_V,t=L z_V+t F_V+b_V,
    b_V=-P[(kappa.u)V+q_V].                          (5)

For a stationary D corrector Lz_D+F_V=0, the full-field Jordan ansatz
z_V=t z_D+r_V requires the ACTUAL equation

    L r_V=z_D-b_V.                                  (6)

It is not enough to differentiate an averaged frequency. A nonzero r_V
is an O(k) microscopic velocity preparation, so it changes gradient
mass and cross phase/energy rows, but not the leading rho kinetic mass.

For the exposing planar input kappa=(0,1,1), V=(0,1,-1), d=2,
the actual b_V has planar stream function d alpha and axial component
d alpha/2. If r_V has planar stream function phi_R, equation(6) entails

    H(phi_D-d alpha)=-T mathcalB phi_R.

Consequently every smooth first integral g(psi) gives the necessary moment

    <g(psi) H(phi_D-d alpha)>=0.                     (7)

This is a necessary test for the stated full-field Jordan preparation,
not a new requirement of per-orientation stationarity for whole-law
physical moment closure. The latter may still use correlated nonstationary
returns whose unobserved components cancel only after law averaging.

The psi, psi³ and psi^5 tests below expose whether the finite polynomial
D construction actually admits this particular velocity partner. Actual
pressure, physical initial means, initial phase, and gradient energy remain
part of the response even if a Jordan construction passes these moments.

Solving the psi and psi³ moments fixes eta-1 and xi. After those controls,
the psi^5 residue is

    5(A^8-8A^6B²-80A^4B^4-8A²B^6+B^8)/(128(A²-B²)),

and the psi^7 residue is

    21(A²+B²)(2A^8-17A^6B²-208A^4B^4-17A²B^6+2B^8)
        /(512(A²-B²)).

Their numerator polynomials have no common root in B/A: the exact
polynomial gcd is constant. Thus the specified cubic stationary family
has no full-field planar Jordan partner at any nonzero unequal contrast.
This verdict does not impose that preparation on whole-law moment closure.

## 3. Why planar high harmonics alone do not repair the positive action

More first-integral harmonics can lower the range energy R=||z_range||².
For example at the small-A/B endpoint, odd Chebyshev coefficients satisfy
c1=-B and sum_(n>=3 odd) n c_n=2B. Distributing that constraint among N
harmonics lowers the exact minimum range energy as O(B²/N); at fixed N
the finite polynomial moment equations continue to small positive A/B.

That reduction is not a positive-action conclusion. Keep the same axial
return r=-d alpha/2+d zeta psi and put t=(E eta+Delta xi)/E,
r0=Delta/E. After minimizing the remaining first-shell controls, the
mismatch at the restoring boundary t=5/3 is

    min(h2-a)|_(a=0)
       = E(27r0²+59)/(360(1+r0²))+R/5>0.

The derivative with respect to t is positive for t>=5/3, so the entire
positive restoring branch retains this defect. The next registered
construction changes the ACTUAL axial first integral g(psi), not merely
the number of planar harmonics. Initial configuration returns remain a
distinct physical phase alternative.
