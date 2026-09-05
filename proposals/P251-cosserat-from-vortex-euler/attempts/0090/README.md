# 0090 — compact positive material stiffness on the same EPS field

Owner `/root`. Active parent P251; this is a construction of the actual
material Jacobi form, not an assertion about the earlier coadjoint Hessian.
The original conditional slow-affine scope, density and Euler dynamics
remain unchanged. No empirical comparator or numerical sign threshold.

Frozen route: use a strictly positive exact periodic divergence-free
Rayleigh witness, localize its vector potential with a finite error bound,
then use the decaying EPS seed plus an arbitrarily small same-eigenvalue
periodic field. Place the material cage sufficiently far from the seed
to retain its sign without disturbing the knotted core. The alternative
single-plane Beltrami material cage is nonpositive by0084; that failure
motivated the multiwave pressure mechanism in0088/0089.

## 1. Exact periodic input and the proposition it supplies

0089 verifies, by both the pressure-Hessian and curl forms, a real smooth
divergence-free mean-zero trigonometric displacement X for

    b=(-sin(2y),cos(2x),cos(2y)-sin(2x)), curl b=2b,
    p_b=-|b|^2/2.

In normalized periodic volume measure,

    kappa=average[X.Hess(p_b).X-|b.grad X|^2]
          =5008301/1250000000 > 0,
    m=average |X|^2=14496029/125000000 > 0.

The finite rational coefficients and exact rounding repair belong to0089.
This input proves a positive material test direction, not spectral stability,
Euler eigenmode closure, or an electromagnetic assumption about the fluid.

## 2. Compact localization with an explicit analytic error hierarchy

Since X is periodic, divergence free and mean zero, its finite Fourier
polynomial vector potential A is defined by

    A_k=i k cross X_k/|k|^2,  k!=0,  A_0=0,
    curl A=X.

Fix a smooth real cutoff chi supported in B_2, equal to one on B_1. For
R>=1 define chi_R(x)=chi(x/R) and

    X_R=curl(chi_R A)=chi_R X+R^-1(grad chi)(x/R) cross A.

This is an actual smooth compact divergence-free material displacement.
All coefficients and derivatives of b,X,A,p_b have explicit finite Fourier
supremum bounds (sum absolute coefficient norms times powers of |k|).
Let X0,X1,A0,A1,b0,H0 bound |X|,|grad X|,|A|,|grad A|,|b|,|Hess p_b|.
Let c0,c1,c2 bound |chi|,|grad chi|,|Hess chi|, respectively, using
operator norms consistently. Put

    E0=c1 A0,
    E1=c1 X0+c2 A0+c1 A1.

Then |X_R-chi_R X|<=E0/R and
|grad X_R-chi_R grad X|<=(E1+c1 X0)/R for R>=1; the deliberately
redundant last term is a harmless conservative bound on differentiating
chi_R X. With D1=E1+c1 X0, the direct Jacobi integrand differs from
chi_R^2 times the periodic integrand by at most C_loc/R, where

    C_loc=H0(2 c0 X0 E0+E0^2)
          +b0^2(2 c0 X1 D1+D1^2).

All support is in B_(2R), so this contribution is bounded by
Vol(B_2) C_loc R^2. There is no unmeasured numerical cutoff error.

Let f=X.Hess(p_b).X-|b.grad X|^2. Its mean is kappa. The zero-mean finite
Fourier polynomial f-kappa has periodic vector primitive B with
div B=f-kappa, for example B_k=-i k f_k/|k|^2. Hence

    |integral chi_R^2(f-kappa)|
      <=2 c0 c1 ||B||_infinity Vol(B_2) R^2.

Define C=Vol(B_2)[C_loc+2c0 c1||B||_infinity] and
J_chi=integral chi^2>0. The COMPLETE compact material stiffness satisfies

    K_b(X_R)/rho >= kappa J_chi R^3-C R^2.

Any finite R>max(1,2C/(kappa J_chi)) gives the strict lower bound
K_b(X_R)/rho>kappa J_chi R^3/2. The cutoff is a declared material ensemble
geometry, not natural localization of an unrestricted eigenmode. This is
an exact compact test field, not a large-box numerical conclusion.

## 3. Preserve the actual EPS core and the compact positive cage

Take an actual decaying EPS Beltrami seed u_E with the desired structurally
stable invariant knotted tube, and rescale its eigenvalue to2 if necessary.
The R^3 existence/decay theorem and scaling are the archived inputs in
0057/0076. The entire smooth field

    u=u_E+epsilon b,  curl u=2u,
    p=-rho |u|^2/2

is an exact stationary Euler field. Choose epsilon>0 first, small enough
in the fixed core C^m norm to preserve the EPS invariant tube. Choose the
finite R above. Translate X_R to a centre a sufficiently far from the core.
Periodicity permits taking a along lattice translations of b, so its
periodic stiffness is unchanged. EPS derivatives through second order tend
to zero uniformly on a+B_(2R) as |a| tends to infinity.

For a fixed compact X_R, the Jacobi stiffness is continuous in the C^2
background norm: expand Hess(-rho|u|^2/2) and |u.grad X_R|^2, retaining
both mixed terms. The error is bounded by

    rho [C_p(delta,epsilon)||X_R||_2^2
          +(2epsilon b0 delta+delta^2)||grad X_R||_2^2],

where delta bounds the C^2 norm of u_E on the cage and C_p is the explicit
product-rule bound for the pressure-Hessian difference, tending to zero
with delta. For example with all derivative tensor norms bounded by delta,
one may take C_p=epsilon delta(2b0+4b1+2b2)+4delta^2, with b1,b2 the
corresponding periodic derivative bounds. Increase the finite separation
until this is less than half of epsilon^2 K_b(X_R). Then

    K_u(X_R translated)>epsilon^2 K_b(X_R)/2>0.

Thus the compact positive material direction and the stationary EPS core
belong to ONE smooth Euler field, with no glued interfaces, singular core,
external wall, or borrowed independent energy. Their finite separation is
part of the declared coherence geometry.

## 4. Attach a physical material core jet without losing positivity

Let Y be a fixed smooth compact divergence-free material displacement with
the required rotational core jet, supported away from the cage. Such curl
extensions are constructed in0084. Since the material Jacobi functional
is local, disjoint supports give exactly

    K_u(Y+t X_R)=K_u(Y)+t^2 K_u(X_R),
    M(Y+t X_R)=M(Y)+t^2 M(X_R).

For a finite t with t^2 K_u(X_R)>max(0,-K_u(Y)), this gives positive
stiffness and positive actual material mass while preserving Y's core jet.
The cage amplitude is a declared kinematic ensemble input, selected by an
analytic sign inequality, not fitted to a desired modulus. Its spin,
mean displacement, boundary motion, gyro terms and gradient masses must
all be carried into the joint reduction. No equality with the old
coadjoint inertia, core angular momentum alone, or a postulated rigid
rotor is made.

Route verdict: established for compact positive MATERIAL Jacobi directions
with a prescribed core jet on an actual EPS-containing stationary field.
The remaining parent construction is the exact varied material/GLM or
fixed-Kelvin reduction and its physical collective-spin/translation map;
0087 and the continuation of0084 execute that join. Positive Rayleigh
stiffness alone does not certify that larger action or the parent goal.
