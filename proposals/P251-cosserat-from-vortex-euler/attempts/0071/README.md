# 0071 — isotropic ergodic Euler samples by finite Gaussian superposition

This is an optional stronger assembly representation, generated from the
explicit distinction in 0057/0060 between ensemble isotropy and spatial
ergodicity. It does not invalidate the reviewed Haar-mixture construction.
The parent remains the declared-ensemble Euler-to-Cosserat objective. The
new route supplies a genuine isotropic ergodic Gaussian law with the same
normalization, full local support and exact stationary Euler samples. No
micropolar coefficient or datum selects the law.

## Finite construction from the archived primary source

Use twelve independent copies u_l of the exact Gaussian Beltrami field of
Enciso--Peralta-Salas--Romaniega, arXiv:2006.15033, Props. 3.2--3.8, archived
and applicability-reviewed in 0057/0060. Let e_l be the twelve normalized
vertices of an icosahedron, and choose any proper rotation R_l taking e1
to e_l. Define

    u_iso(x)=12^(-1/2) sum_l R_l u_l(R_l^T x).

All copies have the SAME curl eigenvalue. Thus every sample is almost surely
smooth and curl u_iso=lambda u_iso. Superposition is legitimate here:
although Euler is nonlinear, every such Beltrami sum is an exact stationary
Euler solution with p=-rho |u_iso|^2/2. No independently glued fluid domains
or interfaces occur. Independence makes the sum Gaussian and makes its
covariance the average of the twelve rotated covariances. Point covariance
is still I; the physical rescaling gives U_*^2 I and e0=3rho U_*^2/2.

The source spectral vector p(n), Eq. (2.4), is a scalar linear polynomial
in n1 times (n1^2-1,n1*n2-i*n3,n1*n3+i*n2). On the unit sphere it spans
the positive-helicity line, so

    p(n)p(n)*=w(n1) H_+(n),
    H_+(n)=(I-n n^T+i[n cross])/2,

where w=|p|^2 is an even polynomial of degree at most four in n1. Its
precise source normalization is not guessed from PDF text: the independently
stated point covariance fixes the final normalization below. The two exact
icosahedral identities are

    (1/12)sum_l(e_l.n)^2=|n|^2/3,
    (1/12)sum_l(e_l.n)^4=|n|^4/5.

Hence the averaged w(e_l.n) is constant on the unit sphere. Since H_+
is rotation-covariant, the entire averaged spectral covariance is isotropic,
not merely its point covariance. Its scalar weight is 3/(4pi), because
tr H_+=1 and integral tr covariance=3. A centered Gaussian law is determined
by its covariance; therefore this finite sum has exact SO(3)-invariant law.

This is not the global Haar mixture of twelve field laws. Independent sums
average covariances BEFORE fixing a Gaussian law; a random choice of one
rotated sample averages laws and need not have this property.

## Ergodicity and the same full-support EPS reconstruction

The resulting spectral measure is atomless with smooth sphere density.
The argument of source Prop. 3.7 therefore applies without a new spectral
hypothesis: the Gaussian translation action is ergodic. More explicitly,
the covariance entries are finite derivatives of the sphere transform
4pi sin(r)/r and tend to zero as r tends to infinity. Finite Gaussian
coordinate blocks at increasingly separated positions become independent;
cylinder approximation extends this mixing identity to bounded measurable
observables. Consequently volume averages of integrable translation-covariant
observables agree almost surely with their ensemble expectations.

Full C^m compact-open support is inherited directly, without assuming that
the source theorem itself named this new law. To approximate any Beltrami
target v on a compact ball, prescribe the first transformed summand within
epsilon/(2sqrt(12)) of sqrt(12)*v and each remaining transformed summand
within epsilon/(22sqrt(12)) of zero. Every one of these events has positive
probability by source Prop. 3.8, and independence gives positive probability
for their intersection. The triangle inequality then puts u_iso in the
desired epsilon neighborhood. The same good EPS prototype plus small ABC
perturbation and finite Gram bounds of 0057/0059 is therefore available.

For an isotropic spatial selection, use the 0059 Poisson candidate rule
rather than a globally rotated cubic grid: intensity tau>0, independent
Haar orientation mark at each point, retain a candidate only when no other
candidate lies within 2R and the field on its R-ball lies in the rotated
good-patch neighborhood. Its intensity is

    nu=tau exp[-tau volume(B_2R)] p_good >0.

Poisson independence, isotropy of the field law and compact-open full support
supply the displayed factors. The joint Gaussian/Poisson system is mixing
and its translation-covariant thinning and generator reconstruction are
factors, hence ergodic. The marks are then frozen and materially transported,
not reselected after a deformation. All accepted patches belong to the ONE
sampled Euler field. The full nonlocal reaction inverse remains unchanged
in meaning; spatial ergodicity does not turn it into isolated-cell inverses.

Reflection still changes helicity, and time reversal still changes KKS
signs. The declared independent reaction variations for paired opposite
samples are retained before averaging actions. This new law makes each
fixed-helicity spatial assembly isotropic/ergodic; it does not claim a
single constant-lambda field contains both helicities, nor that orientation
ergodicity by itself implies zero angular current.

## Evidence scope

The analytic finite-sum proof and source support/ergodicity argument establish
this stronger stationary ensemble construction. The exact verifier derives
the icosahedral moments and helicity projector from coordinates; it does not
claim a finite random sample proves ergodicity or an EPS existence theorem.
The 0057/0069 full-action estimates apply to this law because their only
statistical inputs are the proved local support, finite moments, covariance
normalization and translation stationarity. Physical rotation and gradient
coupling remain the separate 0070/0065 constructions.

The initial exact verifier reached four passing checks and then hit an
implementation-domain error: SymPy inferred integer coefficients for a
Groebner reduction containing rational projector entries. The first log is
preserved; explicitly using QQ repairs that domain without changing the
polynomial identities or their acceptance criteria.
