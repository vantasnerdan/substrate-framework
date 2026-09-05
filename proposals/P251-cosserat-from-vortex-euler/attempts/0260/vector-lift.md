# Two exact cohomological Kelvin preparations on one fixed annulus

Let an oriented axisymmetric action-angle chart have coordinates
(I,theta1,theta2), theta2 toroidal, volume form J(I)dI dtheta1 dtheta2,
and u=Omega^a(I)partial_thetaa. Assume curl u=lambda u, lambda nonzero,
throughout a neighborhood of the compact source support. Fix a nonzero
toroidal Fourier index n and an open action interval with no resonance
m Omega1+n Omega2=0 for integer m. On a compact subinterval, the fixed-n
cohomological inverse exists uniformly as in0250: only finitely many m
can approach resonance, and large |m| denominators grow linearly. The
metric is axisymmetric, so the fixed-n subspace is preserved.

## 1. Both tangent directions are neutral

Every b=b^a(I)partial_thetaa commutes with u. With k=grad I, k.b=0 and
material derivatives obey D_t b=Du b, D_t k=-(Du)^T k. Therefore

    a_b=-k cross b/|k|^2,       k cross a_b=b              (1)

satisfies the full-pressure leading Euler amplitude equation already
derived in0250. No spectral inference is needed: b=partial_theta1 and
b=partial_theta2 are two independent exact commuting tangent fields.
This is an additional polarization, not a retraction of the first one b=u.

## 2. Exact finite-carrier lift, including shear

Put f=A(I)exp(i N I)exp(i ell.theta), ell=(m,n), nu=ell.Omega, and
alpha=(a_b f)^flat. If e_A=partial_A x and det(e_I,e_1,e_2)=J, the
identity e_1 cross e_2=J k gives the exact covariant components

    alpha_1=J b^2 f,       alpha_2=-J b^1 f.              (2)

alpha_I can contain the poloidal metric harmonics, but retains index n.
For w=curl(a_b f), the coordinate curl formula gives

    w^I=J^-1(partial_1 alpha_2-partial_2 alpha_1)
        =-i(ell.b)f,                                    (3)
    w^a=iN b^a f+O(1)                                   (4)

at every fixed spatial derivative order, with the corresponding polynomial
N losses. Here O(1) means in the undifferentiated carrier amplitude; it is
not an N-independent bound for arbitrary spatial derivatives.

Write D=Omega.partial_theta and use its uniform fixed-n inverse. Define

    xi^I=lambda^-1 D^-1 w^I,
    xi^a=lambda^-1 D^-1 w^a+D^-1(xi^I Omega^{a prime}).   (5)

The second term is required by the actual bracket

    [u,xi]^I=D xi^I,
    [u,xi]^a=D xi^a-xi^I Omega^{a prime}.                (6)

Equations(5)-(6) establish lambda[u,xi]=w exactly. Taking divergence of
the bracket gives D(div xi)=div w/lambda=0. The divergence retains the
fixed nonzero toroidal index and the same nonresonant inverse, hence
div xi=0. Compact action support is preserved by D^-1, and the source
extends smoothly by zero off its interior annulus.

Consequently

    curl(xi cross omega)=[omega,xi]=lambda[u,xi]=w,
    P(xi cross omega)=P(a_b f).                           (7)

The second equality uses the whole-space Helmholtz projection: the two
projected fields have the same curl and divergence and are square-integrable,
so their difference is zero. No hidden constant harmonic field is admitted
in that space. Both pre-projection vector fields have compact support.
This is exact isovortical Euler/Lin data on the actual fixed ring, even
though the normalized leading Kelvin symbol at k.omega=0 is singular.

For each fixed derivative order s, (5) gives finite bounds C_s N^(s+1).
The fixed frequency interval, chart, support and inverse constants are
selected first. These are source costs, not a bounded right inverse on the
space of all desired observation histories.

## 3. A normal column and a tangential column at the same frequency

Equation(3) yields exactly

    xi^I=-(ell.b)/(lambda nu) f,
    xi^a=N b^a f/(lambda nu)+O(1).                        (8)

For b=u, this reproduces xi^I=-f/lambda. For

    b_perp=-n partial_theta1+m partial_theta2,             (9)

ell.b_perp=0, so xi_perp^I=0 exactly. The two b directions are independent
because det(Omega,b_perp)=nu is nonzero.

For a stationary material tag chi(I) and a scalar laboratory weight F,
integration by parts gives the literal initial configuration moment

    delta M_F=-rho0 integral F chi'(I) xi^I dx.            (10)

Thus the b_perp preparation has zero initial centroid and covariance
variation for every such F, whereas b=u has the known nonzero normal row.
The tangential displacement of b_perp is not zero and can contribute to G
and full current. Its actual measured gain and subsequent normal history
remain the separate0250 computation. Formula(10) alone is not an all-time
zero-output theorem: full-pressure evolution can generate subleading
normal displacement, and its size must be compared to the same measured gain.

## 4. Distinct-polarization KKS pairing

This is an append-only analytic consequence of the lifted directions, not
a positive-action normalization claim. Let Xi_u(f) and Xi_perp(f) denote
the two exact preparations. Since Xi_perp has zero I component, the
coordinate triple product and (8) give

    omega dot [Xi_u(f) cross Xi_perp(conjugate(f))]
        = (J N/lambda)|f|^2+O(1).                       (11)

The sign follows from det(Omega,b_perp)=nu and Xi_u^I=-f/lambda.
The complex conjugate reverses both N and nu in the second source, leaving
N/nu unchanged. The O(1) amplitude is uniformly integrable on the fixed
support. Therefore the physical KKS bilinear form has leading term

    Omega_KKS(Xi_u(f),Xi_perp(conjugate(f)))
        =(rho0 N/lambda) integral J(I)^2 |A(I)|^2
                                   dI dtheta1 dtheta2 + O(1).       (12)

The two factors of J are distinct: one is the physical triple product in
coordinate components, the other is the physical volume element. For a
nonzero A and fixed sign lambda, the coefficient cannot vanish. In real
quadratures the same-phase cross-polarization entry carries half this
leading value. Terms with doubled toroidal index integrate to zero; this
does not require a delicate radial cancellation.

This nonzero cross-polarization pairing is compatible with0250's exact zero
KKS pairing between the two quadratures of its single span{grad I,u}
polarization. These are different pairs. It supplies an available symplectic
cross direction at sufficiently large finite carrier. The Jacobi energy,
all other Gram entries, positivity and physical action normalization are
still actual computations; none is inferred from nonvanishing of(12).

## Result and continuation

The exact vector cohomological lift, its shear/divergence conditions, both
normal-displacement rows and the nonzero leading distinct-polarization KKS
pairing are established as stated. The next use is0250's literal
same-frequency angle/G history matrix with full-pressure errors divided by
the actual fixed-tag gain, followed by the complete ambient hybrid and
Jacobi/current normalization. This does not transfer to a periodic
approximant without its invariant-action/cohomology construction, and does
not supply0255's compact stationary carrier or parent completion.

## 5. Failure-derived axisymmetric source extension

The fixed-n nonresonance condition was sufficient, not necessary. For n=0
and m nonzero, choose axisymmetric geometry and the same pure poloidal
harmonic f. The angular covector components in(2) have zero poloidal mean.
The coordinate curl components obey

    w^1=J^-1(partial_2 alpha_I-partial_I alpha_2),
    w^2=J^-1(partial_I alpha_1-partial_1 alpha_I).         (13)

Here partial_2 alpha_I=0; averaging partial_1 alpha_I around the poloidal
circle gives zero. Consequently every component of w has zero poloidal
mean, regardless of the poloidal metric harmonics in alpha_I. Equation(3)
also makes xi^I a pure nonzero poloidal harmonic, so its shear forcing in(5)
has zero mean. One can therefore use D^-1=(Omega1 partial_1)^-1 on the
actual mean-zero source, uniformly wherever Omega1 is separated from zero.
No irrational frequency ratio is required for this axisymmetric extension.
The resulting divergence has zero poloidal mean and lies in the kernel of
D, hence is zero as before. All other conclusions(7)-(12) continue to hold
with nu=m Omega1 and b_perp=m partial_theta2.

This provides exact axisymmetric field-changing sources. Toroidal symmetry
forces their transverse centroid and covariance-tilt components to vanish,
but permits an axial centroid and axial angular current. Those are available
symmetry sectors for the separate physical gain calculation. Their existence
does not establish the full hybrid pressure row or its normalized histories.

The candidate expansion is motivated by the remaining independent acoustic
output after the transverse angle/G repair. It preserves the fixed ring,
constant-curl region, action-supported smoothness and finite source costs;
it adds an actual mean-zero cohomology license instead of assuming that the
full n=0 cohomological operator has no kernel.
