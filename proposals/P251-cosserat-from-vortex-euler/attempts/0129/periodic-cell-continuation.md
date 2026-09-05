# Executed continuation: the true mean cell equation and a complete Euler branch

The Bessel crossing's physical-mass mismatch motivates the registered
periodic Jordan-sector alternative. This file executes it on a simple
stationary Beltrami field. It finds a genuine negative acoustic
coefficient, not the desired positive one, and thereby distinguishes a
real dynamical failure of this candidate from the former misidentified
carrier. It is not a no-go for other stationary ensembles or the parent.

## 1. The actual first cell equation

Let a fixed zero-mean periodic stationary Beltrami field satisfy
curl u=lambda u. Put E=|u|^2/2, t_U=-partial_U u, and choose a propagation
direction n with n·U=0. Let B(W)=u·grad W+W·grad u. For a Bloch wave
number K=epsilon*n, seek a mean translation branch with omega=epsilon*c:

    v_K=t_U+i epsilon W-i epsilon c U+O(epsilon^2),
    pi_K=partial_U E+i epsilon Pi+O(epsilon^2).         (1)

The actual Bloch divergence condition yields div W=-n·t_U. In Euler,
the order-epsilon time derivative -i c t_U is exactly canceled by the
convective term of the Galilean mean -i c U, since B(U)=-t_U. What
remains is the independently testable forced cell equation

    B(W)+grad Pi=-[(u·n)t_U+n partial_U E],
    div W=-n·t_U,   <W>=0.                             (2)

The complete mean momentum equation then gives

    c^2 U=-P_n <u tensor W+W tensor u> n.               (3)

Here the actual mean velocity is -i omega U at leading order, so its
mass is the actual total density rho. Formula (3) derives the candidate
stiffness from the actual Euler response; it does not insert one. A
solution of (2) and the remaining spectral solvability are substantive
parts of the construction, not consequences of a static affine action.

For a Beltrami field the forcing in (2) equals u cross (n cross t_U).
The vector identity B(W)=grad(u·W)+u cross (lambda W-curl W) shows that
the sufficient ansatz

    (curl-lambda)W=-n cross t_U                         (4)

would solve it with Pi=-u·W. However (4) is resonant on the source's curl
shell. At a Fourier vector p with |p|=|lambda| and source coefficient u_p,
the projection of its right side onto the same-helicity nullspace is

    (n·p)(U·p) u_p/lambda.                             (5)

This need not vanish. Blind inversion of curl-lambda would delete the
actual response. The general pressure/reaction equation (2), rather
than that shortcut, is the next useful object.

## 2. A cell solution with a fully derived sign

Take the exact stationary positive-curl field

    u(z)=V (cos(lambda z),-sin(lambda z),0),  V>0,
    n=e_x,  U=e_z,  t_U=-u'(z).

Equation (2) has the explicit solution

    W=V cos(lambda z) e_z,   Pi=0.                     (6)

Indeed div W=-lambda V sin(lambda z)=-n·t_U,
B(W)=W_z u'(z), and partial_U E=0. The mean is zero. Hence

    c^2=-V^2/2.                                       (7)

This z-polarized coefficient is unchanged by adding any zero-mean
periodic homogeneous corrector: its divergence-free difference Z has
zero horizontal-average Z_z, because that average is constant in z and
has zero total mean. Therefore <Z_z u_x>=0. No favorable kernel gauge
can change (7). The physical mass is positive; the dynamical stiffness
of this particular stationary shear is negative.

## 3. This sign is realized by genuine Euler modes, not only a formal cell jet

Write Ux(z)=V cos(lambda z), Uy(z)=-V sin(lambda z), and let c be a
complex phase speed away from the range of Ux. Seek an actual perturbation
proportional to exp(i K x-i K c t). With a periodic scalar f(z), set

    v_x=-[(Ux-c) f]',
    v_y=-Uy' f,
    v_z=i K (Ux-c)f,
    pi=(Ux-c)^2 f'.                                   (8)

Direct substitution into ALL THREE Euler components and divergence
reduces them to the single equation

    [(Ux-c)^2 f']'=K^2 (Ux-c)^2 f.                     (9)

The exact material displacement is

    xi=(i f'/K,0,f) exp(iKx-iKct),                     (10)

which is divergence free. Its Lin velocity is precisely (8). Moreover
v-xi cross curl u is the gradient of i(Ux-c)f'/K; the last component of
this identity is exactly (9). Thus the mode is genuinely fixed-Kelvin,
not a forced velocity ansatz.

Here is an analytic existence proof near K=0. Put a_c(z)=(Ux-c)^2,
normalize <f>=1, and split (9) into its mean-zero equation and scalar
mean. For c near c0=+iV/sqrt(2) or -iV/sqrt(2), a_c never vanishes. The
operator L_c f=(a_c f')' is invertible between mean-zero periodic spaces
provided <1/a_c>!=0. This can be checked without a solver: integrate
a_c f'=G+C for a periodic primitive G of a zero-mean right side; periodicity
fixes C=-<G/a_c>/<1/a_c>, and mean zero fixes the remaining constant.
The elementary unit-circle residue integral gives

    <1/(Ux-c)^2>=c/(c^2-V^2)^(3/2),                    (11)

with the analytic square-root branch selected by the usual large-c
continuation. It is nonzero at either c0. This supplies the actual
bounded analytic inverse near those points.

Solve the mean-zero part as
f=1+K^2 L_c^(-1) Q(a_c f). A Neumann series gives f(c,K^2)=1+O(K^2)
in every fixed periodic smooth norm, analytically in c and K^2. The
remaining scalar equation is <a_c f>=0. At K=0 it is

    <a_c>=V^2/2+c^2=0,

with derivative 2c0!=0. The analytic implicit function theorem gives

    c(K^2)=c0+O(K^2),
    omega(K)=+/- i V K/sqrt(2)+O(K^3).                 (12)

These are exact smooth linearized Euler modes for small nonzero K, not
just the necessary mean condition (7). Rational K is realized on a finite
periodic supercell. Real perturbations combine the conjugate modes.
The growing and decaying branches have an actual small-K rate and can
be realized over any fixed finite interval by sufficiently small smooth
Euler disturbances. No global-time claim is needed.

At K->0, (10) tends to the genuine material translation e_z, while
<v_z>=-i K c+O(K^3). The mode therefore carries the appropriate actual
mean momentum and the negative coefficient (7). This is precisely the
physical mean/Jordan sector missing from the neutral Bessel pattern.

## Route conclusion

The periodic single-helical-wave acoustic construction is **refuted as
a positive acoustic candidate**, with its negative stiffness/growing
Euler branch established by (8)-(12). Its exact cell equation, material
displacement and pressure response are useful positive mathematical
results. They do not refute positive acoustic response for a different
multiwave, cellular, vortex-array or other declared stationary ensemble.
The parent objective and that expanded candidate class remain active.
