# Balanced first shape correction and the remaining finite-radius inverse

This calculation uses the activated 0255 contract, the exact 0253 straight
compacton, and the independently derived 0256 translation balance.  It is an
analytic continuation step, not a finite-radius existence theorem.  Write
`epsilon=1/R`, fix the simple balanced member `Phi=Phi_eta*` of the frozen
profile family, let its support radius be `a`, and put

\[
 y(\rho)=\Phi'(\rho)<0\quad(0<\rho<a),\qquad
 L=-\Delta_2-\lambda^2+g'(\Phi).                       \tag{1}
\]

The first curvature equation is

\[
 Lh=-f,\qquad
 f=(y+2\rho g(\Phi))\cos\theta .                       \tag{2}
\]

The singular abstract coefficient `g'(Phi)` is retained in (1).  All
regularity conclusions below concern its spatial compositions and the
defining function, not smoothness of `g` at label zero.

## 1. Exact quadrature in the broken-translation sector

Differentiating the straight radial equation gives

\[
 L_1y=0,qquad
 L_m=-\partial_\rho^2-\rho^{-1}\partial_\rho
       +m^2\rho^{-2}-\lambda^2+g'(\Phi).               \tag{3}
\]

Set

\[
 h(\rho,\theta)=y(\rho)v(\rho)\cos\theta .             \tag{4}
\]

Reduction of order is exact:

\[
 L_1(yv)=-\frac{(\rho y^2v')'}{\rho y}.
\]

Thus (2) is equivalent on `0<rho<a` to

\[
 (\rho y^2v')'=\rho y\,[y+2\rho g(\Phi)],             \tag{5}
\]

and hence

\[
 I(\rho):=\int_0^\rho s\Phi'(s)
             [\Phi'(s)+2s g(\Phi(s))],ds,
 \qquad
 v'(\rho)=\frac{I(\rho)}{\rho\Phi'(\rho)^2}.          \tag{6}
\]

The additive constant in `v` adds a multiple of `Phi_x`; it is precisely the
radial-center gauge.  In polar coordinates the 0256 projection is

\[
 \langle\Phi_x,f\rangle
 =\pi\int_0^a\rho\Phi'(\rho)
       [\Phi'(\rho)+2\rho g(\Phi(\rho))],d\rho
 =\pi I(a).                                            \tag{7}
\]

Consequently the exact scalar balance is not only a Fredholm pairing: it is
the endpoint condition for the quadrature (6).

At the center, the exact Bessel core gives
`Phi'=-kappa rho+O(rho^3)`, `kappa=lambda^2 Phi(0)/2`, and `g(Phi)=0`.
Therefore `I(rho)=kappa^2 rho^4/4+O(rho^6)` and

\[
 v'(\rho)=\rho/4+O(\rho^3).                            \tag{8}
\]

Thus (4) is a smooth Cartesian `m=1` function at the center.

At the flat edge put `d=a-rho` and use the exact common collar

\[
 \Phi=C e^{-1/d},\qquad
 \Phi'=-\Phi d^{-2},
\]

\[
 g(\Phi)=\Phi\left[\lambda^2+d^{-4}-2d^{-3}
                    -\frac{d^{-2}}{\rho}\right].      \tag{9}
\]

The integrand in (6) has the asymptotic expansion

\[
 \rho\Phi'[\Phi'+2\rho g(\Phi)]
 =-2a^2C^2e^{-2/d}d^{-6}[1+O(d)].                      \tag{10}
\]

If `I(a)` were nonzero, (6) would grow like
`I(a)d^4e^{2/d}/(aC^2)`, excluding a physical flat-edge tangent.  At the
balanced profile `I(a)=0`, however,

\[
 I(\rho)=-\int_\rho^a s\Phi'(s)
              [\Phi'(s)+2s g(\Phi(s))],ds
 =a^2C^2e^{-2/d}d^{-4}[1+O(d)],                        \tag{11}
\]

so that

\[
 \boxed{\lim_{\rho\uparrow a}v'(\rho)=a}.             \tag{12}
\]

The correction `h` and every spatial composition appearing in (2) are flat:
polynomial losses in `d^-1` are dominated by `e^-1/d`.  Thus their zero
extensions are smooth.  This establishes a regular first physical correction
in the `m=1` sector.  It does not make the additive approximation
`Phi+epsilon h` uniformly positive up to the moving edge.

## 2. The same quadrature supplies the first free-boundary shape

The nonuniformity is removed in the 0253 defining-function coordinate.  Since

\[
 \frac{d}{dT}e^{-1/T}=e^{-1/T}T^{-2},
\]

the additive correction (4) is the linearization of

\[
 T_\epsilon(\rho,\theta)
 =d-\epsilon v(\rho)\cos\theta+O(\epsilon^2),
 \qquad \phi_\epsilon=e^{-1/T_\epsilon}.               \tag{13}
\]

Indeed `delta phi=Phi d^-2(-v cos theta)=Phi'v cos theta`.
The boundary displacement `v(a)cos(theta)` depends on the radial-center gauge,
as it should.  Its normal derivative does not.  On `T_epsilon=0`,

\[
 |\nabla T_\epsilon|
 =1+\epsilon v'(a)\cos\theta+O(\epsilon^2)
 =1+\epsilon a\cos\theta+O(\epsilon^2).               \tag{14}
\]

The exact transformed finite-radius equation requires

\[
 |\nabla T_\epsilon|=r/R=1+\epsilon x.
\]

At the perturbed boundary `x=a cos(theta)+O(epsilon)`, so (14) matches the
full eikonal boundary condition through first order.  This is a second,
independent interpretation of the value `v'(a)=a`: the translation balance
selects exactly the shape derivative compatible with physical flatness.

## 3. What part of the linear inverse is now controlled

For compactly supported admissible radial amplitudes `w=yq`, integration by
parts and (3) give the ground-state factorization

\[
 \int_0^a\left[(w')^2+
   (\rho^{-2}-\lambda^2+g'(\Phi))w^2\right]\rho\,d\rho
 =\int_0^a\rho y^2(q')^2\,d\rho.                      \tag{15}
\]

Thus the physical `m=1` quadratic form is nonnegative and its regular null
direction is the translation `y`.  The second reduction-of-order homogeneous
solution grows outside the flat physical domain.  For every `m>=2`,

\[
 \mathfrak q_m[w]=\mathfrak q_1[w]
 +(m^2-1)\int_0^a\frac{w^2}{\rho}\,d\rho>0            \tag{16}
\]

for nonzero regular `w`.  Therefore all higher angular modes are coercive at
the quadratic-form level; the `z` member of `m=1` is removed by evenness and
the `x` member is handled by the simple profile border and (6).

At this stage two infinite-dimensional rows were load-bearing:

1. the `m=0` radial operator had not yet been proved to exclude zero outside
   its permitted negative spectrum, and (15) gives no radial inverse;
2. no weighted Banach scale has yet been built in which the nonlinear map in
   `T`, the interior positive-label equation, the boundary graph, and the
   profile border are differentiable with a bounded right inverse.

The first-order construction eliminates neither row.  In particular,
ordinary unweighted implicit-function theorems do not apply merely because
the tangent (13) is smooth: `g'(Phi)` diverges at the edge and a moving flat
support is invisible to relative-error norms in `phi`.

The append-only continuation in `radial-nondegeneracy.md` closes the first row
by an explicit three-border profile construction preserving the core and
`Q=0`.  The minimum Candidate-I repair is therefore now the second row alone:
construct the right inverse in a defining-function/flat-field space whose
`m=1` component agrees with (6) and whose boundary trace enforces (14).  An
interior elliptic bootstrap through the exact `g=0` core must give `C^4` core
convergence on a neighborhood of one fixed regular inner annulus.  That
regularity is sufficient to consume the checked 0258 Bessel-core flux-action
twist; no second profile parameter is required for twist merely because `eta`
is already used by translation balance.

## 4. Exact variational diagnosis

The first variation of the frozen weighted functional is

\[
 \delta\mathcal E_{R,\eta}(\phi)[\xi]
 =\int\left[-\operatorname{div}\left(\frac{R^2}{r}
     \nabla\phi\right)-\frac{R^2\lambda^2}{r}\phi
     +r g_\eta(\phi)\right]\xi\,dr\,dz,               \tag{17}
\]

and multiplication by `r/R^2` gives exactly the finite-radius equation.
This confirms the sign and weight, but it does not supply a minimizer.

There are two concrete costs hidden by a direct-minimization slogan.

First, inverse profile design specifies `g_eta` only on the label interval
reached by `Phi_eta`.  A global variational class also samples labels above
`Phi_eta(0)`.  Coercivity and even the definition of the global functional
therefore require a declared high-label extension.  A box constraint at the
core value is not free: one must prove that it is inactive before the
minimizer is an Euler critical point.  A coercive high-label extension must
likewise be accompanied by an a priori bound keeping the selected critical
point inside the unchanged physical label interval.

Second, the straight two-dimensional limit is not an unconstrained minimizing
problem.  Put

\[
 W(s)=G(s)-\frac{\lambda^2}{2}s^2,
 \qquad
 \mathcal E_0[u]=\int_{\mathbb R^2}
       \left(\frac12|\nabla u|^2+W(u)\right).          \tag{18}
\]

The straight Pohozaev identity is `int W(Phi)=0`.  Near label zero,
`W(s)~s^2 log(1/s)^4/2>0`.  Since the edge occupies a set of positive measure,
the zero integral forces `W(s_0)<0` for some label attained in the transition
or core.  A disk plateau of value `s_0`, with a fixed-width smooth boundary
layer, has negative potential energy of order its area and gradient/boundary
cost of order its perimeter.  Hence

\[
 \inf_{H^1(\mathbb R^2)}\mathcal E_0=-\infty.          \tag{19}
\]

There is also a local obstruction specific to the desired compacton.  For
`Phi_k(x)=Phi(x/k)`, two-dimensional scaling and `int W(Phi)=0` give

\[
 \mathcal E_0[\Phi_k]
 =\frac12\int|\nabla\Phi|^2+k^2\int W(\Phi)
 =\mathcal E_0[\Phi]                                  \tag{20}
\]

for every nearby `k`.  If `Phi` were a local unconstrained minimizer, every
nearby equal-energy `Phi_k` would also be a local minimizer and hence a
critical point.  But direct substitution yields

\[
 -\Delta\Phi_k-\lambda^2\Phi_k+g(\Phi_k)
 =(k^{-2}-1)[\lambda^2\Phi-g(\Phi)](x/k),             \tag{21}
\]

which is nonzero for `k ne 1`; in particular it is nonzero in the retained
`g=0`, `Phi>0` Bessel core.  This contradiction proves that the desired
straight compacton is not a local unconstrained minimizer of (18).

Equations (19)--(21) refute only the direct unconstrained-minimizer subroute.
A Nehari, mountain-pass, or other constrained critical-point construction can
still be viable.  It must declare the high-label extension, prove constrained
compactness for the signed logarithmic potential and cylindrical weight,
show that the constraint is natural, select the balanced member, keep the
support away from the auxiliary wall and axis, and recover the prescribed
core and `C-infinity` flat edge.  No located theorem supplies that chain.

## 5. Exact conditional twist transfer from 0258

Attempt 0258 establishes on the retained core
`Phi=c J_0(lambda rho)` that the rotation number divided by the auxiliary
straight period scale and the velocity-flux action are

\[
 q(\rho)=\frac{J_1(\lambda\rho)}{\rho J_0(\lambda\rho)},
 \qquad J_u(\rho)=c\rho J_1(\lambda\rho),              \tag{22}
\]

with the exact nonzero central coefficients

\[
 \frac{dq}{dJ_u}(0)=\frac{\lambda^2}{8c},
 \qquad
 \frac{dq}{dJ_\omega}(0)=\frac{\lambda}{8c}.          \tag{23}
\]

For an actual finite-radius 0255 solution, the corresponding normalized
quantities on a regular level `C_R(s)={phi_R=s}` and its enclosed disk are

\[
 \frac{\operatorname{rot}_R(s)}R
 =\frac{2\pi}{\lambda s\displaystyle
   \int_{C_R(s)}\frac Rr\frac{d\ell}{|\nabla\phi_R|}},
 \qquad
 J_{u,R}(s)=\frac\lambda{2\pi}\int_{D_R(s)}
                  \frac Rr\phi_R\,dx\,dz.             \tag{24}
\]

These are contour and area functionals on one fixed regular annulus, not a
trajectory comparison over a period growing with `R`.  `C^4` convergence of
the constructed core gives `C^1` convergence of the normalized rotation as a
function of `J_{u,R}`; choosing the annulus where (23) is bounded away from
zero preserves the twist sign.  This is a conditional consequence of a future
0255 solution, not evidence that the solution exists.  It sharpens the
nonlinear repair: `eta` is used only to solve the translation row, while the
unchanged Bessel core supplies twist without further tuning.

## Route verdict and supplier license

`candidate_I_route_verdict: blocked with a substantially reduced missing
construction`.  The simple balance, exact regular `m=1` correction, its
flat-edge asymptotics, the first boundary shape/eikonal match, coercivity of
every `m>=2` sector, and constraint-preserving removal of a radial zero are
established.  The sole remaining load-bearing step is a nonlinear right
inverse in compatible defining-function/flat-field spaces.

`candidate_II_route_verdict: blocked after refutation of its unconstrained
direct-minimizer subroute`.  The weighted first variation is exact; the
remaining viable representation is a declared constrained critical-point
construction with high-label, compactness, topology, matching, and smoothness
costs above.

`evidence_scope: exact first-order finite-radius shape construction and
representation-scoped variational diagnosis`.  There is still no actual
finite-`R` compact Euclidean Euler field, no positive measured density, and no
response/action/current join.  Neither route verdict propagates to the parent
obligation or establishes exhaustion.

If Candidate I is completed with the stated `C^4` core convergence, it licenses
one smooth, zero-far-flow, compact-velocity stationary Euler torus with the
0255 pressure, volume and energy formulas.  The checked 0258 implication and
(22)--(24) then license nonzero flux-action twist on a positive-width regular
inner annulus without another profile parameter.  Disjoint compact-support
copies give the exact positive volume-density normalization once the measured
tube observable is integrated over that annulus.  This still does not license
the 0250 same-field response, physical action/current, coupled histories, or
parent completion.
