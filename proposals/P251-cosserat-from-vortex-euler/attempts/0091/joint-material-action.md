# A positive core-angle action with its actual complete-fluid spin

## 1. One field, one relative displacement and computed spin normalization

Use 0090's single smooth bounded stationary field
`u=u_E+epsilon b`, `curl u=lambda u`, `p=-rho|u|²/2`, retaining its EPS
core and its actual compact positive material cage. Fix a unit marked axis
`e` and a core point `c`. Put `R_e(x)=e cross (x-c)` and choose a smooth
radial cutoff `chi`, with `0<=chi<=1`, equal to one on the unit ball and
zero outside the radius-two ball, nontrivial on the intervening shell.

Define the material rotational lift

    B_L(x)=chi(|x-c|/L) R_e(x).

It is exactly divergence free because the radial cutoff gradient is
orthogonal to `R_e`. For sufficiently large finite `L`, it rotates the
whole observed core rigidly. Define the COMPLETE compact perturbation rows

    A(Z)=rho integral R_e.Z,
    M(Z)=rho integral |Z|²,
    K(Z)=integral Z.Hess(p).Z-rho|(u.grad)Z|².

These integrals include ambient fluid. They are not tube-only angular
impulse or the old coadjoint inertia. With `s=(x-c)/L`, let

    a_chi=integral chi |e cross s|²,
    m_chi=integral chi² |e cross s|²,
    d_chi=a_chi-m_chi=integral chi(1-chi)|e cross s|²>0.

Then `A(B_L)=rho a_chi L^5`, `M(B_L)=rho m_chi L^5`. The difference
`D_L=rho d_chi L^5` is a strictly positive geometrical moment, not a fitted
inertia allowance.

Translate the fixed compact positive cage `X` of 0090 to a sufficiently
distant periodic lattice point, disjoint from `B_L`. Its actual background
stiffness has the fixed strict lower bound `k_*>0`, its mass is `m_X>0`,
and its axial spin row is `a_X`. The latter is translation independent:
every compact divergence-free field has `integral X=0`, by integrating
`div(x_i X)`. The actual value of `K(X)` retains the EPS correction; it is
not set equal to the limiting periodic value.

Set

    t=(a_X+sqrt(a_X²+4 m_X D_L))/(2 m_X),
    Z=B_L+t X.                                             (1)

Disjoint support gives `M(Z)=M(B_L)+t²m_X`, `A(Z)=A(B_L)+t a_X` and
`K(Z)=K(B_L)+t²K(X)`. The defining quadratic equation
`m_X t²-a_X t=D_L` proves the exact equality

    j:=M(Z)=A(Z)>0.                                        (2)

This chooses a material ensemble geometry by an observable identity. No
desired frequency, shear modulus, fitted constant or geometric rotor mass
has been supplied.

## 2. The normalization also retains strictly positive actual stiffness

Let `U0=||u||_infinity` and
`c_grad=integral ||D[chi(|s|) e cross s]||_F²`. Twice integrating the
pressure Hessian by parts, using `div B_L=0`, gives

    integral B_L.Hess(p).B_L=integral p tr[(D B_L)²].

Therefore, with the chosen pressure constant,

    |K(B_L)| <= (3rho/2) U0² c_grad L³.                     (3)

The apparent `L^5` pressure estimate would lose this exact cancellation.
If `D_L>=a_X²/m_X`, formula (1) gives
`t>=sqrt(D_L/m_X)/2`, and hence

    K(Z)>=k_* D_L/(4m_X)-(3rho/2)U0² c_grad L³.

A finite `L` satisfying
`L²>6 U0² c_grad m_X/(k_* d_chi)` and the preceding moment inequality
therefore gives

    kappa:=K(Z)>0,    j=A(Z)=M(Z)>0.                        (4)

The cage can be placed farther away as `L` increases; 0090's uniform
positive lower bound and `m_X,a_X` are unchanged by that selection. The
entire background remains one Euler field. All separations and amplitudes
are declared finite coherence geometry.

The exact identity useful to 0093 is, with `r_X=K(X)/m_X`,

    r_X-kappa/j=(r_X M(B_L)-K(B_L))/j>0                    (5)

for sufficiently large `L`. It follows directly from the complete mass and
stiffness expressions, independently of solving the normalization root.

## 3. Physical angle, actual material spin and common rotation

Because `Z=R_e` in an open core neighborhood, its volume-preserving flow
`h_q=exp(q Z)` is an exact rigid core rotation there for a nonzero local
angle interval. Let `R_beta` rotate the common coherence geometry and the
background about `e`, including its cage and ambient. The material map
`R_beta composed with h_q` gives the observed core angle

    Phi=beta+q.

This is an actual core geometry observation, not a relabeling of a normal
coordinate. The background, material pressure and statistical orientation
law rotate together; a static common rotation is a Euclidean symmetry of
the complete action. Common rigid affine rotations use their actual
boundary/mean terms, not an artificial finite isolated spinning cell.

For the relative material displacement `xi=Z q`, the actual Eulerian
velocity variation is

    v=Z qdot+[u,Z]q.

It therefore transports the core and boundary tags by 0084's exact material
identity. Define its complete physical angular-momentum variation using
any finite material parcel containing the supports over the time interval
under consideration. Such a parcel exists because `u` is bounded. The
displacement vanishes near its boundary; `integral Z=0` cancels the
centroid terms. Thus the actual axial angular-momentum response is

    delta S_e=A(Z) qdot+C(Z) q,
    C(Z)=2rho integral e.(Z cross u).                        (6)

This is independent of adding more unperturbed ambient labels to the parcel.
Its background parcel need not be a stationary spherical domain. No
nonvanishing boundary term has been dropped: the compact displacement
vanishes at that boundary, while all fluid around the core/cage is included.

The full material Jacobi action is varied before pairing `u` with `-u` and
matching their material lift. The pressure, `M` and `K` are even under that
pairing; `C` and every transport gyroscopic coefficient are odd. Their
cancellation is consequently an averaged action and observable identity,
not a deletion in each realization. Average also over rotations about `e`
of the complete marked field/lift to remove transverse spin rows. The
remaining averaged relative spin is precisely `j e qdot`.

For a local affine common velocity `Udot`, translation cross terms vanish
because `integral Z=0`. The symmetric affine kinetic cross also vanishes:
`integral (r_i Z_j+r_j Z_i)=0`, by integrating `div(r_i r_j Z)`.
Only the common rotational rate contributes, giving

    rho integral Udot.Z=j beta_dot

after the axial average. Thus the microscopic/common kinetic block is

    T=T_affine+j qdot beta_dot+j qdot²/2
      =T_affine-j beta_dot²/2+j Phi_dot²/2.                  (7)

This is not an added mass: it is an orthogonal decomposition of the actual
velocity norm. Indeed `A(Z)=M(Z)` means `R_e-Z` is L²-orthogonal to `Z`.
The corresponding velocity is
`(R_e-Z) beta_dot+Z Phi_dot`; its core sees `Phi_dot`, and the latter
component's actual axial angular momentum is `j Phi_dot`. The remaining
affine angular momentum belongs to the retained ambient component. The
finite enclosing affine inertia exceeds `j` by Cauchy–Schwarz, so the
subtracted term in (7) does not contradict positivity of the full metric.

This identifies complete-fluid mode spin, not the isolated EPS tube spin.
It also does not identify the common collective coordinate with a chosen
Eulerian point mean; a physical mean operation retains its own map below.

## 4. Full local averaged action and genuine mean-shear terms

For a marked scalar mode at number density `nu`, average the jointly rotated
field/lift and axis with the isotropic measure. Its second moment is
`E[e tensor e]=I/3`. The relative vector angle is
`q=Phi-beta`, `beta=curl U/2`. Let

    j_bar=nu j/3,    kappa_bar=nu kappa/3.

The paired, isotropic zeroth microscopic-gradient block becomes

    L= rho|Udot|²/2
         +j_bar|Phi_dot|²/2-j_bar|beta_dot|²/2
         -kappa_bar|Phi-beta|²/2
         +L_mean_strain+L_gradient+L_reaction.                (8)

All fluid, including ambient, contributes to the leading translation mass
`rho`: a uniform common translation moves every particle with that velocity
and is orthogonal to the compact relative mode. The density is not multiplied
by a tube filling fraction. The negative second-gradient kinetic term shown
in (8) is retained; 0093 supplies the other actual gradient masses. One
cannot discard it when computing the optical `k²` correction. In the old
micropolar potential convention its positive relative coefficient would be
`alpha=kappa_bar/4`, a newly computed ensemble coefficient, not `L_v T/6`.

The mean material stiffness is not the old coadjoint affine shear. For a
homogeneous isotropic background law, `E Hess(p)=0` and
`E[u_i u_j]=U_*² delta_ij`, so its raw quadratic form is

    K_mean(U)=-rho U_*² integral |grad U|².

For the compatible compact/periodic incompressible macro field, the
boundary-null-Lagrangian identity is
`integral |grad U|²=2 integral |sym grad U|²`. This is compatible with
rigid-rotation objectivity; a nondecaying affine rotation cannot silently
discard its boundary terms. Main's 0094 adds an explicitly declared,
independent positive STF material-cage population and retains the linear
affine/cage cross. Its resulting
`mu_total=-rho U_*²+b_s t_s+a_s t_s²`, `a_s>0`, can be made positive by a
finite geometrically specified `t_s`, with its full strain-rate inertia
retained. Equation (8) imports that material block, not 0057's old shear.

Isotropy removes a quadratic STF-strain/axial-vector local cross: these
inequivalent rotation representations have no invariant bilinear pairing.
When nonlocal, gradient or unaveraged reaction blocks are retained, their
cross terms remain in `L_gradient+L_reaction` until the actual full Schur
reduction. Disjoint cage supports eliminate local mass and stiffness
crosses, not an arbitrary nonlocal reaction inverse.

Spin population removal `nu->0` at fixed positive affine population removes
`j_bar,kappa_bar` and leaves 0094's declared Cauchy sector. Removing every
vorticity/attached structure instead gives the neutral incompressible Euler
limit. These are different physical population limits, not a universal
five-scalar filament formula or an interchangeable limit operation.

## 5. Actual mean and all reaction coordinates remain part of the same action

Use 0087's defined volume-preserving mean if a GLM coordinate is required:
`g_mean=C_F composed with E[g]`, together with its exact defining constraint.
For `g=Xi composed with g_mean`, the exact kinetic blocks are
`rho E[(D Xi)^T D Xi]`, `rho E[(D Xi)^T Xi_,a]`, and
`rho E[Xi_,a.Xi_,b]`. These contain the present material lift and every
retained shape, pressure and mean-gauge variation. The full blocks, not
just `j`, are transformed by the actual derivative of this mean map.

Equation (7) establishes the physical common-rotation/core-angle block in
its stated affine coordinate. It does not prove that the GLM covariance
spin equals `j_bar Phi_dot` by changing its name. That observable is
`rho E[eta cross eta_t]` for the actual registered fluctuation field `eta`;
0087's physical mean observation also contains its covariance-rate term.
Both must be evaluated for the actual evolving registered family. The
full material spin (6) is independently meaningful and remains part of the
angular-momentum balance even when those observables differ.

For any retained linear reaction coordinates `z`, one can explicitly put
`xi=U+Zq+sum Z_a z_a` in the SAME 0084 Jacobi form. Its entries are

    M_ab=rho integral Z_a.Z_b,
    G_ab=rho integral Z_a.(u.grad Z_b),
    K_ab=integral Z_a.Hess(p).Z_b
                  -rho integral (u.grad Z_a).(u.grad Z_b).

The complete finite or operator block, its pressure boundary conditions and
mean gauge determine the constrained Schur action. Equations (2)–(8) do
not set that reaction operator to zero. A Schur change to the angle field
or physical spin row must be propagated to both, not only to its energy.

## 6. Exact Euler/Kelvin license and the precise remaining join

All maps used above are genuine volume-preserving material maps. Their
transported core geometry is exact, and (8) is a pullback of the full Euler
Jacobi action under the declared Cauchy–Born material restriction. This is
stronger than assigning an oscillator energy: the pressure, convection,
mass, spin, core observation and finite-core ambient field are constructed.

However a restriction of variations is not itself the assertion that every
solution of the restricted equations solves all unreduced Euler equations.
In particular fixing each material Kelvin one-form requires

    [xi_t+(u.grad)xi+(Dxi)^T u]=0 modulo exact forms,          (9)

with periods retained. On the present global Beltrami domain this becomes
`xi_t=(lambda P-curl)(xi cross u)`. The q-only material family is not
automatically on that fixed leaf: at `q=0` its arbitrary `qdot Z` generally
has nonzero circulation. Passive relabeling covariance of a specified
reference field cannot repair this missing dynamical equation.

A useful exact bridge, derived without borrowing the old inertia, is

    H_orbit(xi)-K_material(xi)
        =rho||(lambda P-curl)(xi cross u)||² >=0.            (10)

Here `H_orbit=rho[lambda²||P(xi cross u)||²
-lambda<P(xi cross u),curl(xi cross u)>]`. Thus the constructed positive
material cage is also a positive actual orbit-Hessian direction on this
same field. This supplies the next phase-space/relabeling construction with
a positive direction; it does not equate its reduced inertia with `M(Z)`.

For a full fixed-Kelvin dynamical conclusion, the actual phase-space shape
and mean/relabeling reactions must now realize (9), and their constrained
action and physical observation map must be evaluated together. Alternatively
a explicitly conditional material Cauchy–Born model may use (8), while
stating the retained circulation/reaction law rather than claiming it has
already derived unrestricted Euler trajectories. No new all-k or nonlinear
invariant-manifold condition is imposed here. The missing condition is the
specific material circulation equation already required by the Euler model.

## Verdict

`route_verdict: established` for the constructive moment-normalized positive
material lift, complete-fluid physical spin, common/relative kinetic block,
and stated conditional full material action. `evidence_scope:
REPRESENTATION_SCOPED`. Fixed-Kelvin phase-space/GLM reaction closure remains
an explicit parent construction, so the full requested Euler–Cosserat
objective is not declared complete by this action attachment.
