# Actual Euclidean symmetry, its rate response, and complete current

This attachment establishes the exact symmetry/current part of candidate A
and executes its full-pressure rate correction. It does not infer an
acoustic mode from a rigid stationary family. The inputs are the reviewed
0211 ring, 0206's Euclidean impulse/Jordan sector, and 0222's actual Kelvin
optical preparations.

## 1. Static and rate variations are different Euler data

Write T=u.grad, C xi=-T xi+Du xi, and
L w=-P[T w+Du w]. Lin reconstruction is xi_t=C xi+w. For the rigid
generator K_a(x)=a cross x, a constant vector a,

    R_a=[u,K_a]=a cross u-Du K_a,
    L R_a=0,                                             (1)

because rotation of the complete stationary Euler field is a stationary
family. The constant far velocity rotates with it. Equation (1) is not
a claim that rotation of only one component of an interacting field is
stationary.

The proposed rate lift xi=t K_a instead requires w=K_a+t R_a. Its exact
Euler discrepancy is

    L K_a-R_a=-2 P(a cross u),
    curl(2 a cross u)=-2(a.grad)u.                        (2)

Thus a pressure gauge cannot generally remove the discrepancy. These
identities can be evaluated in an affine-plus-decaying class, fixing the
constant far acceleration by a linear pressure, or after a declared
compact divergence-free cutoff and its additional return. No L2 Euler
semigroup on the unbounded affine field K_a is asserted.

The actual rate response is obtained by evolving the initial Euler field
w(0)=K_a in that declared class. Relative to the proposed response, its
decaying/non-affine correction is

    delta w(t)=-2 integral_0^t exp((t-s)L)P(a cross u) ds,
    delta xi(t)=integral_0^t exp((t-s)C)delta w(s) ds.     (3)

In a finite-energy cutoff realization the forcing in (3) is replaced by
the exact cutoff residual, including its ambient return. The formula
retains the entire pressure operator. In particular the local Taylor
rows are xi_tt(0)=-2P(a cross u) and, using the stationary rotational
Jacobi zero-mode equation,

    xi_ttt(0)=4 P T P(a cross u).                         (4)

Time reversal makes the first row odd in u but the second even. Whole-law
time reversal alone therefore does not turn (3) into t K_a.

An exact diagnostic makes this distinction transparent. For the actual
uniform-rotation stationary family u_epsilon=(B+epsilon A)x, B and A
constant antisymmetric matrices, the Euler velocity w=A x is stationary
and its material displacement is

    xi(t,x)=[integral_0^t exp(sB) A exp(-sB) ds]x.         (5)

The pressure is the quadratic potential whose Hessian is
-(B+epsilon A)^2. If B rotates about ez at rate Omega and A tilts about
a perpendicular axis, a stationary axisymmetric material tag sees a
rotation vector integral_0^t Rot_ez(Omega s)a ds. Pairing the two signs
of B gives sin(Omega t)a/Omega, not t a. Its cubic discrepancy is
-Omega^2 t^3 a/6. This is an exact Euler example of the even row in (4),
not a proposed global background for the ring construction.

## 2. Actual ring Jordan response and its impulse

Use the laboratory travelling ring v_*(x-Ut ez) and its relative steady
field u=v_*-U ez. The exact family

    v_{Q,b}(t,x)=Q v_*(Q^T(x-b)-Ut ez)

has, in the original translating chart,

    delta X=b+Ut delta N,  delta N=epsilon cross ez.       (6)

Here delta N is constant. The complete Kelvin impulse is I=I0 N and
the Euclidean phase pairing is
Omega(T_b,K_epsilon)=I0 b.dot(epsilon cross ez). The actual reduced
quadratic action is the 0206 action

    p.dot X_t + J0/(2I0^2)(p_x p_y,t-p_y p_x,t)
                                  -|p|^2/(2M_eff),
    p=I0 delta N,  M_eff=I0/U>0.                         (7)

It supplies an independent translational Jordan sector and constant
physical tilt. It does not supply a linearly increasing whole-ring tilt.
The latter would rotate the conserved nonzero impulse unless another
part of the actual perturbation carries the compensating impulse.

The tag's material mass and mechanical spin are not I0 and J0. In
particular assigning rho to M_eff per cell without the actual cell mass
and its exterior return would change the observable and the action.

## 3. Complete Noether moments versus tag moments

For a smooth compact isovortical generator xi let F=xi cross omega. The
complete vorticity impulse and angular impulse, with these conventions,
are I=rho/2 integral x cross omega and
J=-rho/2 integral |x|^2 omega. Their exact variations are

    delta I=rho integral F,
    delta J=rho integral x cross F.                      (8)

They are also the translation and rotation rows of the KKS form
rho integral omega.dot(eta cross xi). Boundary terms vanish for the
compact vorticity variation; a physical outer boundary instead retains
those terms.

For the actual positive stationary material tag, use the distinct rows

    delta X=M_tag^-1 integral rho chi xi,
    G=integral rho chi x cross xi,
    S=integral rho chi[x cross xi_t+2 xi cross u],
    delta I_shape,ij=integral rho chi(x_i xi_j+xi_i x_j).

They obey G_t-S=-2 integral rho chi xi cross u. The exterior reaction and
shape row are retained when collapsing the tagged parcel to its centroid.
Neither a varying tag spin nor its first spatial dipole identifies it
with the complete conserved J in (8). Equations (3), (6), and these rows
give actual response data for a future common-K preparation, with no
assumed equality between canonical impulse and the measured centroid.

## Route result and continuation

The static Euclidean response, positive impulse/Jordan action, and full
pressure rate-response formula are established. The particular inference
"stationary rigid rotation implies xi=t K is Euler" is refuted by (2)
and the exact example (5). This is not a refutation of an acoustic tag
response with ambient reaction. Candidate C executes a separate
density/accuracy and full quadratic-normalization repair; an actual
common-K acoustic preparation remains the downstream construction.
