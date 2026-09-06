# A global finite-energy twisted Euler initial state and its exact spinor lift

This explicit construction is the registered0017 candidate. It gives a smooth
localized initial state and its smooth-local-Euler time development. It is not
claimed to be a steady solution or nonlinearly stable particle. All quantities
below belong to the same physical velocity after the Hodge projection.

## 1. Smooth profile, Berry field and divergence

Choose a C-infinity nonincreasing F(r) with 0<=F<=pi, F=pi for r<=a,
F=0 for r>=b, and a<b. Choose it strictly decreasing on a nonempty interval.
Smooth flat splicing exists by integrating a nonnegative compact smooth bump
on (a,b) and normalizing its integral to pi. Write n=x/r and define

    z0=(cos F+i n_z sin F, (-n_y+i n_x)sin F), |z0|=1.

It is constant near the origin and outside the outer sphere, so no apparent
1/r term is singular. With kappa>0 define v^flat=-i kappa z0^dagger dz0.
Direct differentiation gives the axisymmetric one-form

    v^flat/kappa=cos(theta) F' dr
       -sin F cos F sin(theta) dtheta
       +sin^2 F sin^2(theta) dphi.                       (1)

Its Cartesian vector is

    v/kappa=F' n_z n+(sin F cos F/r)(e_z-n_z n)
                         +(sin^2 F/r)(e_z cross n).

The vector v is smooth and compactly supported, but is generally not
divergence-free. Specifically

    div v=kappa cos(theta) S(r),
    S=F''+2F'/r-sin(2F)/r^2.                            (2)

Using v itself as Euler velocity would therefore be a false transfer.

## 2. Execute the physical Hodge inverse

Solve Delta phi=div v with decay, writing phi=kappa h(r)cos(theta).
The exact l=1 Poisson Green formula is

    h(r)=-(1/3)[r^-2 integral_0^r s^3 S(s)ds
                          +r integral_r^infinity S(s)ds]. (3)

Differentiation proves h''+2h'/r-2h/r^2=S. Since S vanishes near zero,
h is proportional to r there and h cos(theta) is smooth Cartesian linear
data. Outside b, h=-C/(3r^2), where

    C=integral_0^infinity r^3 S(r)dr
     =integral_0^infinity r[2F(r)-sin(2F(r))]dr>0.        (4)

Both boundary integrations in(4) are valid because r^2 F tends to zero at
zero and F is zero for r>=b. The interval F=pi near the origin contributes
to the second expression even though S itself is supported in the shell.

Define the ACTUAL initial Euler state

    u0=v-grad phi,  omega0=curl v,
    z=exp(-i phi/kappa) z0.

Then exactly

    div u0=0,  u0^flat=-i kappa z^dagger dz,
    omega0 is smooth and compactly supported.           (5)

The exterior velocity is the exact dipole

    u0=(kappa C/(3r^3))[e_z-3n_z n], r>b,
    I=(1/2)integral x cross omega0 dx=-(4pi kappa C/3)e_z. (6)

The last identity follows either by a direct compact-vorticity moment integral
or comparison with the exact leading Biot--Savart formula0007. Thus u0 is smooth
on all R3, decays as r^-3, has finite kinetic energy and belongs to H^s for
every finite s. The exterior velocity is not artificially set to zero.

## 3. Same-field helicity and axial angular momentum

Compute the three-form from(1) rather than assigning a knot number:

    v^flat wedge d v^flat
        =2 kappa^2 F' sin^2 F sin(theta) dr wedge dtheta wedge dphi.

The phase correction changes it only by the exact form -d(phi d v^flat).
The vorticity is compact and the boundary term vanishes, giving

    Hel(u0)=integral u0 dot omega0 dx
           =8pi kappa^2 integral_0^infinity F' sin^2 F dr
           =-4pi^2 kappa^2.                             (7)

The sign is fixed by the particular components of z0 and the outward
Cartesian orientation. Reversing the relevant spatial orientation reverses
helicity; an arbitrary sign must not be inserted after construction.

The Hodge phase has no azimuthal component, so the physical swirl is

    u_phi=kappa sin^2 F sin(theta)/r.

Consequently the absolutely convergent AXIAL component of intrinsic angular
momentum about the construction center is

    j_z=rho integral (x cross u0)_z dx
       =(8pi rho kappa/3) integral_0^infinity r^2 sin^2 F dr>0. (8)

This proves a nonzero finite classical rotation moment on the same helical
flow. It does not assert absolute convergence of the complete vector moment:
the dipolar meridional tail can fail that stronger condition. A full rotational
KKS sphere needs its generator/boundary pairing justified separately. The
degree of z0 on the one-point compactification is -1 in the convention fixed
by(7); the decaying phase homotopy preserves that degree. This is a topology
of an actual velocity representation, not a statistics postulate.

## 4. Exact Euler time transfer of the global representation

Smooth-local-Euler existence on R3 applies to u0. On its smooth existence
interval let X(t,a) be the actual volume-preserving fluid flow and let

    z(t,X(t,a))=exp[(i/kappa) integral_0^t
                        (|u|^2/2-p/rho)(s,X(s,a))ds] z(0,a). (9)

This is driven by the actual Euler solution, not a selected future target.
It solves the nonlinear spinor transport equation0013,

    i kappa D_t z=(p/rho-|u|^2/2)z.

The one-form alpha=-i kappa z^dagger dz therefore obeys

    (partial_t+L_u)alpha=d(|u|^2/2-p/rho).

Euler's u^flat obeys the identical equation with identical initial data, so
transport uniqueness gives alpha=u^flat throughout this interval. This is a
global-in-space spinor lift of actual unsteady Euler, including its pressure.
Vorticity support is the material image of the original compact shell;
helicity is conserved. The statement does not prove global-in-time smoothness
or preservation of a compact shape/particle identity. Those are P2 tasks.

## 5. What scale and topology actually constrain

Keep a dimensionless shape F_*(r/L) fixed. Equations(3)--(8) imply

    H=rho kappa^2 L E_*,
    I=kappa L^2 I_*,
    j_z=rho kappa L^3 J_*,
    Hel=-4pi^2 kappa^2,                                 (10)

where E_*>0, J_*>0 and I_* is a fixed directed shape constant. In particular,
helicity alone supplies no positive whole-space kinetic-energy lower bound:
at fixed kappa and degree, shrinking L sends H to zero. This dilation changes
other transported labels and j, and is not automatically dynamically accessible
on one Euler coadjoint orbit. It is therefore NOT an instability proof.

If both nonzero j and helicity are fixed, (10) fixes kappa and L within this
chosen shape family. That is a classical conserved-sector scale relation,
not a universal action quantum: the values of those invariants were inputs,
and more general profiles carry additional freedom. We have not identified
the texture's degree with electric charge or a prequantum Chern class.

## Route result and continuation

The global smooth compact-vorticity finite-energy twisted initial carrier,
its finite nonzero axial moment, and its exact local-time Euler/spinor lift
are established as stated. This repairs the nonlocalized0003 dynamical example
and the local-only0013 coordinate example in one actual physical field.

The persistence route remains blocked on a controlled same-leaf neighborhood
and a restoring/modulation theorem for this candidate. The stationary
Gavrilov/Cao--Zhan alternatives in0014 are being tested in parallel. Physical
quantum probabilities, exchange, spin-half, electromagnetic and weak currents,
and relativistic dispersion remain explicit parent constructions; none is
inferred from(7), a two-component coordinate or the existence interval.

## 6. Execute the stationary/translation test before calling it a particle

The spherical profile makes the poloidal streamfunction particularly explicit.
With the usual axisymmetric convention
u_r=(r^2 sin(theta))^-1 partial_theta Psi and
u_theta=-(r sin(theta))^-1 partial_r Psi, set

    G(r)=r^2(F'-h')/2,  Psi=kappa G(r)sin^2(theta),
    G'=sin F cos F-h,
    xi=r sin(theta)u_phi=kappa sin^2 F sin^2(theta).     (11)

The identity for G' follows from the FULL radial Poisson equation(3).
For an axisymmetric steady field, angular Euler requires poloidal advection
of xi to vanish. On the transition shell away from its angular critical sets,
the Jacobian of Psi and xi consequently gives

    (sin^2 F)'G-sin^2 F G'=0.                           (12)

Choose the registered strictly decreasing bump profile with 0<F<pi throughout
(a,b). Equation(12) would force G=c sin^2 F there. Flatness at b would then
give G(b)=0. But(3)--(4) give G(b)=-C/(3b), which is nonzero. Thus this exact
initial family is not steady; the failure is the actual angular momentum
transport row, not a numerical residual.

A common axial translation at speed U replaces the streamfunction coefficient
by Q(r)=kappa G(r)-U r^2/2. Angular Euler would instead require Q=c sin^2 F.
Matching its flat edge requires BOTH Q(b)=0 and Q'(b)=0. The exterior formulas
give

    G(b)=-C/(3b), G'(b)=C/(3b^2),
    Q(b)=0 => U=-2kappa C/(3b^3),
    Q'(b)=kappa C/b^2 !=0.                              (13)

No common translation repairs it. This refutes the stationary or uniformly
axially translating realization of this radial Hopf--Hodge ansatz. It does not
refute persistent unsteady motion or arbitrary twisted Euler fields.

The failure-generated geometry repair is explicit: allow the amplitude and
phase surfaces to depend on two meridional coordinates so the swirl can be a
function of the actual translating streamfunction. The genuine Gavrilov and
Cao--Zhan carriers under0014 already satisfy that angular transport row and
are the current stationary candidates. A periodic/deforming version of the
global spinor carrier would instead need a controlled recurrent orbit and
nearby-state estimate. Those remain live routes; neither the initial-state
success nor the stationary-ansatz failure closes P2 or P4.
