# Actual common-V operator and whole-law physical current

The boundary was centrally validated before source reopening. Use0188's
same divergence-compatible lift, not a substituted Kelvin lift. The actual
field is u=(psi,A sinZ,-B sinY), psi=B cosY+A cosZ,
alpha=B cosY-A cosZ, curl u=-u, and

    T=A sinZ partialY-B sinY partialZ,
    H=-Delta, mathcalB=H-1, L=-P(u.grad+Du).

All averages are full normalized cell averages followed by the Haar
orthonormal-pair law for body components of one common lab (kappa,V).
The body field and all preparations transform together under O(3)/TR.
No rotated macro wave vector is selected independently per realization.

## 1. The same-lift forcing, including its pressure derivative

Let T_V=-(V.grad)u, q_V=kappa cross T_V/lambda. Expansion of the exact
Leray projection gives

    L_K V=T_V+ik L1 V+O(k²),
    L1 V=-grad Delta^(-1)(kappa.T_V)-P[(kappa.u)V].

Since div q_V=-kappa.T_V, the solenoidal remainder is

    b_V=L1 V-q_V=-P[(kappa.u)V+q_V].                 (1)

Writing the actual V velocity as
w_V=V+tT_V+ik[tq_V+z_V]+O(k²) yields

    z_V,t=L z_V+t d F0+b_V,
    d=kappaY VY-kappaZ VZ.                           (2)

The exact initial class eta_V=0, w_V=V+ik P_K r_V has z_V(0)=r_V.
The return r_V is a genuine changed initial circulation/velocity, with
mean zero; it is not a freely assigned observed force.

For0188's stationary forced z_D, y=z_V-tz_D obeys
y_t=L y+b_V-z_D. This equation remains the correct interface, but its
whole-law observation is what the physical closure constrains.

## 2. Exact planar and axial-plus-vorticity decomposition

Average in the invariant X direction first; this commutes with L and
the physical current. Write z_V=(b,-phi_Z,phi_Y),
r=b-Hphi. Let beta be the planar stream function of b_V and
r_b=(b_V)_X-Hbeta. The actual full source(1) gives

    H phi_t=-T mathcalB phi+t d Talpha+Hbeta,
    r_t+Tr=-(t d/2)Talpha+r_b.                       (3)

These are both equations of the THREE-component Euler velocity. In
particular r is not passive when its actual forcing is present.

Let e=kappaY VY+kappaZ VZ and phi_bar=E[d phi]. Direct exact pressure
and Haar calculation gives

    E[d beta]=alpha/5, E[e beta]=0, E[d²]=1/5.

Therefore

    H phi_bar,t=-T mathcalB phi_bar+(t/5)Talpha+alpha/5. (4)

The e-weighted psi moment is conserved. For every smooth first integral
g(psi), (4) gives the actual conserved-forcing row

    d/dt <g(psi)H phi_bar>=(1/5)<g(psi)alpha>.        (5)

This applies to arbitrary correlated initial returns; it is not a
per-orientation stationarity assumption.

## 3. Full passive-current forcing, not an omitted axial component

The physical current row is c=(kappa.u)V+(u.V)kappa. The complete
same-cell identity from0184 is

    <c z>=<c_X r>-e<psi phi>-d<alpha phi>.

The derivative forcing proportional to t d in the second line of(3)
drops out of the whole-law current because E[d c_X]=0. But r_b does
NOT drop out. In the scalar basis cosY,cosZ,sinY,sinZ its complete
current/source covariance is

    -(1/5) [[B²,AB,0,0],[AB,A²,0,0],
             [0,0,B²,0],[0,0,0,A²]].                 (6)

Equivalently it is the rank-three tensor
-(psi tensor psi+B² sinY tensor sinY+A² sinZ tensor sinZ)/5.
The exact transport contribution is consequently

    -(t/5)<psi²>
    -(1/5) integral_0^t [B² C_Y(s)+A² C_Z(s)]ds,
    C_Y(s)=<sinY exp(-sT)sinY>, C_Z(s)=<sinZ exp(-sT)sinZ>. (7)

An arbitrary correlated homogeneous axial return contributes its actual
bounded analytic transport output P0(t). The conserved e-weighted psi
moment supplies a constant. Hence, with E=A²+B², the complete physical
mean acceleration coefficient of the V column is

    R_V(t)=E t/15-<alpha phi_bar(t)>
           -(1/5) integral_0^t[B² C_Y(s)+A² C_Z(s)]ds
           +P0(t)+C0.                              (8)

No orientation frequency has been averaged to derive(8). The mean and
transport rows are observed on the actual Euler solution with the same
initial input. In particular(8) contains a genuine continuum of transport
responses until an actual cancellation or finite physical-row closure is
constructed.

## 4. Phase, energy and physical finite-row alternatives

For zero initial material first cell, the same-D,V initial phase remains
rho J. An O(k) microscopic velocity return r_V changes the actual VV
energy by rho k² E||r_V||²/2, and its cross row with the D material rate
b_D is rho k² E<b_D,r_V>. These are retained independently of whether
the observed current(8) simplifies. A full coherent TR choice may cancel
an odd cross row, but it does not erase the even gradient mass.

If a configuration return xi is used, the initial one-form is0180's

    Theta0=rho V.dD+rho k²[-P_kappa<u(kappa.xi)>.dD
                 +<v+u.grad xi+(kappa.u)D, dxi>].

The actual physical mean velocity condition and initial G/tag transport
must accompany this form. A mean-zero return alone does not establish
a physical canonical angle or a matching gradient mass.

The sufficient full-field Jordan construction sets Lr_V=z_D-b_V.
A weaker isolated acoustic construction sets the whole-law output of y
to zero. A coupled construction instead realizes its output through the
actual finite material angle, measured spin and integrated-current rows.
The parent permits this last possibility; autonomous isolated mean
evolution is not an additional acceptance requirement. Each candidate
has to retain(1)–(8), its initial phase and complete energy at the scope
where it is used.
