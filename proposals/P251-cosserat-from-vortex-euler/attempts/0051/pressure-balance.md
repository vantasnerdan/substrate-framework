# Exact pressure-mediated force and couple balance

Let D_a(t) form a material partition of the fluid, with constant density rho,
mass M_a, center X_a, velocity V_a=Xdot_a and intrinsic angular momentum

    L_a=integral_Da rho (y-X_a) cross (u(y)-V_a) dy.

Euler and Reynolds transport, with no cell-rigidity assumption, give

    M_a Vdot_a=-integral_boundary(Da) p n_a,
    Ldot_a=-integral_boundary(Da) (y-X_a) cross (p n_a).

For each unordered neighboring pair (a,b), orient the force on a as
t_ab=-p n_a on the common face and define

    R=X_b-X_a,
    F=integral_face t_ab,
    m_a=integral_face (y-X_a) cross t_ab,
    m_b=-integral_face (y-X_b) cross t_ab=-m_a+R cross F.

This is exact pressure reaction. The opposite faces do not supply opposite
intrinsic torques, because their centers differ. In particular pressure
need not produce a central resultant F parallel to R for irregular cells.

Write delta_a=delta(x-X_a), x_s=X_a+sR, and delta_s=delta(x-x_s). Define
the pair's force stress and couple stress as distributions

    sigma_ab=integral_0^1 F tensor R delta_s ds,
    mu_ab=integral_0^1 [m_a-s(R cross F)] tensor R delta_s ds.

The first tensor index is the force/couple component; divergence acts on
the second. The identity R.grad_x delta_s=-partial_s delta_s yields

    div sigma_ab=F(delta_a-delta_b),
    div mu_ab=m_a delta_a+m_b delta_b-(R cross F) integral_0^1 delta_s ds.

Let ax(sigma)_i=epsilon_ijk sigma_jk. Then
ax(F tensor R)=F cross R=-(R cross F), and therefore

    m_a delta_a+m_b delta_b = div mu_ab-ax(sigma_ab).

The sign is fixed by the stated index convention; transposing stress or
using the opposite axial convention changes both occurrences together.
It is not chosen by comparison with a desired continuum equation.

Sum pairs, including explicit external boundary tractions where present.
Introduce the exact center distributions

    P=sum_a M_a V_a delta_a, S=sum_a L_a delta_a,
    T=sum_a M_a V_a tensor V_a delta_a,
    C=sum_a L_a tensor V_a delta_a.

Their balances are

    partial_t P+div T=div sigma,
    partial_t S+div C=div mu-ax(sigma).

The same identities hold after convolution by any smooth spatial averaging
kernel, since convolution commutes with these distributional derivatives.
No local constitutive approximation was used to obtain them. The cell
labels and their material centers remain part of the coarse-graining.

Finally x cross div sigma=div(x cross sigma)+ax(sigma) in these conventions.
Adding the spin equation cancels the antisymmetric force-stress term and
gives total angular-momentum conservation exactly. Thus nonzero coarse
couple stress and asymmetric center-based stress are compatible with the
symmetric microscopic pressure stress; they record its resolved moment
about different material centers.

The independent rotational coordinate and its inertia still come from the
constructed Euler microscopic action, not these balance identities alone.
For reference the relation between vorticity impulse and material angular
momentum is, on a finite domain D centered at X,

    integral_D rho r cross u
      =-rho/2 integral_D r² omega
        +rho/2 integral_boundary(D) r² (n cross u),  r=y-X.

It follows by integrating curl(r² u)=2 r cross u+r² omega. The boundary
term is generally nonzero and is retained when matching a relative-EPS
orbit's impulse to a material-cell spin. Dropping it is not licensed by
passing the bond-balance checks.
