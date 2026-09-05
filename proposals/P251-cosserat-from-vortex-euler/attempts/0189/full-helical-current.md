# One positive stationary helical tag with its complete actual current

Use the ACTUAL m=2 ground radial Euler mode of0185, not a prescribed
oscillator. Its smooth positive radial streamfunction is phi(r), its
frequency sigma>0, and its full action coefficient h_m>0. Put
d=c²+r², f=C/d, and use the physical phase convention

    phi_mode=phi(r)[a cos(ms)+b sin(ms)],
    a_t=sigma b, b_t=-sigma a, s=theta-z/c.

The exact Lin displacement in0185 gives, in cylindrical components,

    xi_r=m phi/(sigma r)[a cos(ms)+b sin(ms)],
    xi_theta=[c²phi'/(sigma d)-m f_r phi/sigma²]
                                             [b cos(ms)-a sin(ms)]. (1)

The Cartesian basis rotation is already in this displacement. In
particular v_theta=xi_theta,t-f_r r xi_r; neither of the two terms may
be replaced by a guessed time derivative of material spin.

## Literal displacement and mechanical spin on the same tag

Take the stationary material fraction

    w0=chi(r)[1+epsilon b_tag(r)cos(ms)]

of0185, with smooth compact radial chi>=0 and sufficiently small epsilon.
Its nonzero reference quadrupole fixes a unit physical angle

    theta=c_obs b, c_obs=I_obs/(sigma epsilon B_tag),
    B_tag=integral chi b_tag r^(m+1)dr,
    I_obs=integral chi r^(m-1)(m phi+r phi')dr !=0.

The actual scalar phase mass is M=h_m/(sigma² c_obs²)>0. The m=2
reference centroid and first planar momentum vanish, so the following
are centered physical moments with no omitted centroid correction:

    G_z=rho integral w0(r cross xi)_z =G_b b,
    G_b=rho pi Lz epsilon integral chi b_tag
           [c²r²phi'/(sigma d)-m f_r r²phi/sigma²]dr,

    S_z=rho integral w0[r v_theta+(r² f)'xi_r]=S_a a,
    S_a=rho pi Lz epsilon integral chi b_tag J(r),
    J(r)=-c²r²phi'/d+2mCc²r phi/(sigma d²).             (2)

Both the moved position and the velocity variation enter S. The
linked axial mechanical momentum is still S_z=-c delta P_z as in0185;
it is not discarded to isolate an unphysical rotor.

Equations(1),(2) derive the exact connection

    S_z-G_z,t=(2mC rho pi Lz epsilon/sigma) a R_tag,
    R_tag=integral chi b_tag R(r)dr,
    R(r)=r phi/d.                                      (3)

Thus matching the mechanical spin alone does not automatically match
the accumulated physical current. Equation(3) names the additional
moment, rather than adding an unexplained condition to the old result.

## Three independent radial rows give both equalities exactly

Fix any eta>0 and any nonzero reference moment B_*. Impose

    B_tag=B_*, R_tag=0,
    integral chi b_tag J(r)dr
                       =-eta h_m B_*/(rho pi Lz I_obs). (4)

These are THREE linear conditions on the signed marker modulation,
not on the sign of the material mass fraction. The functions
r^(m+1), R(r), J(r) are linearly independent on every open core
interval. Indeed analyticity extends any putative constant linear
relation to all positive r. Exponential decay of the actual eigenmode
forces the polynomial row's coefficient to vanish. If J=k R, its
pointwise definition would give

    phi'/phi=2mC/(sigma r d)-k/(c²r),

whose nonzero solutions have algebraic, not exponential, decay.
This contradicts0185's actual endpoint asymptotic. R is nonzero since
the ground mode is positive. This proves independence without measuring
a nearly singular moment matrix or assuming three generic constraints.

Choose three evaluation radii in a region where chi>0 with nonsingular
row matrix; independence guarantees such a choice. Smooth narrow radial
bumps preserve its determinant, and the exact matrix inverse gives a
smooth bounded b_tag satisfying(4), with support in that core interval.
Choose epsilon>0 afterward so |epsilon b_tag|<1 and chi<=1/2. Then
0<=w0<=1, its nonzero quadrupole persists, and w0 is exactly stationary
because every factor is a function of the material invariants r and s.
The mode, its frequency and its full action have not been adjusted to
fit the marker. The physical observation is the constructed positive tag.

Substitution of(4) into(2),(3) gives on the entire linear mode history

    S_z=eta M theta_t, G_z=eta M theta,
    G_z(0)+integral_0^t S_z=eta M theta(t).             (5)

In particular the initial displacement moment is the ACTUAL Lin one;
no integration constant or independent rotor is appended to force(5).
This supplies precisely the nonzero physical optical-current interface
used by0182, at the exact helical carrier and the stated periodic mode
scope. Eta may equal1 or any fixed positive number.

route_verdict: established as stated.
evidence_scope: actual smooth helical Euler mode with positive action,
one positive stationary material quadrupole and complete all-time
displacement–spin current. The mode action is finite per axial period.
Off-helical carrier-two current, finite-bandwidth packet errors, same-
field acoustic response and stationary Euclidean EPS embedding remain
separate active constructions. Equation(5) is not asserted for an
arbitrary perturbation of this tag or for an arbitrary carrier family.
