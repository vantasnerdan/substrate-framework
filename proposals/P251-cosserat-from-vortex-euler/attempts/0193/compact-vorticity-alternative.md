# A derived generalized-force-free compact-vorticity alternative

This executes the second radial-localization representation registered
in 0193. It is not an assertion that a scalar cutoff preserves the
force-free property. In fact `f h` is force-free only if
`(d f)'=0`; a nonzero compact radial velocity cutoff cannot satisfy
that equation. Here the axial velocity is changed by its actual Euler
Bernoulli equation instead, with the original core preserved exactly.

## Exact stationary profile

Keep the same smooth cutoff and set

    V=C r/d, W=C c/d, V_R=chi_cut(r/R)V,
    p_R(r)=-C^2/(2c^2)+integral_0^r V_R(s)^2/s ds,
    W_R=sqrt(-2p_R-V_R^2),
    u_R=V_R e_theta+W_R e_z.                            (1)

The square root is real and strictly positive, not a formal branch
assumption. Since the uncut pressure is `p=-C^2/(2d)`, direct
subtraction gives

    W_R^2=W^2+(1-chi_cut^2)V^2
                +2 integral_0^r (1-chi_cut(s/R)^2)V(s)^2/s ds. (2)

Every additional term is nonnegative. The positive root is smooth
for each finite R, including the axis and the cutoff edges. Inside
`r<R` the profile is EXACTLY the original helical field. Outside
`2R`, `V_R=0` and `W_R=W_infinity,R>0` is constant. Thus its vorticity
is compactly supported, but its velocity has a uniform axial exterior.

Write `Z_R=V_R'+V_R/r`. Differentiating (1) gives

    W_R W_R'=-V_R Z_R,
    curl u_R=(0,-W_R',Z_R)=(Z_R/W_R)u_R.                (3)

Also `div u_R=0`, its only convective acceleration is `-V_R^2/r e_r`,
and `p_R'=V_R^2/r`. Hence this is exact smooth stationary Euler and
generalized force-free, with constant Bernoulli `p_R+|u_R|^2/2=0`.
The factor `Z_R/W_R` can change sign in the return annulus; it is
not a constant Beltrami eigenvalue and supplies no EPS Runge license.

## Full fixed-sector operator convergence

For an arbitrary smooth axisymmetric column `u=r O(r)e_theta+W(r)e_z`,
the complete fixed-sector linear Euler generator is

    L v=-P_mk[i(mO+kW)v+2O Jv
                         +(rO' e_theta+W'e_z)v_r].     (4)

In (1), `O_R=chi_cut f`, so the angular differences in (4) are bounded
by `C_R-independent/R^2`, as already derived in the primary proof.
Equation (2) gives `||W_R-W||_infty<=C_1/R`.

The derivative bound needs more than pointwise square-root convergence.
For `r>=R/2`, direct differentiation of
`g=W_R^2=-2p_R-V_R^2` gives, with fixed profile constants,

    |g'|<=C_2/R^3, |g''|<=C_3/R^4.

Choose `K=max(C_3,2C_2)/R^4`. At any `r>=R`, the Taylor displacement
`h=-g'(r)/K` lies within `|h|<=R/2`. Nonnegativity of g throughout
that interval and the upper Taylor bound imply

    0<=g(r+h)<=g(r)-g'(r)^2/(2K),
    |W_R'|=|g'|/(2sqrt(g))<=sqrt(K/2).                 (5)

There is no division by a guessed uniform positive lower bound on
W_R near the inner edge. In `r<R`, W_R=W exactly. Thus
`||W_R'-W'||_infty=O(R^-2)`. Substituting the bounds into the FULL
operator (4), including the common norm-one Leray map, proves

    ||L_R-L||<=C_4 |k|/R+C_5(|m|+1)/R^2.              (6)

For fixed `(m,k0)`, this is enough for the very same rank-one Riesz
continuation as the primary compact-velocity route. It retains the
outer axial Doppler shift and the changed helical-charge dynamics.
It is not a statement about a pure constant-charge quotient outside
the unchanged core.

## Actual Kelvin and energy continuation

For the resulting eigenmode write `lambda_R=i sigma_R`, initially
allowing a complex sigma_R, and put

    nu(r)=sigma_R+m O_R(r)+k0 W_R(r),
    g_vec=r O_R' e_theta+W_R' e_z,
    xi_R=v_R/(i nu)+g_vec v_R,r/(i nu)^2.               (7)

By (6), nu is uniformly separated from zero for sufficiently large R:
the uncut value is the positive constant sigma. The varying-frequency
derivative `nu'=m O_R'+k0 W_R'` cancels the divergence of the second
term, and the formula satisfies full Lin reconstruction.

For clarity, the extra radial part of the Kelvin identity is not
discarded. Let `Z=2O+rO'`. The tangential and axial Euler equations
give

    i nu[(rO')v_theta+W'v_z]
                      +[(rO')Z+W'^2]v_r=-i nu' pi.

Consequently the radial correction in `xi cross curl u` is precisely
the derivative of the denominator in `grad[pi/(i nu)]`. Together with
the other components this gives the exact identity

    xi_R cross curl u_R=v_R+grad[pi_R/(i nu)],
    P(xi_R cross curl u_R)=v_R.                        (8)

The full KKS form is continuous by (6),(7) and the bounded vorticity
differences. Its nonzero sign again forces the continued pole to stay
on the imaginary axis; `sigma_R>0` and its physical-clock phase
Hamiltonian `h_R=sigma_R beta_R` remains positive. The uniform exterior
has zero vorticity; there the Doppler frequency is
`sigma_R+k0 W_infinity,R`, retained and nonzero. The actual pressure
tail is the decaying radial Helmholtz branch, not a wall condition.

The background has finite energy DIFFERENCE from its specified uniform
axial exterior per axial period; it does not have finite absolute
background energy on R2 times that period. The subtraction is fixed by
the actual exterior velocity, not used to select the mode's energy sign.
The positive mode action is the full Euler phase Hessian in its original
physical laboratory clock. No Galilean clock is silently substituted.

## Exact physical tag/current transfer and boundary

In the unchanged inner core, nu=sigma_R and the local actual mode obeys
the same radial equation (8) of the primary proof. All its local rows
converge on the fixed tag annulus. The same three fixed controls therefore
solve that proof's actual system (9), with this alternative's full h_R,
phi_R and sigma_R. The positive stationary tag then has, exactly on
its entire linear eigenmode history,

    G_z=eta M_R theta, S_z=eta M_R theta_t,
    S_z=-c delta P_z.

The invariant helical scalar is needed only on the tag support, where
the base is exactly the unchanged original field. It is not claimed
stationary in the outer annulus of this more general column. The full
material variation uses the particular displacement (7), including all
outer reactions in the action and all inner position terms in the spin.

This alternative is established as a genuine fixed-sector pole/current
continuation on an explicit compact-vorticity generalized-force-free
steady column. Its advantages and costs differ from the primary route:
it preserves force-free structure with a variable factor but requires
a uniform axial exterior; the primary route has compact velocity and
finite absolute energy per axial period. Neither gives a constant-curl
EPS tube, a full-channel spectral gap, or an array/bending theorem.
