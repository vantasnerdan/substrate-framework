# Exact same-carrier accessible packet and semigroup lift

## 1. Whole-space group and energy estimate

Let `X=L^2_sigma(R^3)` and let `P_L` be the Fourier Leray projector. For the
fixed smooth compact Gavrilov equilibrium, write

    Gv=-P_L[(u_* dot grad)v+(v dot grad)u_*]=A_0 v+Bv.       (1)

On its maximal transport domain, `A_0=-P_L(u_* dot grad)` is skew-adjoint;
`B=-P_L(grad(u_*) .)` is bounded on `X`. Hence `G` generates a group and

    ||S(t)||_{X->X} <= exp(|t| ||grad u_*||_infinity).       (2)

This also supplies the Duhamel constant used below. No global `H^1` domain is
asserted.

## 2. Every unstable polarization is dynamically accessible

Fix the 0032 resonant ray `(x_*(t),k_*(t))`. Its invariant scalar
`mu=k_* dot omega_*` is nonzero on the selected sufficiently small shell (the
0032 local expansion has a nonzero leading term). If `A` is any physical
polarization with `k dot A=0`, choose the leading displacement amplitude

    a=(k cross A)/(k dot omega_*).                            (3)

Then `k dot a=0`.  The action angle `phi` is periodic, so choose an integer
`N>=1`, set `h=N^-1`, and use the globally single-valued exact orbit tangent

    delta omega=curl(xi cross omega_*),
    xi_N=N^-1 chi_delta Re[a exp(i N phi)].                 (4)

is

    i[(k dot omega_*)a-(k dot a)omega_*].                    (5)

Biot--Savart has symbol `i k cross(.)/|k|^2`; substituting (3) into
(5) gives exactly `A` (up to the common harmless real/complex phase
convention). Thus the expanding eigenpolarization `e_+` of `M_*` belongs to
the principal symbol of genuine dynamically accessible data. Apply a
compactly supported solenoidal correction to `xi_h`; it is one order lower
because `k dot a=0`, and (4) remains an exact coadjoint tangent.

The raw reconstructed velocity `v_N=BiotSavart(delta omega_N)` is `O(N^-1)`
in `L^2`. Define

    q_N=v_N/||v_N||_X.                                       (6)

Scalar normalization preserves accessibility. The leading symbol is nonzero,
so `c_delta/N <= ||v_N||_X <= C_delta/N`. Smooth compact vorticity makes
`q_N` a smooth finite-energy member of the maximal domain and

    ||q_N||_D=||q_N||_X+||Gq_N||_X <= C_delta N.             (7)

Equation (7) records graph size; it is not a uniform graph-norm claim.

## 3. Direct finite-time Euler--Egorov remainder

Let `phi_t=phi_0 composed Phi_{-t}`. Then

    D_t phi_t=0,     k_dot=-(grad u_*)^T k.                  (8)

For the velocity ansatz obtained from (4), apply `P_L` exactly and solve the
leading transverse transport

    D_t A=-(grad u_*)A
      +2 k [k dot (grad u_*)A]/|k|^2,   k dot A=0.           (9)

The second term is the principal symbol of the whole-space pressure
projection. A first corrector cancels the longitudinal and order-zero
residual. Standard symbol composition for the constant-coefficient multiplier
`P_L` then gives, on every fixed interval `[0,T]`,

    sup_{0<=t<=T} ||r_N(t)||_X
       <= C_T N^-2 ||chi_delta a||_{H^3},                   (10)

for the unnormalized accessible velocity, whose norm is `Theta(N^-1)`. Three
profile derivatives suffice: two for the order-zero multiplier expansion and
one for the solenoidal/accessibility corrector. For a tube cutoff of width
`delta`,

    ||chi_delta a||_{H^3} <= C delta^-3 ||chi_delta a||_2.  (11)

After (6), Duhamel and (2) therefore give the relative estimate

    sup_{0<=t<=T} ||S(t)q_N-q_N^WKB(t)||_X
      <= N^-1 C_T delta^-3 T exp(T||grad u_*||_infinity).   (12)

All constants are independent of `N` after `T,delta` are fixed.

The oscillatory source in (4) is supported in an invariant action tube a
positive distance from the flat carrier collars and chart boundary. The
velocity is not compactly supported. The local part of `P_L` is already in
(9). Between the tube and a separated collar/exterior cutoff, its kernel and
the Biot--Savart kernel are smooth; repeated integration by parts with
`|d phi_t|>0` gives `C_{M,T,delta}N^-M` for every fixed integer `M`.
Commutators with
the tube cutoff are order minus one and are included in (10)--(12). Hence the
Hodge tail, collar, and exterior are retained rather than set to zero.

The source theorem audited in `access-inventory.md` proves the same BAS/WKB
mechanism on `T^3`; (8)--(12) are the direct `R^3` fixed-time bridge needed
here and do not import its torus essential-spectrum theorem.

## 4. Every circuit count

Let `lambda_+>1`, `e_+`, and `T_*` be fixed by 0032. For any integer `j>=1`,
put `t_j=jT_*`. Smooth finite-time dependence of the base flow, cotangent
lift, physical frame, density, and cocycle implies that for every `r in (0,1)`
there is an invariant tube width `delta_j>0` on which the WKB output/input
energy ratio differs from the central value `lambda_+^j` by at most
`r lambda_+^j/4`. The axisymmetric rotation after each meridional circuit is
an isometry and returns the carrier, frame convention, and density; derivative
resonance returns the covector.

With `T=t_j` and `delta=delta_j` in (12), choose an integer

    N_j>4 C_tj t_j exp(t_j||grad u_*||_infinity)
             /[r lambda_+^j delta_j^3].                     (13)

Then the exact normalized accessible solution satisfies

    ||S(t_j)q_Nj||_X >= (1-r/2)lambda_+^j.                  (14)

This is the frozen quantifier

    for every j, choose delta_j, then choose integer N_j.    (15)

No constant in a single WKB expansion is claimed uniform as `j` grows.
Taking the supremum over unit inputs in (14), and then `r downarrow 0`, yields

    ||S(t_j)||_{X_DA->X_DA} >= lambda_+^j   for every j.    (16)

Therefore

    limsup_{j->infinity} (1/t_j)log||S(t_j)|| >=
      (log lambda_+)/T_*=gamma_*>0.                         (17)

## 5. Weak packets, quotient, and essential norm

For a fixed `j`, take a diagonal sequence `delta_l downarrow 0` and then
integers `N_l` with `N_l delta_l^3 -> infinity` satisfying (13). The normalized packets concentrate
on the one-dimensional orbit and oscillate with nonzero covector, so
`q_l weakly ->0` in `X`. Their exact images have the same property at `t_j`.
For every compact operator `K` on the accessible closure,

    K q_l ->0 strongly,
    ||S(t_j)-K|| >= limsup_l ||S(t_j)q_l|| >= lambda_+^j.    (18)

Thus

    ||S(t_j)||_essential,X_DA >= lambda_+^j.                (19)

Euclidean symmetry tangents span a finite-dimensional smooth space. Weak
nullness makes their projections vanish at input and output, so (16)--(19)
hold unchanged on `X_DA/sym` when its norm is quotient distance. This proves
the time-domain version of preregistered Route B. It does not prove the
generator Weyl statement `(G-z)q_N ->0`; longitudinal closure around the
orbit would be an additional construction and is unnecessary for (17),(19).

The exact Euler evolution maps smooth coadjoint tangents to smooth coadjoint
tangents by equivariance of vorticity transport. Boundedness of `S(t)` extends
this invariance to `X_DA`; hence all norms above are norms of the restricted
semigroup rather than ambient approximations.

## 6. Verdict and physical separation

Route A is established with the per-circuit diagonal quantifier (15). Route
B's time-domain weak-packet/essential-norm alternative is also established by
(18)--(19); the generator pseudomode subroute remains unneeded and unclaimed.
The corrected hyperbolic cocycle therefore produces genuine exponential
growth of the linear Euler semigroup norm on normalized finite-energy
dynamically accessible data of this one fixed compact carrier.

This does not prove nonlinear Lyapunov instability. Such a claim additionally
requires a differentiable nonlinear Euler flow map in a chosen topology and a
quadratic remainder controlled on a logarithmic/growth-compatible time scale.
It also does not prove carrier breakup or an LP2 stable-particle embedding.
Indeed, the result is adverse evidence for stability of this particular
Gavrilov carrier sector and activates a different-carrier/restoring-sector
route. No action quantization, exchange, Born, relativity, electron, or
neutrino mechanism follows from the linear growth theorem.
