# Flat-tail generator control beyond the physical velocity estimate

This supplements the actual residual in `curved-residual.md`. The
profile and derivative order are fixed before the small positive
spectral parameter c. The assertion here is about the reflected
straight threshold channel, not an assumed curved inverse.

## 1. Exponential flatness costs logarithms, not inverse powers

At the actual 0211 edge use t=delta/phi. The explicit cutoff and its
primitive give

    w(phi)=C phi² exp(−delta/phi)(1+O(phi)),
    |partial_phi^j w| <= C_j w(1+t)^(2j).                  (1)

The same bounds hold for finitely many physical radial derivatives,
because phi has a simple zero with bounded higher derivatives. Let
0<c<c_* and L=1+|log c|, with fixed dimensional reference units inside
the logarithm. For every fixed nonnegative integer M,

    sup_phi [w/(c+w)](1+t)^M <= C_M L^M.                  (2)

For t<=2L+constant this is immediate. In the remaining region,
w/c<=C exp(L−t) times a fixed polynomial in t; exponential decay
dominates t^M, giving the same bound. Changing the positive fixed
reference c_* only changes the constant. Consequently

    sup |partial^j w|/(c+w) <= C_j L^(2j),
    sup w |partial^j [1/(c+w)]| <= C_j L^(2j).            (3)

The second formula follows by the finite Faà di Bruno expansion:
each differentiated denominator is a product of derivatives of w,
divided by an equal number of additional factors c+w. Every such
ratio is covered by (1)-(2). The exterior w=0 contributes zero to
these weighted rows. There is no appeal to smoothness alone.

## 2. The actual forced reconstruction

Use the complete formulas of 0206 with d=c+w, y=s f+i s eta_s/(kd):

    [(d² y'/s)]'+(kappa²−k²d²)y/s
       = Z eta_theta+d(eta_z'−ik eta_s),
    f=y/s−i eta_s/(kd),
    v_s=−ik d y/s−eta_s,
    v_theta=−Z y/s,
    v_z=(d y)'/s−eta_z.                                  (4)

The weighted form Q_c controls d y' and d y/s, with the regular-axis
weights recorded in 0206, uniformly in c down to its positive-sector
edge. The physical v_z control additionally uses (12) of
`curved-residual.md`. The source norm is exactly (13) there.

On a fixed outer annulus away from the axis, set b=d y. Multiplying
the differential equation by s/d gives the particularly useful exact
normal form

    b''−b'/s−[d''/d−d'/(sd)+k²−kappa²/d²] b
       = (s Z/d)eta_theta+s(eta_z'−ik eta_s).              (5)

No first derivative with a singular d'/d coefficient remains. The
coefficients of (5) and their fixed derivatives are bounded by powers
of L: kappa²=−2ww'/s, Z=−ww'/(s Omega), and (1)-(3) apply. Thus the
one-dimensional equation, its derivatives, and the Q_c bound control
b in every fixed annular H^j norm by a finite power of L times fixed
source Sobolev norms. At the axis w is strictly positive, so regular
axis estimates supply the unchanged interior bound without a flat
weight. The matching exterior is the actual decaying pressure return,
not an imposed wall.

The forced particular term in f contains 1/k. Equations (3)-(5) imply
for fixed j and 0<k<=k_* a bound of the form

    ||w f||_{H^j(core)}
       <= (C_j/k) L^(N_j) ||eta||_{H^(j+2)(core)}.         (6)

Here multiplication and derivatives are taken as displayed; additional
finite background derivative rows are estimated by (1)-(3). Formula
(6) is a weighted radial displacement estimate. It is not an
unweighted generator estimate at the boundary, nor a claim that every
other displacement component has already been reconstructed.

For c comparable to log R/R and k=n/R, the actual source from (10) of
`curved-residual.md` therefore produces a radial weighted displacement
O(R^-1 log^(N_j+1) R). This is a useful response scale on the intended
positive bending contour. Its physical velocity remains at the
stronger O(log R/R²) source scale from 0206.

## 3. Exact full coadjoint generator, including the cancellation

The f in (4) is a forced radial Lin variable. It need not be the radial
component of a generator whose coadjoint velocity is the response.
The latter is constructed directly, without identifying these two:

    Xi_s=y/s,
    Xi_z=i y'/(ks),
    Xi_theta=i[eta_theta−2Omega y/s]/(kd).                 (7)

This is exactly divergence-free. For omega=(0,w',Z), set F=Xi cross
omega and

    pi=i[d y'/s−eta_z]/k.

Substitution of the forced Sturm equation gives

    F−grad_k pi=v,                                        (8)

where v is ALL THREE components in (4). Thus VXi=v with the actual
full pressure return. On the source this proves coadjoint-range
membership; a curl-potential completion of Xi outside that source
does not change its induced velocity.

The large Xi_theta is multiplied by Z in F, not by w. The exact
force-free identity Z=−ww'/(s Omega) gives

    F_theta=−Z y/s,
    F_z=w' y/s,
    F_s= i Z eta_theta/(kd)
          −i 2Omega Z y/(ksd)−i w' y'/(ks).              (9)

In terms of b=d y, every coefficient in (9) is bounded by a finite
power of L, apart from the explicit 1/k: for example
2Omega Z/d²=−2ww'/(s d²), and w'/d has the bound (3).
Consequently the complete compact coadjoint-force representative has

    ||Xi cross omega||_{H^j(core)}
       <= (C_j/k) L^(M_j) ||eta||_{H^(j+3)(core)}.         (10)

This includes every velocity-producing generator component and its
actual pressure; it does not assert that Xi itself is uniformly
bounded where the vorticity vanishes. The source factors extend
smoothly by zero. The compact-force representation (10), rather than
an unweighted Xi norm, is the appropriate next pressure-Schur input.

## 4. What this licenses for the remaining calculation

The first geometric correction has poloidal Fourier order one. Its
zero-average diagonal on the m=0 subspace is an algebraic selection
rule; the full actual pressure operator must be expanded in the same
flux-coordinate convention before using that rule as an operator
estimate. A two-step correction through a uniformly invertible
nonzero-poloidal complement would then have size R^-2 rather than
R^-1. Inserting (6) suggests the small factor
R^-2 k^-1 polylog(1/c), which tends to zero at the bending scaling.

This last sentence specifies the next concrete Schur construction.
It does not silently assert the whole curved graph estimate or the
nonzero-poloidal inverse. The exact
additional objects are the complementary velocity/pressure rows and
their bilinear source estimates on the same actual background.
