# Actual common-laboratory spatial jets, pressure return and density

Fix the finite ring and its preparations from `finite-ring-transfer.md`
before choosing an array spacing. In this file P denotes the cubic
period, not the pressure or Leray projector. The radius and every
condition number of the chosen ring are fixed constants at this step.

## 1. High harmonics give an actual exterior multipole estimate

The stationary compact ring is axisymmetric in its global toroidal
angle. Its full Euler/Lin operator preserves angular harmonics in
cylindrical components, including the pressure equation. The initial
displacement has only |n|>=N, N>=12. Its vorticity perturbation is
compact in the fixed ring support. The linear vorticity equation

    eta_t+u dot grad eta-eta dot grad u
                         =omega dot grad v-v dot grad omega

keeps that support: the background is tangent to its cutoff surfaces,
and the remaining sources vanish off the fixed vorticity region.
Outside a ball containing the ring, the actual velocity is a gradient
of a decaying harmonic scalar. This holds initially by Kelvin
preparation and afterwards by the full Euler pressure equation.

A scalar spherical harmonic carrying toroidal order n has degree
l>=|n|. Thus its exterior harmonic expansion and its physical gradient
give, for every fixed time and carrier derivative inventory,

    |partial_x^a partial_p^j partial_t^b v_R(x)|
      +|partial_x^a partial_p^j partial_t^b Xi_R(x)|
                <= C_(R,T,a,j,b) (1+|x|)^(-N-2-|a|).         (1)

The second term follows from Lin outside the same large ball, where
u=0 and Xi_t=v; the initial displacement is compact there. This is
a consequence of actual angular symmetry and exterior pressure, not
a claim that pressure has compact support. Constants are finite for
the already fixed R; no limit in R is taken using this constant.

The same angular selection shows that the complete initial displacement
and Kelvin-force Cartesian moments of degree <=N-2 vanish. Vector
components shift angular order by at most one; multiplying by a
polynomial of degree d shifts by at most d. This statement concerns
FULL Euclidean moments, not the selected tag moments.

## 2. A real common-K preparation on the exact stationary array

Take P sufficiently large that the entire ring and initial displacement
supports in different cells are disjoint. The background

    u_P(x)=sum_(a in P Z3) u_R(x-a)                         (2)

is exactly stationary smooth Euler, with the sum of the actual compact
pressures, shifted to their exterior constants. Its finite derivative
norms are bounded by the single-ring norms, independently of P.

For a common laboratory K, set p(K)=p0+t0 dot K, where t0 is the
declared axis of the observed arc, and prepare the physical field

    Xi^K(0,x)=sum_a exp(i K dot a) Xi_R(p(K),x-a).           (3)

The periodic Bloch amplitude is exp(-i K dot x) times this field.
It is solenoidal for grad+iK exactly. Formula (3) explicitly cancels
the within-cell transverse phase; its only selected internal carrier
variation is p(K), while neighboring cells have their genuine common
laboratory phase. It is an actual initial preparation, not an assumed
axial embedding of a dispersion curve. Real histories use the complex
conjugate -K preparation as well.

Because the initial displacement supports are disjoint, the Kelvin
force is exactly the same phased sum of single-ring Kelvin forces.
Periodization commutes with the full Leray multiplier. Its apparent
zero-frequency ambiguity vanishes for these data by the moment
selection above. Thus the array's actual initial velocity equals the
periodization of the actual isolated-ring initial velocity. The same
is true for all of its existing initial spatial jets.

Define the reference history by periodizing the ACTUAL isolated-ring
solutions in (3), not just their initial eigenvectors. Estimates (1)
make the phased sums and their first two K derivatives converge
absolutely in local Sobolev spaces; a K derivative costs at most a
lattice distance or an actual p derivative. The reference solves the
array Euler/Lin equations except for the explicit cross-copy terms

    sum_(a!=b) [u_a dot grad v_b+v_b dot grad u_a],
    sum_(a!=b) [u_a dot grad Xi_b-Xi_b dot grad u_a].        (4)

These terms are supported where a background copy is nonzero. Their
full projected pressure is retained. In a cell their H^s norms satisfy

    ||partial_k^j residual||_Hs <= C_(R,T,s) P^(-N-2+j),
                     K=k kappa, |kappa|=1, 0<=j<=2,        (5)

uniformly in kappa near k=0. The convergent lattice sum is bounded
by sum_(a/P!=0)|a/P|^(-N-2+j). No cancellation of a conditionally
convergent dipolar image series is required to obtain (5).

## 3. Full Bloch zero harmonic and actual derivative estimates

For nonzero lattice wavevector q, the Bloch projector symbol is

    I-(q+k kappa)(q+k kappa)^T/|q+k kappa|^2.

Its first and second k derivatives have operator norms <=C P and
<=C P^2 on |k|<pi/P. At the zero harmonic retain the physical
projector P_kappa=I-kappa kappa^T for k!=0, with its continuous
fixed-ray extension at k=0. It is NOT replaced by zero. The initial
zero mode vanishes to high order by the moment selection; the exact
zeroth-order Euler mean is conserved and zero. Hence this extension
coincides with the actual k=0 history, although differentiated mean
pressure and current rows need not vanish.

The projected Euler energy identity uses only the uniform coefficient
norm of (2), integration by parts on the cell, and norm(P_k)=1.
Differentiate the exact comparison equation on each fixed ray, using
the displayed bounds and (5). Increase the spatial inventory when a
coefficient derivative acts on a differentiated lower-order solution.
Duhamel and the Lin transport estimate give the conservative bound

    max_(j<=2,b<=3) ||partial_k^j partial_t^b
                [actual-reference]||_Hs(cell)
                    <= C_(R,T,s) P^(4-N).                  (6)

The initial difference is exactly zero, not just small. The loose
four-power allowance covers projector, gauge and spatial-derivative
bookkeeping; (5) gives sharper powers when desired. N>=12 is ample.
The same argument gives a third-order remainder with finite constants
at this FIXED P, uniformly in kappa. This is a genuine finite-window
derivative theorem, not a fixed-K estimate silently differentiated.

Individual realizations have well-defined ray jets; the zero harmonic
need not have a Cartesian second derivative before averaging. Rotate
or reflect the WHOLE array, preparations, tag and K together, with the
actual time-reversal phase embedding when used. A parity-paired
axial-vector response has no odd spatial row. Rotational covariance
forces its homogeneous second-order averaged vector tensor to be
alpha |K|^2 I+beta K K^T. The uniform ray remainder therefore supplies
the actual Cartesian second jet of that averaged response. Identical
reasoning applies to its averaged phase/action blocks at their tensor
ranks. This justifies the common-lab isotropic second jet without
discarding a realization's mean-pressure term.

## 4. Both full forms and actual physical normalization

Compute the Bloch KKS and inherited Euler Hessian per cell from the
actual initial data, using conjugate pairing of K and -K. The KKS
is local in omega and Xi; disjoint initial displacement supports
make it exactly the single-ring KKS at p(K), with all its derivatives.
The Hessian includes the actual projected velocity and pressure return.
Its difference from the single-ring form, and its two ray derivatives,
obey (6). The actual evolved two-column forms are conserved in their
initial phase chart. The observed angle chart and its moving action
use the actual rows (14) of the companion proof, not the isolated
mode's frequency in place of those rows.

The tag lies in a truly quiet cavity of each stationary array copy.
Its angle, G,S, polar centroid and symmetric shape rows are genuine
material integrals. S=G_t remains exact, while the carrier-two
constant-current comparison and the linked polar row have errors
bounded by (5) of the companion proof, (12) there and (6) here.
Divide only after retaining the actual nonzero quadrupole, phase
coefficient, Wronskian and moment inverse. One can choose all three
errors below any predeclared positive physical action/current/jet
margin, with their conditioning constants included.

Let V=P^3. The actual coefficients per physical volume are

    j_raw=J_ring/V, Delta_raw=Delta_tag/V.                  (7)

For the whole-law angle reconstruction `Phi=3 E[n theta]`, the
corresponding phase and current coefficients are j=j_raw/3 and
Delta_density=Delta_raw/3, exactly as 0202. The reconstruction factor
does not multiply the actual averaged action or mechanical current.

All corresponding errors are divided by that SAME V. Thus their
relative errors do not acquire a spurious V factor. At every selected
finite P the coefficients are strictly positive. With one fixed ring
per cell they nevertheless tend to zero as P tends to infinity.
Even one densely prepared ring has action O(R), so P proportional
to R gives density O(R^-2). A per-ring coefficient is not an exact
nonzero density in that limit. The next file executes the registered
packing alternative and states separately which uniform estimates
it earns.

## 5. Complete ambient and shape: high-n does not erase the hybrid row

For the tag let G be its actual centered displacement angular moment,
S its spin and D its symmetric displacement moment. Define the hybrid
current by collapsing this material mass to its actual centroid while
keeping the ENTIRE continuous ambient complement. Taylor expansion of
the defining Fourier integrals, before variation, gives

    J_E=J_H+i K cross S/(2V)-i I_t K/(2V)+O(|K|^2),
    rho(X0-U_H0)=i K cross G0/(2V)-i D0 K/(2V)+O(|K|^2).   (8)

The full point-filter initial first moments vanish by high-n selection,
including the compensated common-K preparation. Actual Euler obeys
`J_E,t=-i K_j T_ij`, with its full symmetric pressure/stress T.
After the complete O(3) axial preparation, the zero-K linear maps
from the axial input to the symmetric stress and shape tensors vanish
by rotational covariance. The point first-gradient current therefore
remains zero on the finite window. Equation (8), not an assumed
canonical-spin identity, consequently gives

    J_H^(1)=-i K cross S/(2V),
    U_H^(1)=-i K cross [G0+integral S dt]/(2 rho V)
            =-i K cross G(t)/(2 rho V).                     (9)

The last equality uses the actual quiet-tag S=G_t. Symmetric shape
rows at higher order, all individual microscopic shape responses,
the polar companion and their energies remain in the complete action.
The mirror pair of 0202 cancels the tag's coherent polar companion,
not its variance or its energy.

Thus the complete high-harmonic packet has vanishing full low moments
while its actual centroid-plus-ambient observation has the nonzero
tag-current row (9). The ambient is precisely what makes the two
statements compatible. This proves an actual physical coupling
observation; it does not by itself supply an autonomous acoustic
constitutive equation or complete the parent continuum.
