# Common-circle multi-CK geometry and full-pressure packet transfer

The parent registered the failure-generated positive multi-CK expansion
and validated263/12 before the construction below. It replaces the
fixed Taylor--Green higher jets by a controlled Lundquist approximation;
the entire fields remain exactly constant-lambda Euler.

## 1. Positive angular moments at one exact common core circle

Fix lambda,U>0 and a finite integer M>=2. The angular identity is

    J0(lambda sqrt(x²+z²))
      =E[cos(sqrt(t)x)cos(sqrt(lambda²-t)z)],                  (1)

where t/lambda² has the arcsine distribution on[0,1]. Its moments
are lambda^(2h) binomial(2h,h)/4^h. Choose M+1 interior positive
Gauss--Chebyshev nodes t_j^0 and weights1/(M+1). They integrate
polynomials to at least degree M, with positive weights.

For large R select a positive zero j_(0,n_j) of J0 nearest to
sqrt(t_j^0)R and set kappa_j=j_(0,n_j)/R. Bessel-zero spacing is
bounded directly by Sturm comparison for f=sqrt(x)J0(x), which
satisfies f''+(1+1/(4x²))f=0. On x>=X consecutive zeros have
spacing between pi/sqrt(1+1/(4X²)) and pi. Hence
kappa_j=sqrt(t_j^0)+O(1/R). The selected nodes
remain distinct and strictly inside(0,lambda). Solve the square
Vandermonde moment system

    sum_j w_j kappa_j^(2h)
       =lambda^(2h)binomial(2h,h)/4^h, 0<=h<=M.               (2)

Its solution is close to the strictly positive original weights;
therefore w_j>0 for sufficiently large R and sum w_j=1. This is an
exact finite linear system, not numerical quadrature used as evidence.

For k_j=sqrt(lambda²-kappa_j²) define the entire axisymmetric field

    u=sum_j A_j(k_j J1(kappa_j r)sin(k_j z),
                 lambda J1(kappa_j r)cos(k_j z),
                 kappa_j J0(kappa_j r)cos(k_j z)),
    A_j=U w_j/[lambda J1(kappa_j R)].                         (3)

Every mode has div u=0,curl u=lambda u. Every term has the SAME
critical circle r=R,z=0, since J0(kappa_j R)=0. Their sum has core
speed U. Equations(2),h=0,1 give its EXACT circular transverse jet
Omega=lambda U/2. No singular Y-Bessel field is introduced.

As R varies in a small interval, keep the integer n_j fixed and
use kappa_j(R)=j_(0,n_j)/R; solve(2) smoothly for w_j(R).
Positivity and distinctness persist. The exact circular core return
number is lambda R/2, so select nonresonance, or Diophantine core
rotation, by an arbitrarily small change of R. No Bessel-zero
arithmetic is presumed. This is a direct many-mode alternative to
0145's two-mode circularity IFT.

## 2. Local accuracy and global bounds are separate, both controlled

Let y_j(x)=J1(kappa_j(R+x))/J1(kappa_j R). It obeys

    y_j''+y_j'/(R+x)+(kappa_j²-1/(R+x)²)y_j=0,
    y_j(0)=1, y_j'(0)=-1/R.                                  (4)

On every fixed interval |x|<=D, comparison of this ODE with cos(kappa_j x)
gives a uniform finite-order C^s error C_(D,s)/R. The constants can
be chosen uniformly for kappa_j in[0,lambda]; this uses the ODE,
not division by kappa_j. Positivity and sum w_j=1 retain this bound
under summation. The accompanying velocity components follow by
the exact curl/streamfunction relations, with the same bound.

The angular moments(2) make the even normal Taylor jets agree with
Lundquist through degree2M. On |x|,|z|<=D their full analytic error,
with any fixed number s of derivatives, is bounded by

    C_s U lambda^s exp(C lambda D)
      (C lambda D)^(2M+2-s)/(2M+2-s)! + C_(D,s)/R.             (5)

The factorial bound follows from the defining cosine series and
positive weights of total mass1. Thus M can be chosen AFTER the
fixed optical packet and any required error margin, then R chosen
after M. This handles anisotropic higher jets quantitatively;
they are not declared negligible merely because the core jet is
circular. On larger annuli |r-R|<=R/2, normalized Bessel asymptotics
give R-uniform velocity and finite derivative bounds. Globally,

    ||u||_(C^s) <= C_s U(1+lambda)^s sqrt(1+lambda R),          (6)

after all selected kappa_j R are sufficiently large. Indeed
|J1(j_(0,n))| is bounded below by c/sqrt(j_(0,n)). To see this
without an unquantified asymptotic, the same f has E=f²+f'²>0
and |E'|<=E/(4x²). Integration gives uniform positive upper/lower
bounds for E on x>=X. At a J0 zero, E=xJ1(x)². Each entire
J0,J1 and fixed derivative are bounded. The positive weight sum
controls the sum. Constants do not grow with later periodic
quadrature accuracy. A fixed finite M has no hidden normalization
chosen after that approximation.

## 3. An actual invariant solid torus

The streamfunction of(3) has an elliptic maximum at the common
circle. Its local normal flux function converges to J0(lambda r_n)
with the chosen high Taylor order, plus finite-R curvature terms.
The straight circular limit has

    poloidal angular speed O=U J1(lambda r_n)/r_n,
    longitudinal speed W=U J0(lambda r_n),
    flux action I=U r_n J1(lambda r_n)/lambda.

Its toroidal-section rotation is R O/W. At the center,

    (1/R)d(R O/W)/dI=lambda³/(8U)>0.                         (7)

This follows directly from the Bessel series through quartic order.
The actual finite-radius return map and its action derivative depend
continuously on these jets and on a small regular flux annulus, as
in0145's coarea/action calculation. Fix a sufficiently small regular
annulus, choose M and R so the strict margin(7) persists, then choose
a Diophantine boundary flux level. The exact circular core is selected
nonresonant using the continuous R freedom above.

The actual positive section measure is u_theta dr dz. The precise
EPS return-map KAM theorem and Moser measure identification already
read and applied in0145 now apply to this SAME actual field. They
supply a robust invariant unknotted solid-torus boundary and its
nondegenerate elliptic periodic core. No packet is assigned to a
distant knot. Arbitrary-knot coexistence, if desired, remains distinct.

## 4. High-angular-harmonic pressure localization removes the norm circularity

A plain global energy estimate would give exp(C sqrt(R)T)/R and
does not prove convergence. The angular harmonic structure supplies
the needed sharper license. In cylindrical components, the
axisymmetric background preserves each integer angular harmonic n.
The packet constructed below uses |n|/R in a fixed interval about p.

For a scalar pressure harmonic the full Laplacian is

    Delta_n=partial_r²+r^-1partial_r+partial_z²-n²/r².          (8)

There is no radial wall. Its weighted energy identity contains the
positive term integral(n²/r²)|pi|² r dr dz. Choose a radial weight
W_R that equals1 near r=R, decays with logarithmic slope at most
a|n|/r through R/2<r<3R/2, and is a constant <=exp(-a_0 pR)
outside a slightly larger core annulus. Here 0<a<1/4 and a_0>0
are fixed. The weighted pressure estimate follows by testing(8)
with W_R² conjugate(pi) and absorbing the derivative-of-weight term:

    (1-4a²)||W_R (n/r)pi||² + (1/2)||W_R grad_(r,z)pi||²
      <= C ||W_R F||²                                        (9)

when the pressure right-hand side is div F; the cylindrical volume
measure r dr dz is used. Vector components shift the effective angular
orders by at most1, harmless for the large |n| in question. The same
argument after commuting a fixed number of derivatives gives the
weighted finite-order Sobolev estimates. This is the full pressure
inverse, not a local replacement for Leray.

For linear Euler, write its pressure source in divergence form

    Delta pi=-2 partial_i[(partial_j u_i)v_j].                (10)

The weighted kinetic-energy identity has pressure term
integral pi v dot grad(W_R²). Equation(9) controls it without a
derivative loss because |grad W_R|<=a|n|W_R/r. The transport-weight
term costs only C_p times weighted energy where the weight varies:
there the field and its derivatives have the R-uniform annular
bounds of§2. In the far region the weight is constant; its possibly
large strain is bounded by(6) times exp(-2a_0 pR) times the global
unweighted energy. Consequently

    E_weighted(t)<=C_(p,T)[E_weighted(0)+integrated local forcing
                 +exp(-a_0 pR+C sqrt(R)T)E_global(0)].         (11)

The global energy bound itself follows from(6). Constants C_(p,T)
may grow with the already fixed p,T, but NOT with R. The pressure
source and forcing tails obey the same estimate. This is precisely
the mechanism absent from exp(C sqrt(R)T)/R. For fixed p,T, the last
term of(11) tends to zero, while the local comparison error is
C_(p,T)/R. Both carrier derivatives obey the same argument: they
differentiate preparation coefficients, not the fixed background.

An equivalent proof uses the Green kernel of(8): separated annuli
have angular barrier exp[-|n||log(r/r')|], with the same absorption
bound. No global physical perturbation is discarded or prevented
from returning by an artificial boundary.

## 5. From the finite straight packet to the actual torus

Fix the finite straight preparation, its full action/marker and all
strict physical margins FIRST. Choose a smooth transverse cutoff
and potential representation for its solenoidal initial displacement,
with all discarded Schwartz tails below the prescribed C² H^s
margin. The potential is obtained from the full-space inverse curl
on the nonzero longitudinal band. Apply the cutoff to the potential,
not directly to the displacement. This keeps incompressibility
exactly and records the tail error. Exact Kelvin preparation follows.

Periodize longitudinally around the circle of circumference2pi R.
Equivalently sample the smooth spectral weights at q_n=n/R with
the corresponding1/R normalization. Poisson summation gives the
original packet plus its remote Schwartz images. Their C² carrier
errors tend to zero faster than every inverse power of R/L.
The background angular harmonic n is then EXACT, not a noninteger
phase imposed on a closed circle. Choose p_*R integer if convenient;
alternatively place the smooth reference marker in a proper arc,
where its prescribed phase is globally smooth on its support.
The carrier parameter p remains continuous through the smooth
sampled weights a_L(n/R-p); it is a preparation parameter, not an
assertion of a continuously variable toroidal eigenvalue.

Use the exact cylindrical Piola/potential pushforward from the
straight tubular chart. Its metric/Jacobian and material frames
are kept, giving C² data errors O(1/R) at the already fixed p,L.
Full action(4) becomes the complete volume KKS of this finite
toroidal displacement; its geometric and image errors are controlled
by the same local coefficients and tail norms.

The pressure comparison itself is explicit. With x=r-R,q=n/R,

    Delta_n=partial_x²+partial_z²+partial_x/(R+x)
                  -q²R²/(R+x)²,
    Delta_n-(Delta_perp-q²)
          =partial_x/(R+x)+q²[1-R²/(R+x)²].                  (12)

On the fixed normal window this coefficient difference is O(1/R),
including its required derivatives. Apply the resolvent identity
to(12), using the full weighted inverse(9) and the straight Helmholtz
inverse. It gives C_(p,D,s)/R times the fixed forcing norm, plus
the already bounded exterior tails. The Piola Jacobian coefficients
are treated in the same divergence-form identity. Thus both the
local pressure change and remote pressure return are derived,
rather than inferred from local velocity closeness.

Compare actual Euler histories using(11). First fix a large normal
window containing the straight packet's uniformly small weighted
tails on the finite observation interval. Then choose M so the
factorial term of(5), multiplied by the finite C_(p,T) and observation
condition numbers, is below the desired physical error. Choose R
last so its local1/R, image, metric and far-axis terms are below
that error as well. Enlarging R does not worsen the local exponential
constant in(11), so this ordering is not circular.

All material sheets and the physical tag are pushed by the actual
flow. On the torus the measured collective angle is the registered
central quadrupole of each material sheet in its declared physical
core frame; the observed spin is the actual TOTAL tag angular
momentum about its actual centroid, projected on the central core
tangent. Across the finite arc the frame difference is O(L/R).
The EXACT frame/centroid terms are retained and bounded below the
chosen margin, rather than calling their limit identities exact.
Thus the packet lies in this actual robust invariant torus and the
finite-time physical angle/spin bounds transfer at their own scale.

More precisely, the localized displacement and marker can be supported
in a proper interior sub-tube, while the induced Euler velocity and
pressure retain their full exterior tails. Initial finite-amplitude
Kelvin preparation pushes the actual vortex-tube boundary by the
volume-preserving material map. Smooth Euler's transported vorticity
then carries that same material tube throughout the controlled time
window. The perturbation is not claimed stationary; its background
is the exact stationary torus field. The measured registered sheet
angle remains distinct from an absolute vorticity-core director.

## 6. Periodic stationary ensemble and remaining scope

Once the finite R,M field and packet are fixed,0145's Herglotz
quadrature applies to each of the finitely many entire CK terms.
Symmetric rational-direction approximation preserves exactly
curl u=lambda u, realness and a common global period. The sum has
a finite global C^s bound independent of approximation accuracy.
Use0145's full-space pressure/Duhamel comparison with its fixed-field
uniform bounds, and choose accuracy below the topology, packet,
angle-chart, spin and second-carrier margins. This supplies an actual
smooth periodic stationary Euler background with this same marked
unknotted torus and controlled finite-action Euler packet histories.

Uniform translation in its cell and Haar rotation of the WHOLE field,
packet and material labels gives a finite-variance stationary isotropic
law with positive density of these tori. No Gaussian or translation-
ergodicity statement is inferred. Distant prescribed knotted tubes
may coexist by0145, but are not the packet's optical tube.

The local spin/registered-angle/action theorem has now acquired a
finite-action actual global geometry and full-pressure transfer. It
does not supply an absolute director, a universal autonomous optical
band, or the parent's coupled acoustic/optical continuum. Those
remain separate physical observation/closure obligations.
