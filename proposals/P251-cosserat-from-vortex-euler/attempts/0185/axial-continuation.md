# Actual axial continuation and finite-action/localization licenses

The exact helical spectrum in helical-mode.md is at k0=-m/c. This file
derives, rather than assumes, what continues away from that carrier.
It does not assign a sign to the resulting carrier curvature.

## 1. The full nearby-carrier Euler operator

At fixed integer azimuthal m and real axial k, let epsilon_k=m+ck. In
cylindrical components the actual linear Euler generator is

    L(m,k)v=-P_(m,k)[i f epsilon_k v+2f Jv+f_r v_r h].     (1)

The background has no radial advection. Consequently(1) is a BOUNDED
operator on the radial solenoidal physical L2 velocity space; for example
its unprojected norm is bounded by (|epsilon_k|+4)C/c^2. The full pressure
projection has norm1. The transverse Fourier representation of the axial
Leray symbol has |(q_perp,k)|>=|k|, so its first and higher derivatives
are bounded for |k-k0|<|k0|/2. No radial wall or discretized inverse is used.

The solenoidal spaces vary with k. They are identified analytically using
the nearby orthogonal projections: P_k P_0 followed by the inverse square
root of P_0 P_k P_0 on Ran P_0. The inverse square root is its convergent
binomial series when ||P_k-P_0||<1. This retains the actual pressure and
gives a bounded analytic generator on one fixed Hilbert space.

At k0, the helical momentum block has tau_t=0, while the perpendicular
optical block is unitarily equivalent to the compact skew operator from
the positive radial problem. Its ground optical eigenvalue i sigma_0 is
simple and isolated from0 and all other radial eigenvalues. The coupling
from the stationary tau block is bounded and triangular, so this nonzero
eigenvalue remains isolated in the FULL velocity generator, not only its
tau=0 restriction.

Choose a small fixed spectral circle around i sigma_0. On that circle the
resolvent Neumann series for L(k)-L(k0) converges for sufficiently small
|k-k0|. Its contour integral is an analytic rank-one projection. This
constructs an actual nearby eigenvector and eigenvalue of the full Euler
operator. Smooth radial regularity follows from the full pressure equation;
all finite Sobolev norms needed below can be retained in this construction.

## 2. Actual Kelvin reconstruction, including the new tau row

Write the continued frequency as i sigma(k) initially without assuming
sigma real. Put nu(r)=sigma(k)+f(r)epsilon_k. In a neighborhood small enough
that |nu| stays bounded away from0, the exact Lin displacement is

    xi=v/(i nu)+f_r v_r h/(i nu)^2.                       (2)

Its divergence vanishes; the radial derivative of nu in the first term
cancels the helical derivative i epsilon_k of the second. The helical
momentum is no longer silently set to0. Taking the h component of the
pressure equation gives

    tau=-epsilon_k pi/nu.                                (3)

The radial component of v cross omega0 differs from 2f v_theta by f_r tau.
Using(3) and nu_r=epsilon_k f_r therefore yields exactly

    xi cross omega0=v+grad[pi/(i nu)].                    (4)

Hence v=P(xi cross omega0) on the actual fixed-Kelvin orbit, including the
off-helical axial momentum and full pressure return.

The continued complex KKS norm is nonzero by continuity from the positive
value in helical-mode.md. Euler conserves this exact pairing. If its
eigenvalue acquired nonzero real part, the mode/conjugate pairing would
grow or decay by exp[2 Re(eigenvalue)t], contradicting conservation.
Thus sigma(k) remains real. Its signed action stays positive in the
same physical lab clock. This is a local analytic spectral branch, not
proof of its curvature, a chosen Floquet winding, or a global band.

## 3. A genuine finite-action packet on R3

Choose a smooth spectral envelope a_L(k-k0) supported strictly inside the
constructed carrier interval, for example the Gaussian L exp[-L^2(k-k0)^2/2]
times a fixed remote cutoff. Form the Fourier superposition of the ACTUAL
continued displacement and velocity modes with their actual frequencies.
Translation invariance in z makes it an exact whole-space linear Euler/Lin
solution. Plancherel gives finite H^s norms and the complete action

    beta_packet=integral |a_L|^2 beta(k)dk>0,
    H_packet=integral |a_L|^2 h(k)dk>0.                   (5)

There is no assigned axial fiber length in(5). The background still has
its declared nonfinite total energy; the entire perturbation action is
finite on R3. The exact two prepared columns may dephase in time because
sigma varies with k; they are not called a single exact neutral eigenmode.

Use finite material label sheets with a nonnegative axial window of width
proportional to L and the stationary helical marker F_m. Its unperturbed
moment is stationary for EVERY f(r), since u.grad F_m=0. Near k0 the actual
perturbation frequency and its radial material transport differ smoothly
from the helical values. Taylor's formula under the spectral integral
therefore gives a controlled O_T(L^-1) difference of the physical angle,
spin and their required first time derivatives after their own packet
normalizations. The whole packet KKS scales as L. The literal tagged spin
and reference moment also scale as L; the angle itself does not.

The Gaussian axial overlap factors are part of those normalizations, not
set to1. The radial two-control system in helical-mode.md can use these
actual limiting factors instead of the periodic ones: its two independent
rows remain independent. At finite sufficiently large L the exact initial
moment/action equation is a small perturbation of that invertible linear
system and has a root by the finite-dimensional implicit function theorem.
All-time spin matching is then controlled to O_T(L^-1), rather than called
an exact packet identity. Nonzero reference quadrupoles and positive phase
action persist on the fixed interval. No carrier-second-derivative or
long-time bound is inferred from this first finite-packet construction.

## 4. Smooth radial localization is an actual Euler change

Choose f_R=chi(r/R)C/(c^2+r^2), with chi smooth, one on[0,1], zero above2.
The background u_R=f_R h is still EXACT stationary incompressible Euler,
with p_R,r=r f_R^2. It is not declared constant-helical-momentum in the
cutoff annulus. There

    delta tau_t=-[(c^2+r^2)f_R]'v_r

is retained in the full operator. Its outer vorticity carries the return
required by vanishing velocity outside2R. No confining wall or external
radial force was added. The velocity and every fixed derivative have
R-uniform global bounds.

For the helical eigenmode, the residual against this new full Euler
generator is supported in r>R before its FULL Leray return; its H^s norm
is bounded by C exp(-aR), a<m/c. The actual Euler energy estimate and
Duhamel therefore give a full-field finite-time bound C_T exp(-aR).
Exact Kelvin initial data for u_R differ from the original by the same
tail bound. For the actual finite packet the compact carrier interval
gives uniform small radial tails; choose L and then R after the fixed
time/physical margins. The identical Duhamel argument retains all pressure
and induced exterior velocities.

This supplies actual smooth radially localized Euler backgrounds with
controlled finite-time physical optical packet histories. It does not
claim an exact eigenmode of the cutoff field, compact support in the axial
direction of the BACKGROUND, or a Euclidean closed EPS tube. Those geometry
and continuum achievements are distinct from the mode/action result.
