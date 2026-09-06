# Direct Euler quadratic solvability coefficient

Author: main-model supervisor, 2026-09-06. Author calculation, not included
in the completed bounded wording correction of review0041. This continues
0037 route B using the exact Euler equations and the column adjoint pair.
It supplies an actual quadratic projection, not the still-open uniform
Euler modulation/remainder theorem.

## 1. Normalize the right and left threshold vectors

Let c=c0, f=f0 with the 0027 normalization int f^2 dr/r=1, and write

    b=2L/r^4,  d=L'/r,  Phi=r^2 b d,
    e=K0^-1 f=b d f/c^2,  g=-d f/c,
    N=int f e r dr=(1/c^2) int Phi f^2 dr/r>0.       (B1)

The radial pairing is <u,v>_r=int conjugate(u) v r dr. The exact column
equations have matrix H(eta,chi)=(b chi,d K0 eta), with generator
H partial_z. Its right-moving threshold vector and normalized dual are

    R_-=(e,g),
    ell_-(eta,chi)=[<f,eta>_r-<b f,chi>_r/c]/(2N).  (B2)

Indeed H R_-=-c R_-, ell_-(R_-)=1, and ell_- H=-c ell_-.
The last identity uses self-adjointness of K0 with the actual matched
exterior: K0(b d f)=c^2 f. The opposite vector R_+=(e,-g) has
ell_-(R_+)=0. Define Q=I-R_- ell_- for this solvability calculation;
then ell_- Q=0 and ell_-(H+c)Q=0. A correction is amplitude-normalized
by ell_-(W2)=0. This Q retains the opposite branch and the full other
radial complement; it does not discard them from the Euler equation.

The common physical density/azimuthal integration factor cancels between
the pairing and its normalization. At zero axial frequency, f is the
threshold streamfunction with constant exterior, not an L2(r dr) radial
velocity eigenfunction. Its derivative energy and the pairings in (B1)-(B2)
are finite. Here the formula acts on smooth threshold profiles and test
complements for which those pairings exist. Its measurable nonzero-k
version, and the remaining localization issue as k tends to zero, are
separated in supervisor-adjoint-and-transfer.md. No uniform k-derivative
bound follows merely by evaluating (B2) at the threshold.

## 2. Project the actual quadratic Euler vector field

For {u,v}=(u_r v_z-u_z v_r)/r, the nonlinear terms in (eta,chi) order are

    F2_eta=2 chi chi_z/r^4-{psi,eta},
    F2_chi=-{psi,chi}.                              (B3)

Insert the leading field (psi,eta,chi)=(f,e,g) A(z). The coefficient
of A A_z in (B3) is exactly

    n_eta=2g^2/r^4-(f' e-f e')/r,
    n_chi=-(f' g-f g')/r.                           (B4)

Applying (B2), with its radial measure retained, gives

    C=ell_-(n_eta,n_chi)
     =(1/(2N)) int [2f g^2/r^3-f f' e+f^2 e'
                          +(b f/c)(f' g-f g')] dr
     =(1/(2N c^2)) int f^3
                          [b' d+2b d'+2d^2/r^3] dr. (B5)

All cancellations here are pointwise; no boundary term has been dropped.
For example, f' g-f g'=d' f^2/c. With b,d as in (B1), the square bracket is

    4[(L'^2+L L'')/r^5-3L L'/r^6].                 (B6)

The exact steady-source coefficient from 0027 is

    J=-(4/c^3)[(L'^2+L L'')/r^4-3L L'/r^5].        (B7)

Thus the bracket in (B5) equals -c^3 J/r, and consequently

    C=-(c/(2N)) int J f^3 dr/r
     =-c beta/N=-2 sigma beta,                      (B8)
    beta=(1/2) int J f^3 dr/r,
    lambda'(c)=2N/c,  sigma=1/lambda'(c)=c/(2N).

This is a direct adjoint projection of Euler, subsequently identified with
the independently obtained steady coefficient. It does not infer the
projection from the requirement that a stationary profile fit KdV.
In the uniform-vorticity interior the bracket vanishes; outside the compact
vorticity support d=0 and the smooth cutoff gives no distributional edge
term. The transition region supplies the nonzero coefficient.

## 3. Scaling, complement and remaining estimate

In lab coordinates, the projected terms calculated above are

    A_t=-c A_z-2 sigma beta A A_z

before dispersion and uncomputed complementary/higher-order terms are
included. For the 0037 long-wave ansatz

    (eta,chi)=mu R_- A(X,T)+mu^2 W2(X,T)+...,
    X=(z-c t)/L_mu,  T=mu t/L_mu,

the quadratic projected term is -2 sigma beta (mu^2/L_mu) A A_X.
The leading contribution of the slaved correction is
(mu^2/L_mu)(H+c) partial_X W2, whose ell_- projection vanishes
algebraically. This explains the complement normalization at this order;
it does not construct W2 or bound it. Frequency-dependent eigenvectors,
the Bessel exterior correction, and possible nonlocal remainder terms
still require the actual time-dependent estimates.

Combining this earned quadratic coefficient with the separately derived
linear dispersion yields the leading solvability equation

    A_T+sigma partial_X(A_XX+beta A^2)=0.            (B9)

The direct coefficient calculation removes the need to assume the nonlinear
KdV coefficient merely from branch compatibility. It does not establish
convergence to (B9), a full approximate Euler solution, uniform control on
times L_mu/mu, scalar supersmoothness, or all-time Euler stability. Those
remain distinct constructions in the unchanged route and parent objective.
In y=X-sigma T, linearization of (B9) about A_* gives
a_T=+sigma partial_y L_* a, consistently with the 0041 sign correction.

The companion verifier checks the symbolic integrand, its exact J relation,
the right/left normalization, and sign-sensitive mutations. These checks
are local algebra evidence, not a substitute for the continuum estimate.
