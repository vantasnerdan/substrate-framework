# An actual isotropic stress memory and its next corrector route

## Exact preparation, not a fitted mean oscillator

Use the actual first-cell Euler/Lin and observed stress equations0170.
The background is the smooth, mean-zero, stationary two-wave field

    u=(A cos z+B sin y, A sin z, B cos y),
    curl u=-u, p=-|u|²/2.

Whole-field translations, rotations and time reversal form the same type
of stationary isotropic law used by the mean construction. The field is
an exact finite Fourier polynomial; every product harmonic in the following
time derivatives is retained. No Galerkin truncation, time step, sampled
orientation or empirical coefficient enters.

For the Kelvin displacement column D, write the first cell
chi=kappa_m D_l chi^{ml}. Its initial rate and constant forcing are

    s^{ml}=P(u_m e_l+u_l e_m), chi_t^{ml}(0)=-s^{ml},
    F^{ml}=P(p_m e_l+p_l e_m).

The first-cell rows are symmetric in m,l. Set
C chi=-(u.grad)chi+(Du)chi and Lw=-P[(u.grad)w+(Du)w].
The actual Euler/Lin recursion, with no inverse operator, is

    w0=-s, w1=Lw0+F, w_n=Lw_(n-1) for n>=2,
    chi0=0, chi1=w0, chi_(n+1)=C chi_n+w_n.

Here a subscript n denotes the nth time derivative at zero, not a Taylor
coefficient divided by n!. The first run independently checks chi2 and
chi3 against the complete second-order Jacobi equation, including its
pressure Hessian. The continuation checks every higher cell's divergence.

## Exact common-direction isotropic contraction

Let T_iljm multiply kappa_j*kappa_m*D_l in the unprojected physical
stress response. Averaging transverse polarization and the one common
direction gives

    R_iso=(2/15) sum_ij T_iijj
                -(1/30) sum_il(T_ilil+T_illi).

This follows from the exact second and fourth sphere moments, not from
averaging independently chosen microscopic wavevectors. Symmetry of the
cell rows gives a separate contraction:

    R_iso=<|u|²>/3 + Z/5-Y/15,
    Z=sum_ml <u_m chi_t,l^{ml}+p_m chi_l^{ml}>,
    Y=sum_il <u_i chi_t,i^{ll}+p_i chi_i^{ll}>.

Only the zeroth time derivative includes the energy term. The trace cell
is exactly sum_l chi^{ll}=-2t*u: its forcing is2P grad p=0,
P(u.grad)u=0, and [(u.grad)²+Hess p]u=0. The continuation verifies
this identity on every computed derivative. Both independent contractions
agree on the overlapping first three time derivatives.

## The first terminal-looking zeros do not give temporal closure

The symbolic-amplitude first run gives

    R_D(0)=-2(1+B²)/15, R_D'(0)=R_D''(0)=0, A=1.

Continuing the actual operator on A=B=1 gives

    R_D^(n)(0), n=0,...,4 = [-4/15,0,0,0,-2/25].

Thus the physical displacement column has

    X_tt=k²[-4/15-t⁴/300+O(t⁶)]D+O(k³)

after the whole-field law. The odd derivatives vanish under time reversal.
The nonzero fourth derivative is an actual response-memory coefficient,
not a stiffness guessed from a high-order frequency fit.

There is also an exact extension to every pair of real amplitudes.
The fourth derivative is a homogeneous polynomial of degree6: each
Euler/Lin operation is linear in u, and pressure is quadratic. Independent
half-period shifts make the whole-field averaged response even in each
amplitude. Proper rotation exchanging the two circular waves makes it
symmetric in A,B. The one-wave result0151 makes it zero on both axes.
The only possible polynomial is c*A²*B²*(A²+B²); the full exact A=B=1
calculation fixes c=-1/25. Therefore

    R_D''''(0)=-A²*B²*(A²+B²)/25.

This is nonzero whenever both waves are present. Background stationarity,
constant curl and full isotropy do not remove this particular physical
mean's temporal memory at second spatial order. Earlier checks of only
the initial coefficient or its first two time derivatives would miss it.

## Route verdict and executed continuation

The universal cancellation candidate is refuted by the displayed actual
Euler family. The stronger useful result is the derived physical stress
jet and its explicit mechanism. This does not refute C-CST-013, which
retains precisely such time connections, or the accepted conditional
Cauchy--Born continuum, or every possible Euler preparation/ensemble.

Method repair was executed: after the first three derivatives vanished,
the full Euler/Lin operator was continued to the fourth derivative rather
than assigning the zeros an all-time meaning. Representation change was
executed: independent fourth-rank and symmetric-cell contractions agree,
with an exact trace identity and independent Jacobi comparison. The next
materially different route is an actual stationary strain/preparation
corrector instead of the bare Kelvin-D/common-V cell data. It must solve
the full cell equations and preserve the measured mean and phase pairing;
declaring the residual zero is not a construction. Optical stationary
observation and band-edge repairs proceed in0176, and the joint physical
action retains this current in0172.

Route verdict: refuted for the claimed general stationary-law cancellation.
Evidence scope: exact full-Fourier initial-time physical response and its
two-wave counterexample, not a campaign no-go or exhaustion certificate.
