# Executed exact Euler memory, with the unresolved initial state retained

This construction supplies an autonomous Euler mechanism for a projected history. No target oscillator, constitutive coefficient, interpolation controller or future history is supplied. The concrete solution is periodic and finite-energy per cell; it is not an isolated localized carrier. Its use is to execute P1's nonlocal/retained-state route and expose what a subsequent particle reduction must preserve.

## 1. Exact full nonlinear split

On a periodic cell let L be the Leray projector and define B(a,b)=-L[(a dot grad)b]. For smooth Euler let Pi be an orthogonal divergence-free Fourier projection commuting with derivatives, Q=1-Pi, v=Pi u and w=Qu. Exactly,

    vdot=Pi B(v+w,v+w),
    wdot=Q B(v+w,v+w), u(0)=v0+w0.                       (1)

These are the actual full Euler equations in resolved/unresolved coordinates; their pressure is recovered from the same u. If the unresolved evolution for a supplied v history is solved, substitution gives history dependence AND a dependence on w0. A general orthogonal-dynamics theorem is not presumed from writing(1).

Their kinetic energies have the exact exchange law

    d[ rho ||v||_2^2/2 ]/dt
      =rho integral w dot [(v dot grad)v]
       +rho integral w_i w_j partial_j v_i
      =-d[ rho ||w||_2^2/2 ]/dt.                         (2)

To derive it, start from rho integral partial_j v_i u_j u_i and use div v=div w=0. The v-v self term and w_j v_i partial_j v_i are exact divergences. A memory model dropping w can therefore discard both energy and the actual acceleration, as0001 independently exposes.

## 2. A full Euler family where the elimination is solved

Use dimensionless coordinates on T^3 and any smooth real periodic U(y). For every smooth real periodic w0(x,y),

    u(t,x,y,z)=(U(y),0,w0(x-t U(y),y)), p=constant.        (3)

This solves ALL nonlinear Euler equations globally: the field is z-independent, the x-velocity is stationary and its only advective derivative is zero; the z-equation is w_t+U(y)w_x=0. Smoothness holds for all finite t, though spatial derivatives can grow. The amplitudes need not be infinitesimal. This invariant triangular class is narrower than arbitrary 3D perturbations.

For one x Fourier harmonic k!=0, put w=Re[e^(i k x)g(t,y)]. Then

    gdot=L_k g, L_k=-i k U(y), g(t)=exp[-i k U(y)t]g0.    (4)

On L2(S1,dy/(2pi)), L_k is a bounded skew-adjoint multiplication operator. This makes the following elimination rigorous without an unproved PDE inverse. Define a=<1,g>, P g=a*1, Q=1-P, mean U=0, b=Q L_k1 and A=Q L_k Q on QL2. Duhamel gives

    a'(t)=F(t)-integral_0^t K(t-s)a(s)ds,                (5)
    F(t)=-<b,exp(t A) Qg0>,
    K(t)=<b,exp(t A)b>.

Here <f,g>=average conjugate(f)g. Equation(5) follows from P L_k P=0 and P L_k Q=-(Q L_k P)^*. The unresolved initial term is essential: it is not a stochastic noise assumption or discarded preparation cost. All operators are bounded, so their exponentials and block Duhamel identities hold globally on the stated Hilbert space.

## 3. Closed-form kernel, continuum band and exact energy transfer

Take U(y)=sin y and k=1. For the natural initial profile g0=1, F=0 and

    a(t)=average exp[-i t sin y]=J0(t),                  (6)
    K(t)=J1(t)/t, K(0)=1/2.                             (7)

The integral representation in [DLMF10.9.1](https://dlmf.nist.gov/10.9.E1) supplies(6). The Laplace transform may also be obtained directly for s>0 by averaging 1/(s+i sin y):

    Ahat(s)=1/sqrt(s^2+1).

Taking the Laplace transform of(5) therefore gives

    Khat(s)=1/Ahat(s)-s=sqrt(s^2+1)-s.                   (8)

To verify(7) without a fitted kernel, write Fhat(s)=Laplace[J1(t)/t]. The derivative identities in [DLMF10.6](https://dlmf.nist.gov/10.6) give Laplace[J1]=1-s/sqrt(s^2+1). Hence Fhat'(s)=s/sqrt(s^2+1)-1; Fhat(infinity)=0 fixes exactly(8). Its t=0 limit is1/2 from the Bessel power series. This is an exact all-time Volterra equation, not a Markov or harmonic approximation.

For a physical shear amplitude U_s and x wave number k, beta=k U_s has frequency units. Then a(t)=J0(beta t), K_beta(t)=beta J1(beta t)/t and K_beta(0)=beta^2/2. A transverse cell length enters U(y)=U_s sin(y/L_y); it is not an empirically chosen particle scale.

The self-adjoint frequency operator i L_k has continuous spectral support [-|beta|,|beta|], while the skew-adjoint evolution generator L_k has spectrum i[-|beta|,|beta|]. The frequency operator has the arcsine observation measure dnu=dx/[pi sqrt(beta^2-x^2)]. It has no nonzero L2 eigenfunction concentrated at one frequency, because a nonconstant sine level set has measure zero. The square-root branch points in Ahat demonstrate that no finite-dimensional constant-coefficient LINEAR system can reproduce this exact transfer for all times: such systems have rational resolvents. This does not forbid every nonlinear representation or a controlled finite-band approximation.

Unitarity gives ||g||_2^2=||g0||_2^2=1 and orthogonality gives

    |a(t)|^2+||Qg(t)||_2^2=1.                            (9)

Resolved apparent attenuation is exact transfer into unresolved shear modes, not dissipative loss. In the full real Euler field the z-energy has a common factor rho*cell_volume/4; the x-shear energy is separately constant. All pressure work remains consistent with(1)–(3).

## 4. Natural unresolved perturbations and the failure of deterministic local closure

Take g0^plus=1+i epsilon sin y and g0^minus=1-i epsilon sin y with any small real epsilon. Both have a(0)=1 and equal total z-energy. Their velocities are real after taking Re[e^(ix)g] and are arbitrarily close in every fixed Sobolev norm as epsilon tends to zero. Their exact projected derivatives are

    a_plus'(0)=+epsilon/2, a_minus'(0)=-epsilon/2.         (10)

Indeed average L_1 g0^plus/minus has exactly these values. More fully, a_plus/minus(t)=J0(t)+/-epsilon J1(t), and F_plus/minus(t)=+/-epsilon J1(t)/t with the limits in(10). For beta=k U_s the forcing is +/-epsilon J1(beta t)/t; its Laplace transform is +/-epsilon Khat(s)/beta when beta!=0. These follow by differentiating the actual sine integral in(6) and applying(8). Their difference is Qg0, the term retained in F of(5). A finite resolved initial observation and its energy therefore do not determine the projected future in this actual Euler family. Setting F=0 is a precise initial-state restriction, not a universal closure.

Equation(6) is obtained from one simple initial field, not from selecting initial data to match an intended J0 history. Equation(10) describes a natural nearby family and quantifies the additional state needed for exact prediction.

## 5. Positive result and next construction

The full nonlinear split/energy exchange and the executed triangular Euler Volterra reduction are established as exact identities on their stated solution classes. They supply a real unsteady Euler action-compatible state evolution with a derived memory kernel and initial-state term. They do not supply a localized carrier, a stable generic 3D background, a quantum probability law, electric/weak currents or a relativistic dispersion.

For P1 on a physical R1/swirl carrier, the next construction is its actual projection with analogous bounded or controlled evolution of the complement, preserving the pressure/current and energy exchange. For P4, the continuous convective band here is a derived structured-background response; it is not a photon or neutrino band and cannot be assigned a finite invariant speed by renaming beta. The same-carrier orbital action in0005 is the variational input to that next projection.

## Independent review correction

Review0011 established the exact Euler, Volterra, kernel and unresolved-state forcing constructions. Its bounded terminology correction above distinguishes the imaginary evolution spectrum from the real frequency spectrum. No equation, module or test changed, and the existing four-test receipt remains applicable.
