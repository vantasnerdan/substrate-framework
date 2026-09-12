# Physical high-harmonic radiation barrier

## 1. Circulation and time normalization

For an axisymmetric no-swirl ring, write

    omega=omega_theta e_theta,
    zeta=omega_theta/r,
    dnu=r dr dz.

The source circulation used by Cao is

    kappa=integral zeta dnu=integral omega_theta dr dz.       (1)

The last integral is the flux of vorticity through a meridional cross-section,
so Stokes' theorem identifies it with the column circulation `Gamma`.  Thus
there is no extra `2*pi` between the two conventions:

    Gamma=kappa.                                             (2)

Reviewed P253/0052 fixes

    F(a)=Gamma/(2 pi),
    Omega_a=F(a)/a^2=Gamma/(2 pi a^2),
    omega_physical=Omega_a sigma.                            (3)

For mean radius `R`, core scale `a`, toroidal harmonic `n`, and normalized
column wave number `k`,

    delta=a/R,
    k=delta n,
    a=R k/n.                                                 (4)

Combining (2)--(4) gives the exact clock and physical frequency

    Omega_a=kappa n^2/(2 pi R^2 k^2),
    |omega|=kappa |sigma(k)| n^2/(2 pi R^2 k^2).             (5)

The density cancels from the frequency only after the physical Hessian/KKS
normalization.  Equation (5) does not provide that normalization; it converts
an already normalized dimensionless eigenfrequency.

## 2. Sharp transverse maximum on the Doppler shell

Use the Fourier/time convention from P253/0085.  A comoving mode of signed
frequency `omega` radiates on

    D(K)=|K|^2-[omega+c_g K_z]^2/c_EM^2=0.                  (6)

Put `s=sgn(omega)` and `t=K_z/|K|`.  The positive radial root is

    |K|=|omega|/[c_EM-s c_g t],                             (7)

whose denominator is positive for `|c_g|<c_EM`.  Hence

    K_perp(t)
      =|omega|sqrt(1-t^2)/[c_EM-s c_g t].                   (8)

Differentiating its logarithm gives

    d/dt log K_perp
      =-t/(1-t^2)+s c_g/(c_EM-s c_g t),

so the only interior critical point is

    t_*=s c_g/c_EM.                                        (9)

The endpoints give zero.  Therefore (9) is the global maximum and

    K_perp,max=|omega|/sqrt(c_EM^2-c_g^2).                 (10)

The weaker bound `|K|<=|omega|/(c_EM-|c_g|)` remains true but is not the
sharp Bessel argument.  This distinction materially enlarges the possible
finite radiation window.

## 3. Exact vector Bessel orders

Let `K_perp=(rho cos phi,rho sin phi)` and use
`exp(-i K dot x)`.  For

    j=e^(i n theta)(j_r e_r+j_theta e_theta+j_z e_z),

define fixed-Cartesian helicities

    j_+=j_x+i j_y=(j_r+i j_theta)e^(i(n+1)theta),
    j_-=j_x-i j_y=(j_r-i j_theta)e^(i(n-1)theta).          (11)

Jacobi--Anger gives the exact angular integral

    integral_0^(2 pi) e^(i m theta)
       e^(-i rho r cos(theta-phi))dtheta
      =2 pi(-i)^m e^(i m phi)J_m(rho r).                  (12)

Thus `j_+`, `j_z`, and `j_-` have orders `n+1,n,n-1`, with their exact
phases.  This conclusion survives the transverse projection covariantly.
Indeed `K_- jhat_+`, `K_+ jhat_-`, and `K_z jhat_z` all have scalar phase
`e^(i n phi)` in `K dot jhat`; multiplying that scalar by
`K_+,K_z,K_-` in `K(K dot jhat)/|K|^2` returns precisely phases
`n+1,n,n-1`.  No `n+/-2` or `n+/-3` order survives.  The real `-n` member is
the conjugate calculation.

This is an order ledger, not a current-size estimate.  The physical current

    delta J=g(chi_g delta u+u_g delta chi)                 (13)

still needs its radial/axial profile, KKS normalization, constraint rows, and
Leray/Hodge response.

## 4. A rigorous fixed-margin Debye bound

The generating-function contour representation permits a vertical shift by
`alpha>0` and gives, for `0<x<m`,

    |J_m(x)|<=exp[-m alpha+x sinh alpha].                  (14)

Optimizing with `cosh alpha=m/x` yields

    |J_m(x)|
      <=exp{-m[acosh(m/x)-sqrt(1-(x/m)^2)]}.              (15)

For a fixed `eta in (0,1)` and `x/m<=1-eta`, define

    c_eta=acosh[(1-eta)^(-1)]
             -sqrt[1-(1-eta)^2]>0.                        (16)

Then `|J_m(x)|<=exp(-c_eta m)`.  For the vector ledger the smallest order is
`m=n-1`.  If the current support obeys

    r<=R+C_core a,                                        (17)

the sufficient uniform barrier is

    K_perp,max(R+C_core a)<=(1-eta)(n-1).                 (18)

At ratio one the natural scale is Airy `n^(-1/3)`.  Above one it is
oscillatory with `n^(-1/2)` envelope away from caustics.  Neither regime earns
exponential suppression.

## 5. Fixed-J modes give only a finite ceiling

Substituting (5), (10), and `a=Rk/n` into (18) gives

    Q_fixed(n;k,sigma)
      ={kappa |sigma(k)|

         over 2 pi R k^2 sqrt(c_EM^2-c_g^2)}
        {n(n+C_core k) over n-1}
      <=1-eta.                                            (19)

For each member of a rational-ray pair, (19) defines an exact integer ceiling
`N_rad` after `n_i=N p_i`, the exact crossing frequency, and the carrier speed
are inserted.  Since its left side grows like `n`, while the Cao speed also
shrinks the square-root denominator, no fixed-constant `N->infinity` Bessel
barrier exists for fixed radial index and fixed massive `k`.

The useful positive theorem is finite:

    N_lower<N<min(N_EM,N_rad),                             (20)

with every graph, IVT, KKS, response, and gate threshold included in
`N_lower`.  Whether (20) contains an integer is a physical parameter question,
not a consequence of taking the ring thinner.

## 6. Fixed-ell/high-J has a positive structural regime

Consider the adjacent pair

    (J,n_1=J ell,k_1=k_*),
    (J+1,n_2=(J+1)ell,k_2=k_*(1+1/J)),                    (21)

on a common ring with `delta=k_*/(J ell)`.  Reviewed P253/0066/0073 gives, at
fixed positive massive `k`,

    sigma_J(k)=k L_Phi/(pi J)+O(k/J^2),                   (22)

where `L_Phi` is the dimensionless normalized-column coefficient.  For either
member of (21),

    sigma=delta ell L_Phi/pi+O(delta^2),

    |omega|
      =kappa ell L_Phi/(2 pi^2 R^2 delta)+O(1),

    n=k_*/delta+O(1).                                     (23)

At an eventual exact crossing the two modes have one common physical
frequency.  For an actual finite subluminal integer `J`, define

    q_0(J)=C_0/sqrt(c_EM^2-c_g(J)^2),

    C_0=kappa ell L_Phi/(2 pi^2 R k_*).                   (24)

The joint high-`J`/thin-ring transfer must quantify the simultaneous support
test, with `min(|n_1|-1,|n_2|-1)`, as

    Q_i(J)=[K_perp,max(R+C a)]/(|n_i|-1)
          =q_0(J)+R_i(J),
    |R_i(J)|<=C_asym/J,                 i=1,2.            (24a)

Because `c_g(J)` grows logarithmically and subluminality ends at finite
`J_EM`, (24) is a finite-`J` leading predictor, not a literal fixed-`c_EM`
limit as `J` tends to infinity.  Fix `eta_rad in (0,1/2)`.  The exact
buffered Debye criterion is

    q_0(J)<=1-2 eta_rad,
    J>=J_asym,
    J_asym>=C_asym/eta_rad.                               (25)

Then (24a) gives `Q_i(J)<=1-eta_rad` for both pair members.  The present
fixed-`k` graph input does not supply `C_asym` or `J_asym`; those constants are
part of the open joint transfer.

Equivalently, first require `C_0<(1-2 eta_rad)c_EM`, then

    c_g(J)^2
      <=c_EM^2-[C_0/(1-2 eta_rad)]^2.                     (26)

Together with (25), this leading speed bound leaves the exact margin
`1-eta_rad` only at an actual integer for which the quantified joint-transfer
remainder is small enough.

The Cao value law has

    c_g(J)=alpha[log J+C]+o(1),
    alpha=kappa/(4 pi R).                                 (27)

Thus (26) defines a finite ceiling `J_rad`, in addition to the subluminal
ceiling `J_EM`.  The actual charged target is

    J_lower=max(J_graph,J_Riesz,J_KKS,J_gate,J_asym,...)
      <J<min(J_EM,J_rad).                                 (28)

Equation (24) is a structural finite-`J` predictor improving the fixed-index
scaling because the high radial index reduces `sigma`.  Equation (27)
prevents a fixed-`c_EM` subluminal `J`-to-infinity limit.  One or finitely many
protected gate scales may exist only when (28) is nonempty after the joint
transfer supplies `C_asym` and `J_asym`.  The joint high-`J`/thin-ring Riesz
transfer and exact crossing are separate missing constructions; the column
asymptotic does not create them.

## 7. From a shell trace to gate loss

Suppose the physically normalized current profiles, their required
derivatives, and the transverse projection contribute at most `C H^M`, where
`H=N` in the fixed-`J` route and `H=J` in the fixed-`ell` route.  Equations
(12), (15), and compact radial/axial support then give

    ||P_T delta Jhat||_(shell trace)
      <=C H^M exp(-c H).                                  (29)

This is still not an `L2` Maxwell eigenmode: any nonzero real-shell trace
violates the exact dark-current condition from P253/0085.  To obtain a
finite-gate statement, an outgoing limiting-absorption/Fermi-golden-rule
identity or a time-domain flux identity must prove

    P_rad<=C_trace ||P_T delta Jhat||_(shell)^2            (30)

with the coarea weight and physical constants, and the analyzer construction
must prove

    T_gate<=C H^G.                                        (31)

Then the fractional action loss satisfies

    T_gate P_rad/A_mode
      <=C H^(2M+G+M_A)exp(-2cH),                          (32)

provided the physical mode action has an inverse polynomial bound
`A_mode^(-1)<=C H^M_A`.  Define `H_loss` as a quantitative threshold above
which the right-hand side of (32), with all proved constants retained, is
below the selected tolerance.  The finite-gate conclusion requires an
integer

    max(H_graph,H_Riesz,H_KKS,H_gate,H_asym,H_loss,...)
      < H < min(H_EM,H_rad).                              (32a)

Exponential asymptotics alone do not prove that this intersection is nonempty
because its upper endpoint is finite.  An uncontrolled gate time, missing KKS
normalization, absent flux conversion, or empty intersection in (32a) leaves
(32) open.

Exact cancellation of the shell trace is a distinct BIC/dark-current route.
It may yield a real eigenmode only after the coupled weighted resolvent and
fluid feedback are solved.

## 8. Route verdicts

**Route A (fixed-J massive pair): established as a finite-ceiling theorem.**
The exact physical clock and sharp transverse shell maximum give (19)--(20).
They refute an asymptotic Bessel barrier for this route while preserving a
possible finite protected window.  No nonempty pair or gate is established:
crossing, response, and the integer window remain open.

**Route B (fixed-ell/high-J pair): blocked at a named same-carrier
construction, with a positive structural inequality established.**  Equations
(21)--(28) derive the finite-`J` predictor and show that exponential
large-order suppression follows only if the joint transfer supplies the
uniform remainder and an actual integer clears `J_asym` below the joint
`J_EM,J_rad` ceiling.  The actual Cao pair still requires the joint
high-`J`/thin-ring Riesz crossing, physical KKS normalization, mixer response,
and quantitative lower thresholds.

**Route C (shell trace to resonance/gate loss): blocked at a named analytic
construction.**  Equations (29)--(32) are the exact implication once the
physical current polynomial bound, outgoing flux/LAP identity, and gate-time
bound are supplied.  A nonzero exponentially small trace supports only a
quasimode/resonance or finite-time leakage estimate, not an exact eigenmode.

The next failure-derived route after an empty Bessel window is the complete
source-specific shell trace followed by a one-block outgoing Maxwell
LAP/Feshbach calculation.  The full Issue #203 electron and neutrino objective
remains active.
