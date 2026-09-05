# Actual elliptic material observation and its pressure-order repair

New analytic continuation within0155; no earlier attempt is changed.
The full positive target in README remains active. The exact conclusions
below concern the principal core and its first pressure order. A completed
high-order finite-action transfer is not inferred from these identities.

## 1. The Euclidean observation is an exact momentum shear

Use the conventions of `elliptic-core-algebra.md`, fix d and odd m>=3,
and put r=((1-sqrt(d))/(1+sqrt(d)))^m. Write

    D(t)=1+r²+2r cos(2m Omega t),
    a(t)=2r sin(2m Omega t)/D(t), omega=2 Omega.

The ACTUAL Euclidean moment angle factor previously derived is exactly

    F(t)=exp(i omega t)[1-i a(t)].                         (1)

Consequently Re(conj(F) exp(i omega t))=1. The pressure-spin observation
at its first nonzero order is in quadrature with exp(i omega t).
Their observation determinant is constant even though the Euclidean
angle has counter-rotating harmonics. It was therefore unnecessarily
restrictive to demand that r be smaller than every high-carrier error.

If radial normalization makes that physical spin P, and the underlying
principal oscillator coordinates satisfy qdot=P/M, Pdot=-M omega² q,
then the measured angle is theta=q-a P/(M omega). This is an explicitly
time-dependent *physical observation map*, not a Floquet log selection.
The actual pulled-back Hamiltonian is

    H(theta,P,t)=1/2 [M omega² theta² + 2 omega a theta P
                  +(1+a²-adot/omega) P²/M],
    det Hess H=omega²-omega adot.                          (2)

The adot term is essential. Eliminating P gives, with
B=1+a²-adot/omega,

    L=M/(2B)[(thetadot-omega a theta)²-omega² B theta²].     (3)

This retains the physical cross connection and its time derivative.
For instance |adot|<=4m Omega r/(1-r)², so the transparent sufficient
condition 2m r/(1-r)²<1 makes H positive. Finitely many further time
derivatives are bounded by C_j Omega^j m^j r/(1-r)^(j+1).
For any fixed d>0, an odd m satisfying all finitely many strict bounds
can be chosen FIRST. All m-dependent moment constants then remain
fixed as the carrier p increases; there is no unproved uniform-m
conditioning assumption.

For a carrier-dependent omega, the shear and its coefficients must be
differentiated as written. In particular det H=omega²-omega adot is
not simply omega². A statement about the scalar equation following
(3) also retains Bdot and the derivative of the cross term. No optical
gap identification is obtained by dropping those terms.

## 2. An independently defined calibrated material observable

The Euclidean mth moment is not the only actual measurable material
shape. A new registered candidate is the fixed-tensor moment

    Q_E=integral_tag (X+iY)^m rho dV,
    theta_E=(2/trG)(1/m) arg Q_E,
    (X,Y)=S(a,b).                                         (4)

S is determined by the stationary core strain. It is a physical tensor
that rotates with the entire field/marker in the ensemble, not a
canonical coordinate assigned a physical name. A pure m-lobed tag has
Q_E(t)=Q0 exp(-im Omega t) in the exact linear core. Its calibration
has unit response to a rigid physical rotation at the reference tag:

    delta Z=i[(trG)/2 Z +(G22-G11)/2 conjugate(Z)] delta phi,
    delta theta_E=delta phi.                              (5)

The second moment in (5) integrates to zero by the pure-m angular
selection for m>=3. Thus (4) measures a registered material-shape angle
with an explicit physical rotation response. It is not automatically
the absolute director required by any later continuum consumer.

The isotropic part of this tag sees only angular sector l=m-1 of the
displacement. Therefore its leading linear angle row has the SINGLE
phase exp(i[omega_0+m Omega]t)=exp(i2 Omega t), rather than (1).
The metric/shear factors are fixed real nonzero constants. They enter
the angle normalization and KKS, not the phase frequency.

## 3. First pressure order and what has actually been computed

Let T=trG, DeltaG=G11-G22, e_+=(1,i)/sqrt(2),
F=R^(m-1) exp(-R²/2)L_n^(m-1)(R²), and D_+=partial_X+i partial_Y.
Keeping the full longitudinal pressure equation gives the resonant
operator (cD/2)(Delta_R-R²), cD=Omega sqrt(T/p). The nonresonant
first-order velocity components are

    delta_e DeltaG/(4T) e_+ D_+² F,
    delta_e DeltaG/(4T) e_- Delta_R F,
    delta_e/8             e_- D_+² F.                     (6)

Their principal frequency denominators are respectively 2Omega,
2Omega and4Omega. The first and third have angular sector l+2 and
do not contribute directly to the isotropic part of (4). The second
does; it is not discarded. The Kelvin inverse has an additional
pressure projection correction

    eta=-GJ V/h + J grad div V/(h p² ell²)+O_d(delta_e²),   (7)

where eta=S xi_perp and ell^4=T/p³. In particular its l-sector
contribution is -i delta_e e_+ Delta_R F/(2hT). These real first-order
amplitude changes enter the angle normalization. Their material-time
phase contribution starts at second order on fixed Omega T_time.

For radial chi proportional to R² on the mode support, the exact
Laguerre Laplace calculation gives

    mean(R²)/2=2n+m+m/(2n+m),
    gamma_E=2Omega+[m/(2n+m)]cD+O_d,m(delta_e²),
    p² partial_p²(gamma_E²)
       =3 Omega² [m/(2n+m)]delta_e+O_d,m(delta_e²)>0.        (8)

This is the principal natural-scale positive carrier curvature of
the *calibrated physical* observation, not the unobserved eigenvalue.
The tag can be smoothly cut off outside a fixed large scaled radius;
its Gaussian-tail error must be made small relative to the displayed
coefficient before the carrier is selected.

The actual pressure torque is computed in physical coordinates:

    r_phys cross grad_phys P
       =(T/2) partial_theta P
        +(G22-G11)/2 (X partial_Y+Y partial_X)P.            (9)

The pure-m tag selects the first term of the leading angular-m
pressure; the second term has angular m-2 and m+2. The nonzero
oscillating material spin therefore starts at delta_e and has a
radial functional proportional to

    exp(-x/2) x^(m/2) P_n(x),
    P_n=L_n^(m-1)-2 partial_x L_n^(m-1).                  (10)

The full material spin includes position motion; the leading
pressure Hessian dI gives no torque, which explains the missing
order-one oscillating spin. Equation(10), rather than a kinematic
inertia, is the candidate normalization functional. The reference
moment Q0 may be made small by signed radial lobe weights while the
complete tag density remains nonnegative. This is allowed because
the unmarked odd moment vanishes. Its fixed-shape moment solve and
the higher pressure orders still require an explicit rank receipt.

## 4. Finite-time alternative, before any claim of a trapped spectrum

There is a concrete reason to prefer the registered finite-time
candidate over asserting an isolated normal mode. The opposite
polarization e_- with sector l-2 shares the principal frequency, but
its pressure correction has the opposite Laplacian sign. A second-
order coupling may encounter its inverted-oscillator continuum.
Absence of a first-order direct angular coupling is not an all-order
Fredholm theorem.

The available repair is the finite-time Euler propagator expansion:
solve each inhomogeneous order with the exact principal transport
propagator, retaining its polynomial time terms, rather than invert
an unproved localized spectral operator. At every fixed finite order
the sources are finite polynomial-Gaussian angular families; the
opposite polarization and any secular terms remain explicit. The
original third candidate in README licenses this representation.
Its next exposing calculation is the actual high-order observation/
spin row rank including these terms, and the full pressure remainder
in the fixed periodic cell. Equations(1)--(10) alone do not supply
that construction or complete0155. No production numerics or source
comparison has been used to choose between these candidates.
