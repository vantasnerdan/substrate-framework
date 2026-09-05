# Positive physical D energy from an actual finite Euler preparation

This is the positive construction reached after the preceding route
repairs. It establishes the specified displacement-column energy/response
interface, not a complete acoustic or Cosserat action. In particular its
common-V and physical tag rows are not inferred from a D-column sign.

## 1. A degree-seven stationary Euler corrector with small range energy

Take B=1 and A=1/100 in the actual unequal two-wave field. Start at the
A=0 endpoint to choose a FINITE polynomial, then solve the actual moment
equations at the selected positive A. For n=3,5,7 put

    w_n=n²(n²+1)/[2(n²-1)²], S=sum n²/w_n,
    c_n=2n/(w_n S).

Thus sum n c_n=2. At the endpoint, use
f(cosY)=-cosY+sum c_n cos(nY). The first kernel moment fixes the
coefficient of cosY to -1. The second, divided by A before taking the
endpoint, fixes <f'(cosY)>=1, equivalently sum n c_n=2. The range
velocity energy is sum w_n c_n²=4/S. All these numbers are rational.

At finite A, keep c5,c7 fixed and solve for c1,c3 in

    f(psi)=c1 psi+c3 T3(psi)+c5 T5(psi)+c7 T7(psi),
    <cosY[alpha+f(psi)]>=<cosZ[alpha+f(psi)]>=0.      (1)

Tn are Chebyshev polynomials, not additional Euler modes or prescribed
frequencies. The desingularized moment Jacobian at A=0 is
[[1/2,0],[1/2,3/2]], so this finite polynomial choice continues. More
strongly, the verifier solves(1) at A=1/100 by exact rational algebra.

Define phi=(H-1)^(-1)[alpha+f(psi)] by inversion on the remaining
FINITE Fourier support, and

    z_* = (-alpha/2+Hphi,-phi_Z,phi_Y).

Its zero first shell was imposed before inversion. The full Euler
identity Lz_*+F0=0, including axial velocity and Leray pressure, is
checked directly. Thus z_D=d z_* is an actual stationary first-cell
corrector for EVERY body-frame (kappa,D), d=kappaY DY-kappaZ DZ.

Let R be its range velocity energy with the actual whole-law weight:

    R=(1/5)||z_*, |Fourier wave|>1||².

The rational finite-field calculation gives

    R=636165619333494275181329617802709521563769 /
       112592022217345440000000000000000000000000000,
    0<R<13E/1280, E=1+A².                           (2)

No floating-point stability or spectral threshold selects this field.

## 2. Complete correlated first-shell Euler kernel

Let q_D=kappa cross T_D/lambda, lambda=-1, a_u=kappa.u, and use the
ACTUAL initial material-rate row and mean-current test row

    b=q_D-a_u D+d z_*,
    c=a_u D+(u.D)kappa.

Let Pi_- be the orthogonal curl=-1 projector on the first Fourier shell,
Pi_- v=(P v-curl v)/2 there, and zero elsewhere. Set

    z_return=-Pi_- b+t Pi_- c,
    b_final=b+z_return.                              (3)

The coefficient of each stationary mode depends on the whole laboratory
input expressed in the body frame. This is not a d-only restriction and
is not an independent rotation of the macro wave vector. Since
curl z_return=-z_return, the actual full Euler generator kills it.
Consequently z_D=d z_*+z_return remains a stationary forced corrector.

Direct Fourier/Haar contraction, or orthogonality of the two helicities
on each of the two first-shell waves, gives

    h2=E t²/15-47E/240+R,
    R_D=E(8t+13)/120.                               (4)

Here h2 is the coefficient of the COMPLETE inherited Euler Hamiltonian
H_DD=rho k² h2/2, computed from ||b_final||²-||a_u D||².
The second line is the physical mean acceleration from the full Euler/Lin
current row. These are independent field contractions in the verifier.

For an actual positive restoring coefficient a=-R_D, their difference is

    h2-a=R+(120a-E)(120a+19E)/(960E).                (5)

The two separate stationary preparations that previously retained a
positive mismatch could not access this correlated helicity cancellation.

## 3. Exact positive energy match and its initial phase

Put s=sqrt(100-960R/E) and choose the EXPLICIT preparation

    t=-(4+s)/8, a=E(s-9)/120.                       (6)

The rational inequality(2) gives 19/2<s<10, hence

    E/240<a<E/120, h2=a>0, R_D=-a.                 (7)

This is an algebraic construction parameter fixed by the actual finite
microgeometry, not an empirical fitted modulus or assigned inertia.
The entire field, its corrector, and the current in(4) are evaluated
before using(6).

For an exact finite-K preparation take eta_D(0)=D and its actual velocity
equal to the divergence-compatible translation lift plus ik P_K z_D. Its physical
mean velocity is zero. The common velocity column at its initial point
has eta_V(0)=0 and w_V(0)=V. Equation(1) of `initial-phase-energy.md`
then gives mean pi_D=0, mean pi_V=rho V and the exact initial averaged
phase rho J. The leading mass is rho; the construction has not renamed
any microscopic kinetic norm as a supplied mass.

All fields are real, and the first-K corrector is odd in kappa. The
finite-K family therefore has the proper conjugate Bloch sidebands.
Whole-field translations, rotations, time reversal and mirrors transform
all preparation rows together. On a mirrored realization the helicity
sign changes with curl; Pi_- above is the body-frame projector, not an
unchanged pseudoscalar imposed after reflection.

## 4. Precise scope and the next common-V equation

On the fixed finite cell, the stationary first-cell Euler velocity and
the current formula imply the whole-law physical D mean

    X_D(t)=D-(a/2)k²t² D+O_T(k³),

with its actual D energy matched by(7). Fixed-time spatial jets are meant;
no acoustic-time uniform remainder or all-K invariant manifold is asserted.

The ACTUAL velocity column still obeys

    z_V,t=L z_V+t F_V+b_V,
    b_V=-P[(kappa.u)V+q_V].

Writing y=z_V-t z_D gives

    y_t=L y+b_V-z_D, y(0)=r_V.                       (8)

The isolated acoustic realization would require the WHOLE-LAW observed
current of y to vanish, not necessarily y=0 per orientation. In semigroup
form its exact response is

    E<c, exp(tL)r_V+integral_0^t exp((t-s)L)(b_V-z_D)ds>.

A full coupled realization may instead retain its actual finite angle,
spin and current response here. The required closure then acts on those
physical rows together; the isolated mean need not be autonomous.

A full-field Jordan solution Lr_V=z_D-b_V is one sufficient candidate,
not the only permitted one. The earlier cubic psi^5/psi^7 no-go does not
classify this degree-seven, generally correlated preparation. Nonzero
r_V changes actual gradient mass by rho k² E||r_V||²/2 and introduces
the cross Hamiltonian row rho k² E<b_final,r_V>; both survive elimination.
Nonzero initial material configuration returns change the initial0180
phase one-form as well and cannot be treated as mean-zero velocity returns.

The present two-wave geometry contains no constructed finite EPS optical
packet or physical tag in this artifact. Same-field angle/spin closure,
moving initial-G/current connections, and the common-V response(8) remain
the active parent obligations. The positive result here removes a named
D-energy defect and supplies explicit physical initial data for them.
