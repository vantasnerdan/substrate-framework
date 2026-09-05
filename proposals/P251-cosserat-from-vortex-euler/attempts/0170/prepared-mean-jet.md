# Uniform-period transfer of the actual prepared physical second jet

This is the pressure/current-cancellation candidate registered in0170.
The base is v0.177.0. The field insertion and local EPS/optical construction
are the frozen0153/C-CST-011 inputs at their actual scope. No general
uniform C2 bound for the Bloch Leray operator is claimed or needed.

## 1. Exact object and the three different projections

Let u be any smooth, stationary, mean-zero periodic incompressible Euler
field and p its pressure divided by rho. On its full cell write

    A=u.grad, a=kappa.u, |kappa|=1, K=k kappa,
    P=the microscopic Leray projector at k=0,
    P_K=the full Bloch Leray projector,
    P_kappa=I-kappa tensor kappa.

P preserves harmonic constants; the first cell functions below have
zero mean. On the slow harmonic P_K=P_kappa. It is not differentiated
as though it had a direction-independent value at K=0. All directional
jets are on the physical transverse bundle D,V perpendicular to kappa,
uniformly in kappa. After whole-field isotropy their quadratic symbol
has its usual transverse interpretation.

Prepare the actual material displacement and independent velocity phases
of0154. At zero k, exact Galilean/translation symmetry gives

    eta0=U0=D+tV

for EVERY stationary u, not just the elementary one-wave field. At finite
k use the full Kelvin-D preparation and independently prescribed common-V
circulation class. The exact Jacobi equation is

    eta_tt+2P_K A_K eta_t+P_K(A_K²+Hess p)eta=0,
    div_K eta=0, A_K=A+ik a.                             (1)

It is the variation of the full pressure-constrained material action,
not a selected positive oscillator.

For each selected finite cell, smooth dependence on k exists near0 on
this transverse bundle: nonzero microscopic projectors are smooth in a
neighborhood smaller than that cell's first nonzero wave number. The
constants for higher remainders may depend on this gap. The result below
controls the first two OBSERVED coefficients independently of that gap.

## 2. A gradient identity removes the apparent inverse wave number

The exact Fourier identity, including the harmonic convention, is

    P_K grad f=-ik P_K kappa f.                          (2)

Consequently the first derivative of P_K acting on(Hess p)U0 is
computed by(2), since(Hess p)U0=grad(U0.grad p). It contains no
uncontrolled inverse microscopic wave number.

Write eta=U0+ik chi+O(k²). The first cell is solenoidal, mean zero,
and satisfies

    chi_tt+2PA chi_t+P(A²+Hess p)chi
      =-2P(a V)+P[(kappa.grad p)U0+kappa(U0.grad p)],
    chi(0)=0,
    chi_t(0)=-P[aD+kappa(u.D)].                          (3)

The two pressure terms follow from A a=-kappa.grad p and(2). They
would be lost by differentiating a projector estimate in isolation.
Every P in(3) is order zero; no Poincare or inverse-period estimate has
been used. The initial row is the actual Kelvin row, not fitted data.

There is also no unbounded second-order propagation hypothesis hidden
in(3). Define the actual solenoidal Euler velocity associated with chi,

    w_chi=chi_t+[u,chi], [u,chi]=A chi-(Du)chi.

Stationarity gives A(Du)+(Du)²=-Hess p. Direct substitution turns(3)
into the forced first-order pair

    w_chi,t=-P[A w_chi+(Du)w_chi]+F,
    chi_t+A chi=(Du)chi+w_chi,                           (4)

where F is precisely the right side of(3). Both fields have zero mean.
Their pressure is the full Euler pressure, with P of norm one.

## 3. The second unknown cell drops out of the actual mean

Write eta=U(k,t)+ik chi+k² zeta+O(k³), with U=<eta>, <zeta>=0.
The first mean correction vanishes; the divergence constraint gives

    div zeta=kappa.chi.                                 (5)

Average the UNPROJECTED form of(1), then project with P_kappa. The
unknown pressure multiplier's mean gradient is removed by that slow
projector, not set to zero. In the second coefficient,

    <Hess p zeta>=-<grad p (kappa.chi)>,
    <a A chi>=-<(A a)chi>, A a=-kappa.grad p.             (6)

Thus zeta does not have to be inverted or estimated. The result for
the mean material displacement is

    U_tt=k² P_kappa{
       <a²>U0+2<a chi_t>
       +<(kappa.grad p)chi+grad p(kappa.chi)>}+O(k³).     (7)

The actual Eulerian Fourier velocity mean is different. Lin's formula
v_pert=eta_t+[u,eta], integration by parts and div_K eta=0 give EXACTLY

    m=<v_pert>=U_t+ik<a eta-u(kappa.eta)>.

Therefore its second jet is

    m=U_t-k²<a chi-u(kappa.chi)>+O(k³),
    m_t=k² R_u(t;kappa,D,V)+O(k³),
    R_u=P_kappa{<a²>U0+<a chi_t+u(kappa.chi_t)>
                      +<(kappa.grad p)chi
                                     +grad p(kappa.chi)>}.       (8)

This is the physical observation needed by the action consumer. Replacing
m by U_t would leave the wrong factor2 in the first stress row. Define
the physical displacement X=D+integral_0^t m; its initial velocity is
exactly V. Equations(3),(8), not a bare projector derivative, determine
its complete fixed-time second spatial jet. For constant-pressure
one-wave u0 they reduce to the actual0151 formulas.

## 4. A uniform estimate with no minimum microscopic frequency

Let u=u0+h be the0153 exact constant-curl insertion. Its complete-cell
normalized norms obey

    ||h||_H2=epsilon,
    ||u||_C3+||u0||_C3<=M3,
    ||p-p0||_H2<=C_M3 epsilon, p=-|u|²/2.                (9)

These constants are independent of the quadrature period:0153 bounds
global derivatives by the fixed angular total variation, and all h
frequencies lie on the fixed lambda sphere. The pressure product estimate
uses those uniform C2 bounds and normalized L2 derivatives, not a
period-dependent Sobolev embedding. It does not require ||h||_infinity
to tend to zero. For the elementary wave, grad p0=0.

The full linear Euler group in(4) has L2 and H1 bounds exp(C_M3 T)
independent of cell size. Multiply by the field (and its first spatial
derivatives), integrate the divergence-free transport by parts, and
use the order-zero pressure return; P commutes with spatial derivatives
and has H1 norm one. There is no inverse Laplacian without its gradients.

The one-wave reference chi0,w0 has bounded C2 norms on every fixed time
interval by its explicit0151 formula. Subtract its equation(4) from the
equation for u. The velocity forcing is

    -P[h.grad w0+(Dh)w0]+F_u-F_0,

bounded in H1 by C_T epsilon. The Lin difference has forcing
-h.grad chi0+(Dh)chi0+(w_chi-w0), also bounded in H1. The initial
velocity difference is the order-zero projector applied to
-[(kappa.h)D+kappa(h.D)]; the initial chi difference is zero. Gronwall
for the actual Euler/Lin pair consequently gives

    sup_[0,T](||chi-chi0||_H1+||w_chi-w0||_H1
                         +||chi_t-chi0,t||_L2)
      <=C_T,M3 epsilon (|D|+|V|).                       (10)

The final term follows from the Lin equation; estimating chi in H1 is
what makes that step legitimate. All norm products involving a reference
column use its C2 bound. No small-frequency denominator has entered.

Apply Cauchy--Schwarz to the ACTUAL averaged rows in(8), retaining both
pressure rows. In particular grad p0=0 and
||grad p||_L2<=C_M3 epsilon. Equation(10) gives

    sup_[0,T]|R_u-R_0|
      <=C_T,M3 epsilon (|D|+|V|),                       (11)

uniformly over kappa and phase. Time integration supplies the same bound
for the complete second coefficients of X,m and their needed time
derivatives. First spatial coefficients vanish for the exact mean-zero
preparation; zero-order coefficients are D+tV for both fields.

This estimate includes arbitrary low-frequency insertion self-beats.
They enter through grad p, where the frequency numerator in(2),(6)
is retained, and through the bounded Euler/Lin group in(4). No estimate
of the form period² times epsilon has been erased. The potentially
singular quantity was eliminated from THIS response before it was
bounded. The full Bloch operator may still have unbounded derivatives
as the period increases.

## 5. Selected finite same-EPS cell and the actual action

Freeze the local optical/EPS target and its strict margins from0153 and
C-CST-011, including its optical time window. Choose the angular smoothing
there first; its total variation and global derivative constants are
fixed. Then choose a sufficiently fine rational quadrature so that both
the local C^r error and epsilon in(9) meet the prescribed margins.
Its exact cell identity

    epsilon²=(1+lambda²)² sum_j |a_j|² ->0

permits this choice. NO spectral separation from the dominant wave is
required for(11). The selected finite cell retains the actual optical
object on its own normalization; its per-volume inertia remains positive
at this finite choice, not at a claimed fixed-density infinite-cell limit.

Use one whole-field O(3)/phase/time-reversal law and the same prepared
D,V in every realization. Average the actual phase actions first, as in
0154. Exact initial pairing is rho[[0,I],[-I,0]], since u has zero mean
and the initial Kelvin momentum occupies nonzero microscopic sidebands.
It is not obtained by averaging separately eliminated inertias.

Let f_N,g_N be the actual physical mean rows in either real transverse
polarization after the whole-law average. Their zero jets are1,t, and
their second jets differ from those of the one-wave reference by(11).
The exact physical chart therefore remains

    W_N=f_N g_N,t-g_N f_N,t,
    M_N=rho/W_N,
    K_N=rho(f_N,t g_N,tt-g_N,t f_N,tt)/W_N².              (12)

Put a0=2v0²/15>0. Choose epsilon so that its coefficient bound is smaller
than a0/4. For the selected finite cell and each fixed T, equations(11),
(12) imply

    M_N/rho=1+O_T(epsilon k²)+O_N,T(k³),
    M_N,t/rho=O_T(epsilon k²)+O_N,T(k³),
    K_N/(rho k²)=a0+O_T(epsilon)+O_N,T(k).                (13)

Choose k nonzero LAST, below the selected cell's analytic neighborhood
and its now finite remainder bound. Then the mass and acoustic second-
jet stiffness are strictly positive. The physical/canonical momentum
difference rho(1-1/W_N)m and the full time-dependent connections remain
in the action. The coefficients in(13) may depend on time; stationarity
of u does not turn the prepared mean into a Markov or autonomous closure.

## Scope

This establishes actual, uniformly controlled prepared MEAN/action C2
coefficients and hence a selected finite same-EPS field with a positive
physical acoustic second-jet action and the separately retained optical
margins. It repairs0153's derivative transfer gap by the observed pressure/
current cancellation. It does not establish full-operator C2 convergence,
an analytic radius independent of period, acoustic-time homogenization,
or the parent coupled continuum. Those are not implied by(11)--(13).
