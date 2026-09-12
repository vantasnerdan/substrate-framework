# A shared compact-vorticity column for exact waves and spectral control

## Primary theorem and its exact scope

Gallay--Smets, arXiv1805.05064v3, Theorem1.3 and Remark1.9, establish imaginary
spectrum for each azimuthal integer m and nonzero axial wavenumber k in their
solenoidal enstrophy space X_mk. Remark1.9 extends this to the C_b^1 closure
of their strict profile class. The class requires decreasing positive
vorticity W, finite circulation, and decreasing J=2Omega W/Omega'^2, with
the stated endpoint regularity. This is fixed-mode spectral control. Their
separate1811.07584v2 energy-space theorem gives subexponential full-group
bounds for its strict class; no uniform limiting constants are supplied here.
Neither theorem is a solitary-wave stability theorem.

Source: https://arxiv.org/pdf/1805.05064v3 , PDF SHA256
081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797. The primary body was opened after the0030 source inventory.
Author-PDF timeout preceded successful versioned arXiv access. No failed
access is counted as theorem evidence.

## Direct construction inside the closure and the0027 existence class

Fix0<a<b and Omega_0>0. Define on(a,b)

    J(r)=exp[1/(r-a)-1/(b-r)],
    alpha(r)=4/[1+sqrt(1+4J(r)/r^2)].                  (1)

Extend alpha=0 forr<=a and alpha=2 forr>=b. At the inner endpoint alpha
is a smooth flat multiple of exp[-1/(2(r-a))]; at the outer endpoint
2-alpha is a smooth flat multiple of exp[-1/(b-r)]. Consequently

    Omega(r)=Omega_0 exp[-int_0^r alpha(s)/s ds],
    W(r)=Omega(r)[2-alpha(r)]                         (2)

are smooth radial functions. W=2Omega_0 near the axis and W=0 forr>=b,
with all derivatives matching at both joins. On the annulus J'<0. Also
J/r^2 decreases strictly, so alpha'>0. Therefore

    W'=Omega[-alpha(2-alpha)/r-alpha']<0               (3)

there. The column is precisely a smooth nonnegative compact-vorticity
pure-swirl background admitted by0027. Its circulation function is
L=r^2 Omega, and L'=r W. In the transition annulus,

    2Omega W/Omega'^2
       =2r^2(2-alpha)/alpha^2=J.                     (4)

Thus the auxiliary function is obtained from the actual field rather than
assigned independently after choosing it.

## Explicit strict approximants and C_b^1 convergence

For epsilon>0 define

    J_e=J/(1+epsilon r^2 J)+epsilon/r^2,              (5)

using1/(epsilon r^2) for the first term whereJ=+infinity onr<=a,
and0 whereJ=0 onr>=b. Let alpha_e,Omega_e,W_e follow(1)-(2) withJ_e.
The reciprocal of the first term is1/J+epsilon r^2, strictly increasing;
its derivative is negative, as is that ofepsilon/r^2. ThusJ_e'<0
throughout r>0. Moreover J_e/r^2 strictly decreases, alpha_e'>0,
and W_e'<0. Near0, J_e=(epsilon^-1+epsilon)/r^2, soalpha_e is a
smooth function ofr^2 withalpha_e=O(r^2), andW_e'(0)=0. Forr>=b,

    J_e=epsilon/r^2, 2-alpha_e=O(epsilon/r^4),
    Omega_e=O(r^-2), W_e=O(epsilon/r^6).              (6)

The derivative bounds from these explicit formulas give r^3 W_e'->0,
r J_e'->0 and finiteint r W_e dr. After choosingOmega_0=1, W_e(0)=2,
so every strict profile meets the source normalization and hypotheses.

For completeness, the convergence is in the required topology, not merely
pointwise. Away from a,b it follows by differentiating(1),(5). Near a use
the flat coordinate q=J^-1/2, extended as0 forr<=a. The reciprocal square
root ofJ_e is

    q_e=sqrt(q^2+epsilon r^2)
         /sqrt(1+epsilon q^2/r^2+epsilon^2).           (7)

On a fixed neighborhood ofa, q_e->q uniformly. In its derivative, the
only nonuniform-looking factor isq/sqrt(q^2+epsilon r^2) multiplyingq'.
Split intoq>=delta andq<delta: the first region converges uniformly for
fixed delta, while sup|q'| on the second tends to0 asdelta->0 by flatness
and monotonicity ofq near a. The remaining derivative contribution is
O(sqrt(epsilon)). Thusq_e->q inC^1 there, andalpha is smooth inq atq=0.
Nearb, J is itself smooth and flat, and(5) converges inC^1; the alpha map
is smooth inJ at0. On0<=r<=a/2 the explicit inner expression yields
alpha_e=O(sqrt(epsilon) r^2), alpha_e'=O(sqrt(epsilon) r).
Finally(6) and its differentiated version control the whole tailr>=b.
Integrating(alpha_e-alpha)/r and using(2)-(3) therefore proves

    ||W_e-W||_infinity+||W_e'-W'||_infinity ->0.       (8)

This explicitly places one complete family of0027 backgrounds in the
source's closure class. It does not assert that every smooth compact
profile belongs to that class.

## Earned transfer and next actual object

For every column(1)-(2), and every fixedm inZ,k!=0, the full pressure/Hodge
vorticity linearization has no spectrum off the imaginary axis inX_mk.
Together with the direct0030 positive axisymmetric physical metric, this
provides a background on which the exact solitary existence construction
and controlled background dynamics have compatible hypotheses.

The remaining solitary problem concernsU+v_c, notU. The fixed-mode
spectral statement supplies neither a uniform resolvent bound as(m,k)
vary nor robustness of the localized excitation. The exact wave breaks
axial translation invariance of its coefficients, coupling axial Fourier
modes. Its full Leray operator and zero-frequency modulation must therefore
be retained in the next propagation argument.
