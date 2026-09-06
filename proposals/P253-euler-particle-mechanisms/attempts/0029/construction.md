# P253/0029 — full-label block scheme and principal positive metric

## 1. Strongest result

The scalar inversion used to construct the García--Hassainia--Hmidi orbit
does not transfer to the physical two-label stability space: the two copies
of every principal shape frequency are exactly degenerate.  Retaining those
copies as (2\times2) complex label blocks removes that false divisor.  On
the resulting block decomposition this attempt establishes four exact facts.

1. Every autonomous off-block principal divisor is uniformly separated for
   sufficiently small \(\epsilon\).
2. Every time-dependent principal divisor is transversal in the actual
   source parameter on its zero set; the derivative is written explicitly
   below for all sign cases.
3. For each finite Fourier cutoff there is consequently a nonempty
   finite-stage Cantor set and an exact Hamiltonian label-block homological
   solve.
4. The diagonal high-mode evolution has an explicit uniformly positive
   metric.  With \(\eta=\epsilon^2|c_3|\le3/4\), its weights lie between
   \(1/32\) and \(67/96\).

These are constructive advances, but not the requested all-power theorem.
The published estimates do not define the actual full-label block frequency
shifts, their parameter-Lipschitz decay, or a convergent operator-norm lift of
the finite-step residual.  Consequently neither the finite-cutoff Cantor
sets nor the principal metric can yet be passed to the actual reduced
monodromy.  The annular two-interface source shows exactly why low blocks
cannot be waved away: its Hamiltonian system has hyperbolic low modes, and
the paper avoids them by an \(m\)-fold symmetry that is unavailable here.

## 2. Actual object and inherited scope

Let \(L=|\log\epsilon|\), \(\phi=\omega Lt\), and

\[
 \alpha(\lambda)=\epsilon^2L\omega(\lambda).
                                                               \tag{1}
\]

The actual input from P253/0021 is the full independent-label generator

\[
 A_*(t)=B(t)\partial_\theta+C(t)\mathcal H+R(t)
 \quad\hbox{on}\quad
 X^s=H^s_0(\mathbb T)\oplus H^s_0(\mathbb T),                 \tag{2}
\]

including both first-mode sectors and every self and cross term.  It gives a
bounded invertible one-period map and the invariant flag

\[
 G=\operatorname{span}\{v_{\rm ph},v_z\}
 \subset K=\ker d\mathscr H_c\cap\ker dP,
 \qquad \mathcal N^s_{\rm red}=K/G.                          \tag{3}
\]

The only unconditional transfer from P253/0024 is the labelwise diagonal
principal symbol and the compactness of separately identified
order-minus-three, separated-cross, and finite-center pieces.  The identity

\[
 \sigma_{\rm ess}(M_{\rm red})=\mathbb S^1                  \tag{4}
\]

remains conditional on P253/0024's `CT`; it is not used anywhere below.
In particular, the source residual \(E_{4,n}\) is not silently inserted into
(2) as a compact operator.

All normalized Hamiltonian quantities carry the same physical factor

\[
 (E_{\rm kin}^{\rm phys},I_z^{\rm phys},\Omega^{\rm phys},
   \mathscr A^{\rm phys})
 =2\pi\rho_m(E,P,\Omega,\mathscr A),                         \tag{5}
\]

where material density \(\rho_m\) is distinct from \(\rho=r^2/2\).  The
positive common factor rescales every metric below but cancels from the
Hamiltonian generator, divisor, sign, and condition number.

## 3. The correct label-block homological equation

On \(|j|\ge2\), the published scalar diagonal frequency is

\[
 \mu_j=c_1j+c_2\operatorname{sign}j+\epsilon^2c_3d_j,
 \qquad
 d_j=-\frac{j}{4(j^2-1)|j|}.
                                                               \tag{6}
\]

For the two physical labels, define

\[
 \mathcal B_j
 =\operatorname{span}_{\mathbb C}
   \{e^{ij\theta}\otimes e_1,e^{ij\theta}\otimes e_2\},
 \qquad D_j=i\mu_j I_2.                                     \tag{7}
\]

The real symplectic block is
\(\mathcal W_j=\mathcal B_j\oplus\mathcal B_{-j}\); the
reality condition relates its \(j\) and \(-j\) matrices.  No reversibility
or half-period label relation is imposed.  Write a periodic smoothing term
as \(R=(R^{ab}_{\ell jk})\), or equivalently as \(2\times2\) label matrices
\(R_{\ell jk}:\mathcal B_k\to\mathcal B_j\).

With the convention

\[
 e^{-X}(\alpha\partial_\phi+D+R)e^X
 =\alpha\partial_\phi+D+R+[\alpha\partial_\phi+D,X]+O(XR,X^2),
                                                               \tag{8}
\]

the first homological equation is

\[
 i\Delta_{\ell jk}X_{\ell jk}=-R_{\ell jk},
 \qquad
 \Delta_{\ell jk}=\alpha\ell+\mu_j-\mu_k,
 \qquad
 X_{\ell jk}=i\frac{R_{\ell jk}}{\Delta_{\ell jk}}.        \tag{9}
\]

Thus the sign in (9) cancels the off-block term in (8).  Hamiltonian
projection is taken together on \((\ell,j,k)\) and its real/KKS partner, so a
Hamiltonian \(R\) produces a Hamiltonian \(X\), and \(e^X\) is symplectic.

For \((\ell,j,k)=(0,j,j)\), (9) is never used.  The complete matrix
\(R_{0,jj}\) is retained in \(\mathcal W_j\).  At a later iteration the
retained generator is a Hamiltonian matrix \(D_j^{(m)}\), and the correct
off-block equation is the Sylvester equation

\[
 i\alpha\ell X_{\ell jk}^{(m)}
 +D_j^{(m)}X_{\ell jk}^{(m)}
 -X_{\ell jk}^{(m)}D_k^{(m)}=-R_{\ell jk}^{(m)}.             \tag{10}
\]

If its block eigenfrequencies are \(\nu^{(m)}_{j,p}\), equation (10) asks
for the same-branch and cross-branch divisors

\[
 \alpha\ell+\nu^{(m)}_{j,p}-\nu^{(m)}_{k,q},
 \qquad p,q\in\{1,2\}.                                    \tag{11}
\]

This is precisely the useful structural transfer from the annular
two-interface Proposition 7.5.  That proposition does not transfer its
estimates: its interfaces, transport speeds, external parameters, symmetry
class, and residual are different.

The first/center modes and the energy/impulse companions are not hidden in
(7).  They form a finite Hamiltonian block \(\mathcal F\).  A coupling to
\(\mathcal F\) is inverted only after its own Sylvester separation is proved;
otherwise it is joined to the relevant retained cluster.  The quotient (3)
is applied only after this construction.

## 4. Exact autonomous gap calculation

The source bounds are

\[
 -\frac23\le c_1\le-\frac13,
 \qquad \frac13\le c_2\le\frac{13}{24},
 \qquad |d_j|\le\frac1{12}.                                \tag{12}
\]

Put \(n=j-k\), \(s=\operatorname{sign}j-\operatorname{sign}k\),
and \(e=d_j-d_k\).  If \(j,k\) have the same sign and \(j\ne k\),

\[
 |\mu_j-\mu_k|
 \ge \frac13|n|-\frac16\epsilon^2|c_3|.                   \tag{13}
\]

Hence \(\eta:=\epsilon^2|c_3|\le1/2\) gives
\(|\mu_j-\mu_k|\ge1/4\).  If the signs are opposite, then
\(|n|\ge4\), \(|s|=2\), and \(n,s\) have the same sign.  Therefore

\[
 |\mu_j-\mu_k|
 \ge \frac13|n|-\frac{13}{12}-\frac16\eta
 \ge \frac14-\frac16\eta.                                \tag{14}
\]

For \(\eta\le1/2\), (14) is at least \(1/6\).  Equations (13)-(14)
prove that the only autonomous principal degeneracy is

\[
 j=k,\quad \ell=0,
\]

including the cross-label copy already retained inside \(\mathcal B_j\).
They do not classify the perturbed retained matrix.

## 5. Exact principal transversality

For every principal pair define

\[
 F_{\ell jk}(\lambda)
 =\alpha(\lambda)\ell+c_1(\lambda)n+c_2(\lambda)s
       +\epsilon^2c_3(\lambda)e.                            \tag{15}
\]

Its derivative is exactly

\[
 F'_{\ell jk}
 =\alpha'\ell+c_1'n+c_2's+\epsilon^2c_3'e.                 \tag{16}
\]

On a zero of (15), put \(\beta=\omega'/\omega=\alpha'/\alpha\).
Eliminating \(\ell\) gives the more useful identity

\[
 F'_{\ell jk}
 =n(c_1'-\beta c_1)+s(c_2'-\beta c_2)
   +\epsilon^2e(c_3'-\beta c_3).                            \tag{17}
\]

The source proof has \(\inf|\beta|>0\) on its compact parameter interval,
and its coefficient expansions have

\[
 c_1=-\frac12+r_1,\qquad c_2=\frac12+r_2,\qquad
 \|r_1\|_{C^1}=O(\epsilon^2L^{1/2}),\qquad
 \|r_2\|_{C^1}=O(\epsilon^2),\qquad \|c_3\|_{C^1}=O(1).
                                                               \tag{18}
\]

Substitution into (17) yields

\[
 F'_{\ell jk}
 =\frac{\beta}{2}(n-s)
 +n(r_1'-\beta r_1)+s(r_2'-\beta r_2)
 +\epsilon^2e(c_3'-\beta c_3).                             \tag{19}
\]

This resolves every frozen sign case:

| case | \(s\) | exact leading derivative on resonance | conclusion for small \(\epsilon\) |
|---|---:|---:|---|
| same sign, \(j\ne k\) | 0 | \(\beta n/2\) | \(|F'|\ge C|n|\) |
| opposite signs | \(\pm2\) | \(\beta(n-s)/2\) | \(|n|\ge4\), \(|n-s|\ge2\), hence \(|F'|\ge C\max(1,|n|)\) |
| same mode, different label, \(\ell=0\) | 0 | 0 | average retained in \(\mathcal B_j\), not divided |
| same mode, different label, \(\ell\ne0\) | 0 | no zero since \(F=\alpha\ell\) | invertible temporal harmonic |
| \(\ell=0\), off block | any | not needed | separated by (13)-(14) |
| \(\ell\ne0\), off block | any | (19) | transversal on every principal resonance |

Here \(C>0\) is obtained by choosing \(\epsilon\) so the explicit remainder
in (19) is at most one half of its leading term.  This uses the actual source
coefficient asymptotics, not a generic parameter assumption.  It is a new
pairwise principal calculation; the source Lemma 10.1 proves only the scalar
first-Melnikov analogue.

## 6. Finite-cutoff Cantor construction and its exact limit obstruction

For a finite cutoff \(N\), let \(\mathfrak I_N\) contain all selected
off-block and nonzero internal-time tuples \((\ell,j,k,p,q)\) with
\(|\ell|,|j|,|k|\le N\), and define

\[
 \mathcal C_N(\gamma)=
 \bigcap_{\mathfrak I_N}
 \left\{\lambda:\
 |\alpha\ell+\nu_{j,p}-\nu_{k,q}|
 \ge\gamma\langle\ell\rangle^{-\tau_t}
               \langle j-k\rangle^{-\tau_x}\right\}.      \tag{20}
\]

At the principal step \(\nu_{j,p}=\mu_j\), (19) and the one-dimensional
sublevel estimate give, resonance by resonance,

\[
 |\mathcal R_{\ell jk}|
 \le \frac{2\gamma}{C\max(1,|j-k|)}
       \langle\ell\rangle^{-\tau_t}
       \langle j-k\rangle^{-\tau_x}.                       \tag{21}
\]

Because \(\mathfrak I_N\) is finite, (21) makes
\(\mathcal C_N(\gamma)\) nonempty for all sufficiently small \(\gamma\),
and (9) or (10) gives an exact finite-cutoff symplectic homological solve
there.  This is stronger than a bibliography transfer: it constructs the
correct physical-label first reduction step and removes every selected
nonresonant coefficient to that order without ever dividing by the label
degeneracy.  Higher commutators remain for the next finite iteration.

The limit \(N\to\infty\) is not licensed.  Same-sign affine frequencies
have infinitely many pairs with the same difference \(n=j-k\); the small
terms \(d_j-d_k\) cluster at zero.  A convergent proof therefore needs
actual iterated block shifts \(r^{(m)}_{j,p}\) satisfying, on the full label
space, estimates of the form

\[
 |r^{(m)}_{j,p}|+|\partial_\lambda r^{(m)}_{j,p}|
 \le C_m\langle j\rangle^{-\sigma},
 \qquad \sum_m C_m<\infty,                                 \tag{22}
\]

together with hybrid off-diagonal smoothing sufficient to sum the clustered
resonance widths.  The annular paper proves its own analogue using two
external parameters, distinct interface asymptotics, and \(m\)-fold spaces.
The García paper gives neither (22) for two independent labels nor a
full-label sequence in which \(E_{4,n}\) vanishes in operator norm.  A naive
pair count in (21) grows with \(N\), so finite-stage nonemptiness is not a
substitute for this missing estimate.

## 7. Exact shape-block KKS matrix

The source shape coordinate does allow the KKS matrix to be computed
exactly.  For one label,

\[
 \gamma=P+i\epsilon L^{-1}V+\epsilon wZ,\qquad
 w^2=1+2\epsilon f,\qquad \det(Z,Z_\theta)=1.                \tag{23}
\]

At fixed center, a shape variation satisfies

\[
 \delta\gamma=\epsilon^2w^{-1}\delta f\,Z,\qquad
 (\delta\gamma\cdot n)\,ds
 =\det(\delta\gamma,\gamma_\theta)d\theta
 =\epsilon^3\delta f\,d\theta.                              \tag{24}
\]

If \(\chi\) generates the normal displacement \(k=\delta\gamma\cdot n\),
then \(k=-\partial_s\chi\), hence
\(\partial_\theta\chi=-\epsilon^3\delta f\).  Substitution in the normalized
product-leaf KKS form gives, on zero-mean shape variables,

\[
 \Omega_f(u,v)
 =-\epsilon^4\sum_{a=1}^2
   \int_{\mathbb T}(\partial_\theta^{-1}u_a)v_a\,d\theta.
                                                               \tag{25}
\]

For \(j>0\), use the \(L^2\)-normalized real basis
\((\pi^{-1/2}\cos j\theta,\pi^{-1/2}\sin j\theta)\) on each label and put
\(J=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)\).
Then the exact normalized and physical KKS matrices are

\[
 W_j=-\frac{\epsilon^4}{j}\operatorname{diag}(J,J),
 \qquad W_j^{\rm phys}=2\pi\rho_m W_j.                      \tag{26}
\]

For the strict untransformed leading self operator, put
\(\mu_j^{\rm lead}=-j/2+\operatorname{sign}(j)/2\).  Its physical-time
generator is
\[
 A_{0,j}^{\rm lead}
 =-\epsilon^{-2}\mu_j^{\rm lead}\operatorname{diag}(J,J).
\]
With \(i_X\Omega=dH\), so that \(A=-W^{-1}L_H\), its normalized
Hamiltonian Hessian is therefore

\[
 L_{H,0,j}^{\rm lead}
 =\epsilon^2\frac{\mu_j^{\rm lead}}{j}I_4.                 \tag{27}
\]

For \(j\ge2\),
\(-1/2<\mu_j^{\rm lead}/j\le-1/4\), so all four real directions of
this strict leading block have the same negative energy sign.
Multiplication of both (26) and (27) by \(2\pi\rho_m\) leaves
\(A_{0,j}^{\rm lead}\) unchanged.

This exact Hessian statement is deliberately not extended to the final
multiplier \(\mu_j\).  García Proposition 8.3 calls its maps bounded,
reversibility-preserving isomorphisms on the restricted periodic space; it
does not prove that their full-label lifts preserve (25).  Therefore
\(-\mu_j/j\) below is an exact positive metric for the transformed diagonal
model, but it is not identified with the actual orbit Hessian without a
separate symplectic-lift calculation.

## 8. Route B executed: an explicit positive principal metric

The physical-time diagonal generator on each label copy is

\[
 A_{0,j}=-i\epsilon^{-2}\mu_j I_2.                          \tag{28}
\]

Define the diagonal weight

\[
 g_j=-\frac{\mu_j}{j}
 =-c_1-\frac{c_2}{|j|}-\epsilon^2c_3\frac{d_j}{j},
 \qquad |j|\ge2.                                           \tag{29}
\]

Since \(|d_j/j|\le1/24\), equations (12) and
\(\eta=\epsilon^2|c_3|\le3/4\) give

\[
 \boxed{\frac1{32}\le g_j\le\frac{67}{96}}.               \tag{30}
\]

Let \(G_0\) have Fourier blocks \(g_jI_2\), with the finite sector assigned
any fixed positive metric.  Then, exactly,

\[
 A_0^*G_0+G_0A_0=0,
 \qquad
 \frac1{32}I\le G_0\le\frac{67}{96}I                      \tag{31}
\]

on the high-mode sector.  The physical metric is
\(G_0^{\rm phys}=2\pi\rho_mG_0\); its high-mode condition number is at most
\(\sqrt{67/3}\).  This is a genuine positive conserved norm for the full
two-label diagonal evolution, not a claim about the actual remainder.

For \(A=A_0+\mathcal R\), write \(G=G_0+K\).  The exact periodic metric
equation becomes

\[
 \dot K+A_0^*K+KA_0
 =-(\mathcal R^*G_0+G_0\mathcal R)
   -(\mathcal R^*K+K\mathcal R).                            \tag{32}
\]

Its Fourier homological operator has the same divisor as (9), up to an
overall sign.  The \((0,j,j)\) component is instead a retained block
Lyapunov equation.  Thus Route B is genuinely constructed through its
positive seed and exact nonlinear metric equation, but it does not evade the
pairwise divisor or the unresolved block monodromy.  Small norm of
\(\mathcal R\) alone cannot solve (32), because a periodic Hamiltonian
perturbation can have a parametric resonance.

## 9. Route C executed: actual Krein domain and the source obstruction

Restrict the normalized KKS form to \(\mathcal W_j\), denote its real
skew matrix by \(W_j\), and let \(L_j(t)=L_j(t)^T\) be the normalized
Hamiltonian Hessian there.  With the convention \(i_X\Omega=dH\), the
physical factor multiplies both \(W_j\) and \(L_j\) and the actual retained
Jacobi block has

\[
 A_j(t)=-W_j^{-1}L_j(t),
 \qquad A_j(t)^TW_j+W_jA_j(t)=0,                            \tag{33}
\]

Equation (33), not label exchange or reversibility, is the correct finite
Krein problem.  Equations (23)-(27) give the exact KKS and strict-leading
Hessian matrices, while (29)-(31) give a positive transformed-diagonal
metric.  The strict leading block has one Krein sign and a semisimple
unitary return.  For the actual retained
block, however, one must compute the period map \(M_j\) and either find

\[
 M_j^TG_jM_j=G_j,\qquad 0<mI\le G_j\le MI,               \tag{34}
\]

uniformly in \(j\), or exhibit its reciprocal off-circle pair or unit Jordan
chain.  Positivity at each instant is not enough unless the metric also
satisfies (32).

The annular source gives a concrete warning and a concrete model for (34).
Its exact \(2\times2\) Fourier matrix has a discriminant \(\Delta_j(b)\):
some low blocks are hyperbolic, while high blocks are elliptic and admit a
symplectic diagonalizer \(Q_j\) because
\(0\le a_j(b)<1-\delta\).  It then discards the hyperbolic modes by imposing
\(m\)-fold symmetry with \(m\ge m_*(b_*)\).  P253's physical perturbation
class contains all modes, so that exclusion cannot be imported.  The actual
ring matrices \(L_j(t)\), especially the finite center/first-mode clusters,
must be calculated from (2); nontranslation of the base orbit supplies no
Krein sign.

## 10. Exact conclusions and next construction

The following ladder has now been executed:

1. scalar inversion was replaced by exact label blocks;
2. all autonomous and temporal principal divisors were calculated;
3. finite-cutoff symplectic block reduction was constructed;
4. a materially different positive-metric route was constructed at principal
   level and reduced to (32); and
5. the remaining resonances were formulated on their actual KKS/Krein
   blocks, with the published low-mode obstruction kept in scope.

The next executable construction is now sharply localized.  Compute from the
actual two-contour Hessian the matrices \(L_j(t)\) for the first/center block
and for \(2\le |j|\le J_0\), while proving a full-label version of (22) for
\(|j|>J_0\).  Then either:

- solve (10) iteratively with a summable clustered-resonance measure bound
  and an operator-norm-zero final residual, which would also earn `CT`; or
- solve the periodic block Lyapunov equations (32), (34), retaining every
  finite block whose Krein form is mixed.

Until that calculation exists, P253/0021's map remains only bounded and
invertible for one period.  Nothing here licenses (4), reduced all-power
boundedness, nonlinear or global Euler stability, swirl or helicity, quantum
spin, an electron or neutrino identification, parent completion, or a global
no-go.
