# Retained mathematics, with hypotheses and limits

This note extracts existing arguments from the stopped P256 campaign. It is an archival mathematical harvest, **not** a new physical theory, a priority claim, or scientific acceptance of the campaign. The uniform orbital model and the later determinant/clover transfer are different constructions; their results cannot discharge one another's obligations by sharing terminology. No new numerical calculation or scientific replay was performed for this publication.

## 1. Exact uniform-field Gram/Clifford reduction

Let \(A_0,\ldots,A_3\) be arbitrary **uniform**, traceless Hermitian \(4\times4\) matrices, and let \(x_p\in\mathbb R^4\) be an **equal-weight** family of \(N\) prescribed scalar momentum coefficients. Set

$$
G_{\mu\nu}=\tfrac14\operatorname{tr}(A_\mu A_\nu),\qquad
h_p=\sum_\mu x_{p\mu}A_\mu,\qquad
V(G)=\tfrac\kappa2\|G-g_0I\|_F^2,
\quad \kappa>0,\ g_0\ge0.
$$

The fermion energy is the actual lowest energy with **total occupation \(2N\)** among all \(4N\) one-body levels, divided by \(N\). Two occupied states per momentum are not imposed on competitors. In the uniform identity-link sector,

$$
\min_{A}\mathcal E(A)
=\min_{G\succeq0}F(G),\qquad
F(G)=\tfrac\kappa2\|G-g_0I\|_F^2
-\frac2N\sum_p\sqrt{x_p^TGx_p}.
$$

**Lower bound.** If a traceless four-level block has \(k\) negative eigenvalues with sum \(-S\), zero trace and Cauchy–Schwarz give

$$
\operatorname{tr}h^2\ge S^2/k+S^2/(4-k),\qquad
S\le\sqrt{k(4-k)\operatorname{tr}h^2/4}
\le\sqrt{\operatorname{tr}h^2}.
$$

The zero block is immediate. Any occupied subset has energy at least the sum of all negative levels. Since \(\operatorname{tr}h_p^2=4x_p^TGx_p\), summing proves \(\mathcal E(A)\ge F(G)\), including occupation redistributed between momenta.

**Attainment.** For every PSD \(G=LL^T\), take four Hermitian Clifford generators \(\Gamma_\nu\) in dimension four and set \(A_\mu=\sum_\nu L_{\mu\nu}\Gamma_\nu\). Then

$$
\{A_\mu,A_\nu\}=2G_{\mu\nu}I,\qquad
h_p^2=(x_p^TGx_p)I.
$$

Each nonzero block has two negative and two positive levels; zero modes can complete the required half filling. These witnesses attain the bound for every \(G\). They are witnesses, not restrictions on admissible competitors.

**Selection.** The negative square-root terms are convex; the elastic term makes \(F\) coercive and \(\kappa\)-strongly convex. Hence its Gram minimizer \(G_*\) is unique. If the coefficient vectors span \(\mathbb R^4\), a putative null direction \(z\) is excluded by adding \(\epsilon zz^T\): the sea decreases strictly at order \(\epsilon\) or \(\sqrt\epsilon\), while the elastic change is \(-\kappa g_0\epsilon+O(\epsilon^2)\). Thus \(G_*>0\), also at \(g_0=0\). Cubic sign/permutation symmetry, when present in the supplied family, forces

$$
G_*=\operatorname{diag}(g_M,g_T,g_T,g_T),\qquad g_M,g_T>0,
\qquad
\mathcal E(A)-\mathcal E_*\ge\tfrac\kappa2\|G-G_*\|_F^2.
$$

For the stronger claim that **every** matrix minimizer is Clifford, spanning alone is insufficient. The ten-column quadratic design \((x_{p\mu}x_{p\nu})_{\mu\le\nu}\) must have rank ten. Equality in the spectral bound then makes \(h_p^2=(x_p^TG_*x_p)I\) at every momentum; full design rank forces each square/anticommutator coefficient to vanish. For the supplied even cubic mass law, the ten momenta

$$
0,\quad\pm ke_1,\pm ke_2,\pm ke_3,\quad
k(e_1+e_2),k(e_1+e_3),k(e_2+e_3)
$$

have quadratic-design determinant magnitude \(8m^2|\mu|^3|t|^{15}\), where \(m=\mathrm{mass}(0)\), \(\mu=\mathrm{mass}(ke_i)\), \(t=v\sin(ak)/a\), provided \(m\mu t\ne0\).

Sources: [soft-budget derivation, §§1–3](soft_budget.md), [rank-ten proof, §3](sea_selection.md), and [uniform periodic global-occupation transfer](bulk_global_gram.md). Unequal quadrature weights require their own occupation interpretation. Uniformity, scalar hopping, four orbitals, exclusion/half filling and the chosen Gram-only potential are supplied. This selects an **orbital algebra**, not spacetime, physical spin, Lorentz invariance, a fluctuating quantum ordered phase or a propagating particle. Historical effective spin/current language in those sources does not remove these limitations.

## 2. A fixed output can prevent coherent cancellation

Pointwise invertibility and positive determinants do not imply a nonzero coherent integral. The retained source-specific argument instead uses one output functional independent of the integrated hidden field.

At fixed output field \(W=I\), the original pretransfer source and both actual endpoint masses are retained. A fixed invertible output transformation, defined using a midpoint reference root \(h_r\), gives

$$
P_r^{-1}=F_U L_U^{-1},\qquad
F_U=\ell_{*,U}^{-1}\ell_u^{-1},\qquad
L_U=P_-h_r^{-1}h_w[Y_u;I],\qquad
w_U=c_wc_u\det P_r>0.
$$

The reference changes the **output test**, not the source, preparation, actual masses or field-dependent weight. A sufficient uniform hypothesis is

$$
\operatorname{Re}(e^\dagger F_Uv)\big|_{v=e}
\ge c\|e\|\|F_U^\dagger e\|,\qquad
\sigma_{\min}(F_U)\ge s>0,\qquad
\|L_U-I\|\le r<c\le1.
$$

Here the angle hypothesis is explicitly for \(F_U^\dagger e\); the retained native argument bounds both adjoint orientations. Put \(v=L_U^{-1}e\). The angle between \(v\) and \(e\), in the real Hilbert inner product, is at most \(\arcsin r\), and \(\|v\|\ge\|e\|/(1+r)\). Combining the two angle bounds yields the correctly normalized inequality

$$
\operatorname{Re}(e^\dagger F_UL_U^{-1}e)
\ge\frac{s}{1+r}
\left[c\sqrt{1-r^2}-r\sqrt{1-c^2}\right]\|e\|^2>0.
$$

The same functional therefore separates every integrand from zero and survives integration against \(w_U\,dU\). This is stronger than pointwise nonvanishing and does not replace a complex physical amplitude by its modulus.

### Why a positive Gramian is legitimate here

The retained correlated seam reduction has \(R(t)=(g_r+tI)^{-1}\) and a Hermitian, generally **signed** small matrix \(K(t)\), with

$$
h_M-h_r=-\frac1\pi\int_0^\infty\sqrt t\,R(t)VK(t)V^\dagger R(t)\,dt,
\qquad \|K(t)\|\le k.
$$

Introduce only the mathematical Gramian

$$
H=\frac1\pi\int_0^\infty\sqrt t\,R(t)VV^\dagger R(t)\,dt\succeq0,
\qquad h_rH+Hh_r=VV^\dagger.
$$

Then \(-kH\preceq h_M-h_r\preceq kH\). Its positive blocks and Cauchy–Schwarz, together with the retained root symmetry, give the specific row bound

$$
\|L_U-I\|\le k\sqrt{A_+}\bigl(\sqrt{A_-}+y\sqrt{A_+}\bigr),
\quad
A_\pm=\frac1\pi\int_0^\infty\sqrt t\,\|P_\pm R(t)V\|^2dt,
\quad \|Y_u\|<y.
$$

The physical signed seam amplitude has not been declared positive. Full-height boundary Schur composition retains all 178 layers while reducing the fixed reference matrices to size at most 288; no layer truncation is licensed by the reduction.

For the specified open \(4\times3\times3\) construction, retained native/root bounds give \(c>0.258\), \(s>1/14\), and the corrected conservative \(r<31/200\), implying a transformed cofactor real part greater than \(3/500\) for unit \(e\). **These are conditional inherited bounds, not an enclosure replay in this publication.** Their native source/Hodge identification, complete root bounds, stated IEEE-operation error model and outward scalar arithmetic are load-bearing; the assembled reproducer was not independently rerun here.

Applying the same estimate on the two original raw-spin coordinates gives a strictly accretive \(2\times2\) integrated output matrix. It is invertible at \(W=I\); continuity gives a positive-measure neighborhood and positive coherent norm for every nonzero positive density on those two original sources. There is no useful global residue bound or detector/field covariance here. Independence of midpoint-multiplied products does not prove independence of propagated input-multiplied products. This is not a photon, scattering or particle theorem. Operator-angle and square-root/Gramian methods are established mathematics: see [Gustafson](http://w2.ict.nsc.ru/jct/content/t4n3/Gustafson.pdf) and [Del Moral–Niclas](https://arxiv.org/abs/1705.08561).

## 3. Conditional fixed-box neutral/charged sector separation

The retained ideal-law proof concerns a compact gauge-field transfer \(T=A^\dagger A\), its physical neutral Ref0 sector, and its charged sector with **one fixed scalar reference hole**. It is conditional on the native chart/domain identities, a proper acute generating middle-exterior cone and the original-source survival premises above. It does not concern an arbitrary transfer or the accuracy of a stored compressed-root implementation.

In the fixed chart, let \(h=g^{1/2}>0\), \(h\Sigma h=\Sigma\), \(g_{++}=B^{-1}\), and

$$
B>b_0>41/40,\qquad \gamma^{-1}I\preceq g\preceq\gamma I,
\quad\gamma=515/12,\qquad J_c=\operatorname{diag}(-I,I).
$$

Writing \(F=\log h=\begin{psmallmatrix}X&Y\\Y&-X\end{psmallmatrix}\), compression Jensen for operator-concave log gives

$$
X=\tfrac12(\log g)_{++}\preceq\tfrac12\log(g_{++})
\preceq-\tfrac12\log b_0\,I.
$$

Congruence and integration then yield

$$
h^\dagger J_ch-J_c
=\int_0^1e^{sF}(FJ_c+J_cF)e^{sF}ds
\succeq\log b_0\,\frac{1-\gamma^{-1}}{\log\gamma}I
>\frac1{200}I.
$$

Unitary centered gauges preserve \(J_c\). Hence an arbitrary \(m\)-root word \(H_m\) expands the full minus subspace by at least \((201/200)^{m/2}\). Reciprocal singular-value pairing gives

$$
\|\Lambda^{N+1}H_m\|
\le(200/201)^{m/2}\|\Lambda^NH_m\|.
$$

The neutral cone argument supplies a fixed interior test Slater \(e\) and one finite constant \(C_N\), independent of fields and word length, such that each centered neutral word obeys

$$
\|L_0\|\le C_N\langle e,L_0e\rangle,
\qquad \langle e,L_0e\rangle\ge0.
$$

The scalar reference hole adds only unit-modulus phases in the charged histories. Complete unreset histories in \(T^n\) contain \(4n\) roots. Applying the triangle inequality **only as an upper bound** on the charged coherent integral, and comparing the neutral scalar kernel with the isometric compression of \(T_0^n\), gives

$$
\|T_{1,r}^n\|\le C_N\bigl[(200/201)^2\lambda_0\bigr]^n,
\qquad
0<\lambda_{1,r}\le(200/201)^2\lambda_0.
$$

Taking \(n\)-th roots removes \(C_N\). Strict charged nonzero support uses the separate source-survival argument, not neutral positivity. With **blocked clock \(2a_t\)**,

$$
E_{1,r}-E_0\ge\log(201/200)/a_t.
$$

The retained neutral proof uses strict graph-ball invariance, centered determinant characters and an interior cone sandwich to obtain a unique neutral top level with nonzero original ideal-trial overlap. These are substantive premises: their complete inherited proof chain has not received a fresh independent audit here. A correctly bound antiunitary squaring to \(-1\) on the odd charged CAR grade gives Kramers degeneracy, **not physical spatial spin**.

The strict \(B>1\) margin depends on the selected open box; no uniform-in-volume or continuum result follows. The static reference, absent charged dispersion, missing field dictionary and lack of guaranteed source overlap with the charged top prevent an electron-mass interpretation. See [Hansen–Pedersen](https://arxiv.org/abs/math/0204049) for the established Jensen tool.

## 4. Exact adverse examples worth retaining

### Fourier truncation does not inherit angle positivity

For one rotor with noncyclic Fourier cutoff \(R=1\), the free heat kernel at separation \(\pi\) is

$$
k_{1,t}(\pi)=\frac{1-2e^{-tU_E/2}}{2\pi}<0
\quad(0<t<2\log2/U_E).
$$

This invalidates importing untruncated angle-space positivity. The stronger two-adjacent-plaquette witness in [the existing derivation](photon_neutrality.md#parent-executed-smallest-cutoff-counterexample) uses \(R=1\), \(U_E=1/100\), complete charge \(q=(0,-1,0,0,1,0)\), and a stated real multi-winding Wilson potential. Its neutral/charged dimensions are 7/6. The retained rational data give

$$
H_{0,1}-\lambda I\succ0,\quad
\lambda=-1771569/500000,\qquad
\mathcal R_q=-28043796471897/7692314248600<\lambda.
$$

Positive exact LDL pivots prove the first inequality; the explicitly recorded rational charged trial proves the second. Thus \(E_{q,1}<E_{0,1}\). This refutes the general **cutoff-neutral** comparison, not the untruncated theorem or the positive fundamental-plaquette subclass. The linked source supplies graph, fields, coefficients and trial vector; no calculation was repeated here.

### A complete charged half-kernel need not be reflection-Gram

For the retained native pencil

$$
g(B,C)=\begin{pmatrix}B^{-1}&-B^{-1}C\\-C^\dagger B^{-1}&B+C^\dagger B^{-1}C\end{pmatrix},
\quad L=\det B\,\Lambda^{n+1}g,
\quad Jv=\Omega_+\wedge v,
$$

where \(\Omega_+\) fills all \(n\) plus coordinates, bordered-minor Schur algebra gives **exactly** \(J^\dagger LJ=B\). The \(C\) terms cancel by the Schur complement, not by omission.

The existing two-endpoint construction sets one SU(3) edge to \(U_e=\operatorname{diag}(-1,-1,1)\) and \(W_e=\operatorname{diag}(1,i,-i)\), all others identity. With the stated clover defect \(\mu(P)=|\det((I-P)/2)|^2\), spatial masses vanish, while the mixed temporal defect is \(1/4\), giving endpoint masses \(m_x=m_y=1/24\). The two roots coincide on each selected pair, so \(K(U,U)=K(W,W)=L(0)\) and \(K(W,U)=K(U,W)=L(m)\). The seam increment is \(\delta B=\begin{psmallmatrix}0&D/2\\D/2&0\end{psmallmatrix}\). For \(v=(e_{\mathrm{first},x}+e_{\mathrm{last},x})/\sqrt2\), the two-point quadratic form on \((Jv,-Jv)\) is

$$
-2v^\dagger\delta Bv=-1/24<0.
$$

Therefore this **complete charged half-kernel** cannot have a universal pointwise positive-Gram factorization. The witness is not a replacement physical source and does not settle a projection onto the original prepared-source subspace. A nonpositive half-kernel does not make \(A^\dagger A\) negative. Fermionic reflection positivity has graded/twisted hypotheses; compare [Jaffe–Janssens](https://arxiv.org/abs/1607.07126).

## 5. Source-correct, singular-safe determinant methods

For left/right sea columns \(L,R\), original source vectors \(x,y\), and detector \(h\), put \(M=L^\dagger R\). With a common wedge ordering, the charged overlap and normal-ordered detector insertion are

$$
\langle L\wedge x,R\wedge y\rangle
=\det\begin{pmatrix}M&L^\dagger y\\x^\dagger R&x^\dagger y\end{pmatrix},
$$

$$
\langle L\wedge x,a^\dagger(h)a(h)(R\wedge y)\rangle
=-\det\begin{pmatrix}
M&L^\dagger y&L^\dagger h\\
x^\dagger R&x^\dagger y&x^\dagger h\\
h^\dagger R&h^\dagger y&0
\end{pmatrix}.
$$

The zero detector corner and the displayed sign encode the normal ordering. Projecting \(h\) onto an intermediate occupied sea is not this CAR identity. The determinants remain polynomial when \(\det M=0\); replacing them by \(\det M\) times an inverse formula, or a pseudoinverse, can destroy nonzero bordered amplitudes.

The retained determinant jets use column multilinearity. If \(M_i\) denotes column \(i\), then the mixed derivative is

$$
\partial_a\partial_b\det M
=\sum_i\det M[i\leftarrow\partial_{ab}M_i]
+\sum_{i\ne j}\det M[i\leftarrow\partial_aM_i,\ j\leftarrow\partial_bM_j].
$$

These are derivatives, not factorial-scaled Taylor coefficients. Safe reference blocks may be Schur-eliminated; singular or unsafe blocks stay in the augmented determinant. This is a coordinate reduction, not a discarded sea direction. Complex phases, both endpoint scalars, source Gram normalization and one common logarithmic scale remain owned by the actual consumer. Floating condition estimates are not certified native/integration error bounds.

These formulas transcribe the retained all-mode bordered-minor and physical-endpoint CAR-jet implementations; no new code or dependent data factories are published here. The underlying exterior/Wick algebra is established, including zero-overlap nonorthogonal matrix elements: [Burton (2021)](https://arxiv.org/abs/2101.10944).

## Scope of the archive

The useful content is the identities, hypotheses, proofs and counterexamples, not process receipts. Private memory, unsanitized archives and operational records are deliberately not reproduced. Historical numerical reports above remain reports with stated assumptions, not newly exercised certificates. No novelty priority has been established.

No retained argument supplies a common interacting ordered state with propagating charged excitations, physical Noether spin, photon dispersion/residue and controlled scattering on shared parameters. No measured electron mass, charge, fine-structure constant, material ranking or nuclear rate is predicted. Those scientific obligations remain unmet despite administrative issue closure and archival PR integration. **The campaign remains stopped.**
