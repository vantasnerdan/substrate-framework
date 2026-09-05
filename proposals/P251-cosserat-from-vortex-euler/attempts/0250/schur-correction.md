# Same-frequency Schur coefficient (bounded correction)

The parity argument in `axisymmetric-hybrid.md` removes the even cross rows,
but it does not make the odd rows a small perturbation of the acoustic block.
Writing the retained normalized gain in acoustic/optical order gives

\[
 M(K)=\begin{pmatrix}K^2 a_2+o(K^2)&K b_1+o(K)\\
 K c_1+o(K)&B_0+O(K^2)\end{pmatrix},
\]

where `B_0` is the nonzero three-row optical determinant from the reflected
law, `b_1` is the optical-to-hybrid coefficient from the retained odd part of
(24), and `c_1` is the acoustic-to-(theta,G,S) coefficient from the retained
odd curl/current rows in (24)--(25).  Since `B_0` is invertible, the actual
acoustic Schur coefficient is

\[
 a_{\rm Sch}=a_2-b_1B_0^{-1}c_1,
 \qquad
 S_A(K)=K^2a_{\rm Sch}+o(K^2).
\]

Thus `O(K)` on both odd sides contributes at the same order as the raw
acoustic gain; a Neumann-smallness claim is invalid unless a separate norm
inequality is supplied.  The current text records only the orders and the
formal sources of `b_1,c_1`; it does not retain their normalized numerical
rows (nor a proof that `c_1=0` on the prepared acoustic column).  Therefore a
nonzero lower bound for `|a Sch|` cannot be derived from the present source
receipts.  The established local hybrid gain, optical determinant, and parity
cancellations are unaffected.

Minimum repair: evaluate the normalized first odd jets of (24)--(25) on the
actual prepared columns, form the displayed scalar/matrix product, and record
an exposed bound `|a_2-b_1 B_0^{-1}c_1| >= c_* > 0` on the common frequency
interval (or construct the acoustic column with `c_1=0` and prove that row).
Only that bound licenses the claimed finite-`K` joint inverse.

## Failure-derived continuation: retained-row evaluation

For the normalized hybrid row `H`, write the first odd jet of the optical
source as

`b_1 = partial_K H_opt|_{K=0}`.

The retained formulas identify this jet as the sum of the optical mechanical
spin row (21) and the optical symmetric first-moment/centroid row from
(24)--(25), after the same phase, tag, and reflected-law normalization used
to define `B_0`.  For the prepared G-cancelled acoustic column `A`, write

`c_1 = partial_K (theta,G,S)_A|_{K=0}`;

its entries are respectively the first-K angle row and the first-K current/
spin rows of (24)--(25), including the Bloch phase and material-centroid
term.  Hence the coefficient to expose is exactly

`a_Sch = a_2 - b_1 B_0^{-1} c_1`.

The present 0250 receipts give only the orders and parity of these jets, not
their normalized integrals (the displayed spin formula (21) is local and no
full reflected integral is tabulated).  Therefore an evaluated nonzero
number or lower bound cannot honestly be supplied from the retained source
formulas alone.  This is a concrete missing-data branch, not a claim that
the coefficient cancels.  The next source calculation must integrate (21),
the symmetric first moment, and all three acoustic odd rows with the same
`theta,phi` and tag normalization, then report `a_Sch` on the common band.
If it vanishes, the prescribed repair is a second axisymmetric
polarization/profile or an additional fraction/current-null combination,
followed by the same integral.

## Bloch-phase expansion (the actual available integral)

For any prepared column, let `F_H(x)` denote the Cartesian integrand in the
full material variation (8), after inserting `D_t xi=Du xi+v`, including the
centroid subtraction, and let `F_O(x)` denote the three optical rows in
(24)--(25).  The physical Bloch source is `exp(i K.x) xi_K`; therefore the
first jets that enter the Schur product are not derivatives of the
`K=0` coefficient alone but

\[
 b_1=\Pi_H\!\left[i\!\int x\,F_H^{\rm opt}(x)\,dx
                  +\int \partial_KF_H^{\rm opt}(x)|_0\,dx\right],
\]
\[
 c_1=\Pi_O\!\left[i\!\int x\,F_O^{\rm ac}(x)\,dx
                  +\int \partial_KF_O^{\rm ac}(x)|_0\,dx\right].
\]

Here `Pi_H,Pi_O` include the common theta/phi measure, tag dual, and the
same reflected-law normalization as `B_0`; `partial_K F` contains the
first-K change of the prepared cohomological lift and of the Leray/pressure
response.  Substitution of (8), (21), and (24)--(25) gives these integrals
directly and shows why the material centroid term belongs in `b_1`.

For the circular ring, the already retained integrands (11)--(12) have only
the displayed `m=+-1` poloidal harmonic.  Multiplication by `x` produces
the `m=0, +-2` harmonics; hence the first-K coefficient is determined by
the corresponding radial envelopes and by the (unspecified) finite-carrier
pressure terms.  No cancellation follows from theta/phi parity alone for an
oblique `K`, because the transverse hybrid components are precisely the
`K_perp K_parallel` projections retained in (15).  The coefficient is thus
an explicit finite integral of the two displayed source profiles, but its
value and a nonzero lower bound require the actual `A(s),V(s),W(s),chi(s)`
and pressure correction functions, which are not present in the 0250
receipts.  The route remains open at this concrete integral evaluation;
it is not licensed to assert either cancellation or invertibility.

## Leading Beltrami substitution

Using `V=lambda A J1(lambda s)` and `W=lambda A J0(lambda s)` from 0211,
with `nu=mV/s` and `C=Ns/(lambda*m*V)`, the retained leading tensor becomes

\[
\delta B_{zzz}^{(u)}=i\pi^2N s^2V/(2\lambda)+O(1),\qquad
\delta B_{zxx}^{(u)}=\delta B_{zyy}^{(u)}
 =i\pi^2NV\left(R^2/\lambda+3s^2/(4\lambda)\right)+O(1).
\]

The exact G-cancelled p column is only `O(1)`, so these `O(N)` terms survive.
After localizing the selectable envelope and tag near any fixed `s_*>0`,
the raw quadratic coefficient is therefore nonzero, with
`h_perp/h_parallel=2R^2/s_*^2+3/2+O(s_*/R)`.

The first-K rows still contain the selectable radial envelope, tag
derivative, and finite-pressure corrections in addition to these angular
integrals.  Hence `a_2` is now explicit and nonzero, but the sign or
cancellation of `b_1 B_0^{-1}c_1` remains a radial-profile choice.  If one
profile cancels, a second disjoint envelope in the same fixed annulus gives
the registered alternative column while preserving the positive leading
tensor.
