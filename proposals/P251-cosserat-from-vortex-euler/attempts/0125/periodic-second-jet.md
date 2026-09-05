# Full periodic Euler/Lin second jet, physical current and retained memory

This derives the finite-time response directly from 0116's actual
periodic Euler operator. It supplies no autonomous constitutive closure
by assumption. The underlying smooth same-EPS periodic background and
actual transported tags are those of 0116/0117. All expansions below are
coefficients, so the second coefficient is one half the second derivative.

## 1. The actual operator and its directional pressure derivatives

Fix a nonzero direction n and write the Bloch vector `kappa=epsilon n`.
The cell-periodic velocity w and Lin displacement eta satisfy

```
w_dot = A_epsilon w,
A_epsilon=-P_epsilon[C_0+i epsilon (u.n)],
C_0 w=u.grad w+(grad u)w,
eta_dot=-u.D_epsilon eta+(grad u)eta+w,
D_epsilon=grad+i epsilon n.
```

The physical perturbed material velocity is `w+(grad u)eta`, not w
alone. On a nonzero cell Fourier mode ell, set `L=|ell|²`, `a=ell.n`.
The pressure projector coefficients are

```
P_0=I-ell ell^T/L,
P_1=-(n ell^T+ell n^T)/L+2a ell ell^T/L²,
P_2=-n n^T/L+2a(n ell^T+ell n^T)/L²
    +( |n|²/L²-4a²/L³ )ell ell^T.
```

On the zero cell mode, use `P_n=I-n n^T/|n|²` for all epsilon,
including the directional epsilon=0 limit; both derivatives there are
zero. Replacing it at the endpoint by the unrestricted mean projector
would introduce a false discontinuity and admit an unphysical longitudinal
mean for the nonzero Bloch fiber. The retained mean displacement and
velocity are transverse to n; pressure is still present at every order.

It follows, with multiplication by `u.n` written in its actual order,

```
A_0=-P_0 C_0,
A_1=-P_1 C_0-i P_0(u.n),
A_2=-P_2 C_0-i P_1(u.n).
```

For `z=(eta,w)`, let `L_j` be the corresponding block coefficients.
The displacement diagonal has first coefficient `-i(u.n)` and zero
second coefficient. Its upper-right block is I at order zero and zero
at higher orders. This includes the Lin/velocity relation throughout.

Choose initial `z_epsilon(0)=E_epsilon a` from actual Bloch Kelvin
generators, together with any explicitly declared harmonic circulation
columns needed for mean translation. No arbitrary oscillator is inserted.
The Leray derivatives of this preparation are part of E_1,E_2; dropping
them would change the prepared physical state.

## 2. A genuinely closed finite-time spatial jet in the Euler state space

Writing `z=z_0+epsilon z_1+epsilon² z_2+O(epsilon³)` gives the exact
triangular equations

```
z_0_dot=L_0 z_0,
z_1_dot=L_0 z_1+L_1 z_0,
z_2_dot=L_0 z_2+L_1 z_1+L_2 z_0,
z_j(0)=E_j a.
```

These close without an infinite hierarchy of *spatial derivatives*, but
each retained z_j is still an actual cell field. They are not a finite
constitutive pencil. With `U_0(t,s)` the actual unperturbed Euler/Lin
propagator, their explicit ordered solution is

```
z_0(t)=U_0(t,0)E_0 a,
z_1(t)=U_0(t,0)E_1 a+integral_0^t U_0(t,s)L_1(s)z_0(s) ds,
z_2(t)=U_0(t,0)E_2 a
       +integral_0^t U_0(t,s)[L_1(s)z_1(s)+L_2(s)z_0(s)] ds.
```

The double ordered `L_1 U_0 L_1` term is essential. Identical cells
remove an ensemble distribution of backgrounds but do not remove this
within-cell propagation into and back from the complement.

On a fixed torus with smooth u, the nonzero Fourier-mode gap makes the
projector expansion analytic for `|epsilon n|` smaller than that gap.
`P_j` has order -j on the nonzero modes. Thus A_1 has order zero and
A_2 has order -1; higher coefficients and the finite-mode mean terms are
bounded perturbations at the regularity of a smooth Euler solution.
The usual finite-time transport energy estimate and Duhamel formula
give the displayed `O_T(epsilon³)` remainder for smooth prepared data,
with enough Sobolev regularity for the requested differentiated tag/core
observations. This is a local-in-time response statement on the given
smooth Euler existence interval. It neither extends that interval nor
assumes spectral stability. Unlike the continuous axial comparison, a
fixed periodic cell has no nonzero-mode infrared accumulation at zero.

## 3. Actual mean stress and the second hybrid-current moment

Let `m_j=<w_j>`. The exact cell-average Euler equation yields

```
m_0_dot=0,
m_1_dot=-i P_n <u tensor w_0+w_0 tensor u> n,
m_2_dot=-i P_n <u tensor w_1+w_1 tensor u> n.
```

The same prepared cell response therefore supplies translation, not an
independently chosen spring force. Mean velocity is not yet the hybrid
centroid current. For each tag, let r denote displacement from its
actual centroid and define its central momentum moments

```
C_ij=integral_tag rho v_i r_j,
T_ijl=integral_tag rho v_i r_j r_l.
```

In the response below C and T mean their full first variations, summed
over tags per unit volume with centroid Fourier phases; tag displacement,
centroid displacement, material velocity and phase variations are all
included. The untagged ambient current is retained unchanged in the
hybrid definition. Taylor expansion of the actual current gives

```
J_E-J_H=-i kappa_j C_ij-(1/2)kappa_j kappa_l T_ijl+O(|kappa|³).
```

Consequently the second-jet observation is

```
J_H,0=rho m_0,
J_H,1=rho m_1+i C_0 n,
J_H,2=rho m_2+i C_1 n+(1/2) T_0:(n tensor n).
```

The angular-momentum decomposition is exactly
`C_ij=I_dot,ij/2-epsilon_ijl S_l/2`, reproducing 0117 at first order.
The second central momentum T is an additional physical row at the
second spatial order. Isotropy can remove specified tensor components
only after their full response is included; first-order isotropic
spin/shape separation does not by itself set T_0 to zero. Integrating
`J_H/rho` supplies the hybrid displacement with its actual prepared
initial value. Dropping T would be an unsupported second-order join.

The angular observation is the actual transported-tag/core row. For
isotropic transverse preparation at order zero it uses
`(E P_n)^-1 E q_n`, not the raw mean. Its higher derivatives, and those
of spin and C,T, act on z_0,z_1,z_2 with the usual ordered coefficient
convolution. In particular, no averaged scalar frequency supplies these
rows.

## 4. Exact finite physical observations carry a complement memory

Here is a formulation that keeps the retained coordinates physical.
Let `O(t)` be a finite full-rank list of actual independent observation
rows (for example hybrid displacement/current and angle/rate; append
independent shape rows when desired). Spin rows already equal to a rate
row on a proposed prepared family are outputs, not duplicate coordinates.
Choose a smooth right inverse E with `O E=I`, and put

```
x=O z,  y=z-E x,  O y=0,
B=O_dot+O L,
A=B E,
F=L E-E_dot-E A,
G=L-E B.
```

Direct differentiation of the actual Euler/Lin equation gives

```
x_dot=A x+B y,
y_dot=F x+G y.
```

The evolution generated by G preserves the moving null space of O:
`O_dot y+O G y=0` whenever `O y=0`; also `O F=0`.
Let V be that null-space propagator. Then the exact physical equation is

```
x_dot(t)=A(t)x(t)+B(t)V(t,0)y(0)
         +integral_0^t B(t)V(t,s)F(s)x(s) ds.
```

This equation is derived from the actual Euler operator; it does not
postulate a memory kernel to fit the response. Changing the lifting
changes its representation, not the physical transfer. The source
`y(0)` is fixed by the preparation and is not silently reset after the
initial time. A causal full-field response and a closed autonomous
finite moment dynamics are therefore different achievements.

For a proposed invariant base embedding, `F_0=0` is a concrete property
to prove, not a consequence of positive finite action. At the first
spatial order the complement source includes

```
F_1=L_0 E_1+L_1 E_0-E_1_dot-E_1 A_0-E_0 A_1.
```

In a time-independent invariant base this is the actual Bloch Sylvester
equation. With physical normalization `O_epsilon E_epsilon=I`, a
successful `F_1=0` preparation removes the first-order complement
excitation. Its second-order analogue, and the observation rows, decide
the second-order closure. A vanishing retained block or zero first-order
frequency splitting does **not** imply `F_1=0`. If F_1 is nonzero, its
complement can already be observed at first order when the base output
couples to that complement. Even when symmetry removes this first-order
response (or an invariant base projection makes that coupling zero),
the ordered complement return can produce order-epsilon² memory despite
a vanishing first-order retained block.

The polynomial 3 Omega preparation in `physical-sideband.md` is an
explicit local solution for one observed second-order residue of this
normal-form problem. The two-row 0128 interface can repair its spin
simultaneously. It does not solve the full periodic Sylvester problem:
the actual periodic complement can contain resonant, continuous-spectrum
or neutral responses, and the observation map also contains stress and
quadrupole rows. A fixed finite-time inverse is not automatically a
bounded time-independent spectral inverse.

## 5. A finite prepared-family equation, with its precise scope

There is also a useful nonautonomous finite-time construction. For d
independent prepared columns, compute from sections 1–3

```
T_epsilon(t)=O_epsilon(t) U_epsilon(t,0) E_epsilon,
x(t)=T_epsilon(t) a.
```

On any interval where the *physical* observation matrix T_0 is
invertible, the measured prepared state obeys the exact finite equation
`x_dot=G_epsilon(t)x`, `G_epsilon=T_epsilon_dot T_epsilon^-1`.
This is computable from the actual Euler response, with second jet

```
G_0=T_0_dot T_0^-1,
G_1=(T_1_dot-G_0 T_1)T_0^-1,
G_2=(T_2_dot-G_0 T_2-G_1 T_1)T_0^-1.
```

It is predictive on the declared prepared family and includes the
physical hybrid translation when the corresponding columns/rows are
included. It may need physical chart changes and still requires the
full cell propagator to compute its coefficients. It is not a proof of
a time-independent Cosserat pencil or a universal moment closure off
that prepared family. No invertibility receipt for a full P251
translation/angular chart is claimed in this attempt.

## Child result and next construction

The established result is the closed field-valued second spatial jet,
the complete second-order hybrid-current observation, and the exact
finite physical-observation memory equation. The measured two-branch
comparison and the visible Euler complement identify a concrete repair,
not a universal no-go. The next constructive achievement is to combine
the infrared-regular moment-matched profile and localized two-row
control with actual periodic E_1,E_2 and prove the required Sylvester/
observation identities, or retain the physically defined branches and
shape moments they excite. Identical fixed-lambda cells alone do not
prove those identities. The autonomous same-EPS parent remains active.

`jet_verify.py` checks the moving physical-null-space identities,
reconstruction of the full operator, the complete varying-centroid/
material quadrupole identity, and the noncommuting prepared-transfer
coefficients. `jet-first-run.txt` records all nine exact checks.
