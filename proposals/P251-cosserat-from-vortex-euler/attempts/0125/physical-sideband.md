# Actual Euler polarizations, a visible complement, and its repair interface

This is a replacement-route derivation, not a review or alteration of
accepted C-CST-008/009/010. The uniform-rotation calculation is exact on its
stated domain. Its transfer to compact EPS cells requires the explicitly
retained localization and moment achievements below.

## 1. Zero frequency variance does not remove physical polarization memory

Take `u=Omega e_z cross x`, `Omega>0`, on the whole uniform-rotation core.
An axial drift can be restored by the common convective phase. In rotating
coordinates, a Kelvin wave has fixed wavevector `k=(epsilon,0,N)`,
`N>0`, and solves the full pressure-eliminated Euler equation

```
b_dot = -2 P_k A b,   A v = Omega e_z cross v,   k.b=0.
```

Put `c=N/sqrt(N²+epsilon²)`, `s=epsilon/sqrt(N²+epsilon²)`,
`e_a=(c,0,-s)`, `e_b=(0,1,0)`. Then
`b=e_a cos(2 Omega c t)-e_b sin(2 Omega c t)`.
The actual core-vorticity angle is `Phi=e_z cross curl(v)/(2 Omega)`.
Its normalized horizontal components in rotating coordinates are
`(cos(2 Omega c t),-c sin(2 Omega c t))`, not a unit circular
polarization at nonzero epsilon. Returning to the laboratory gives

```
F(t) = (1+c)/2 exp[i(Omega-2 Omega c)t]
     + (1-c)/2 exp[i(Omega+2 Omega c)t].
```

Pairing actual opposite rotations takes the appropriate real component:

```
q(t) = a_- cos(omega_- t) + a_+ cos(omega_+ t),
a_-= (1+c)/2, a_+=(1-c)/2, omega_±=Omega(1±2c).
```

Consequently `q''(0)=-Omega²` exactly, whereas
`q''''(0)=Omega⁴[1+16 c²(1-c²)]`. The second spatial jet is

```
q(t)=cos(Omega t)
 + epsilon²/N² {Omega t sin(Omega t)
               +[cos(3 Omega t)-cos(Omega t)]/4}
 + O_T(epsilon⁴).
```

This is new contrary evidence to the *route* “axial zero first-order
frequency splitting alone supplies an autonomous physical-angle pencil.”
The zero-order paired rotor remains intact. The strongest exact closure
for this prepared scalar is the two-branch equation

```
(D_t²+omega_-²)(D_t²+omega_+²)q=0.
```

Its branch coordinates are physical differential observations:
`q_+=(q''+omega_-² q)/(omega_-²-omega_+²)`, `q_-=q-q_+`.
Thus the branch state is not an arbitrary input oscillator. At epsilon=0
the fast residue is zero; its order-epsilon² residue matters for the
second spatial jet. There is no contradiction with a conditional action
whose physical observation is separately specified.

## 2. A genuine Euler complement has the missing core observation

Let `zeta=x+i y`, `e_±=(1,±i,0)` and use the axial/time factor
`exp(i N z-3 i Omega t)`. The following polynomial profiles solve the
complete linearized Euler equation, including its pressure:

```
v_perp = a zeta² e_- + (2a/N²)e_+,
v_z = 4 i a zeta/N,
pi = 8 i Omega a zeta/N².
```

Here pi is pressure divided by density. Direct differentiation proves
`div v=0` and `v_t+u.grad v+A v+grad pi=0`. The actual Lin displacement
`xi=curl(v)/(2 i Omega N)` satisfies
`xi_t+u.grad xi-A xi=v`; thus this is a Kelvin-tangent mode, not merely
a forced solution. Its physical core angle is

```
Phi_3 = i a/(Omega N) e_+ exp(-3 i Omega t) != 0.
```

Adding its opposite paired phase with coefficient
`-epsilon²/(4N²)` cancels the observed fast line in section 1. Adding
`+epsilon²/(4N²)` of the original fundamental restores the prepared
initial angle. The corrected local trace is

```
cos(Omega t)+epsilon² Omega t sin(Omega t)/N²+O_T(epsilon⁴),
```

corresponding to squared frequency `Omega²-2 Omega² epsilon²/N²`.
This is an exact *local polynomial* preparation result, not a positive
spatial stiffness result. The polynomial grows transversely and does not
itself have finite action on the whole transverse plane. A fixed finite
cutoff makes its order-epsilon² amplitude have order-epsilon⁴ quadratic
action, but its localization error still multiplies epsilon². A small
nonzero collar error cannot be relabeled an exact second-jet coefficient.
The ongoing 0128 construction owns the controlled localization upgrade.

## 3. Mechanical spin is an independent row

For a small isotropic material tag, write
`integral rho r_i r_j = I_tag delta_ij`, and `D=grad xi` at its centre.
Keeping both displacement and velocity variations gives, at leading
tag-size order,

```
delta S = I_tag [curl v - Omega (D+D^T)e_z],   div xi=0.
```

Indeed it is the integral of
`xi cross A r + r cross (A xi+v)`. For the polynomial mode above,

```
S_3 = a I_tag/N e_+ exp(-3 i Omega t)
    = (I_tag/3) Phi_3_dot.
```

The first equality and the complete Euler/Lin equations are checked in
`verify.py`; the spin-rate identity is asserted, not inferred from a
canonical momentum. This is a small-tag coefficient, not the exact
finite-tag inertia of 0123. It proves why angle-only cancellation is
insufficient.

The 0128 independent finite-tag calculation supplies actual modes with
the same 3 Omega frequency and distinct spin-rate rows `J_1 != J_2`.
For target angle coefficient a and spin-rate coefficient b, the exact
control amplitudes are

```
c_1=(b-J_2 a)/(J_1-J_2),
c_2=(J_1 a-b)/(J_1-J_2).
```

The same linear system applies separately to both phases, or directly
to complex amplitudes. These coefficients match measured rows without
changing the angle definition. Their existence/localization and full
action cross terms remain dependencies of that upgrade, not assumptions
smuggled into this derivation. Matching these two rows removes this
specific fast observed residue. It does not prove invariance of the
whole periodic Euler complement or remove the stress/shape rows in the
second spatial response.

## 4. Compact axial profiles require an infrared repair

For the actual coherent compact profile h, a single carrier N cannot be
substituted into the spatial calculation. Use the convention
`hhat(k)=integral exp(-ikz)h(z) dz`. If h is smooth compact and odd,

```
hhat(k)=-i B_total k+O(k³),   B_total=integral z h(z) dz.
```

The actual rotating-Euler projector and propagator have second transverse
derivatives of order `k^-2` at an axial wavevector. For example,
`P_(epsilon,0,k)e_x` has horizontal entry
`k²/(k²+epsilon²)` and mixed entry `-epsilon k/(k²+epsilon²)`.
The horizontal second-jet coefficient applied to h is `-hhat(k)/k²`.
When `B_total!=0` this is not square integrable near zero. Hence the
full prepared field already lacks this second jet at initial time for
this preparation. A vorticity/angle row can hide the divergence by its
additional k factor. This distinguishes an observed jet from a finite-
energy state jet and does not refute alternative preparations.

A constructive shape extension removes it. Choose an odd smooth f
supported in two intervals away from zero, with
`B_f=integral s f(s) ds !=0`. For L large enough that its support is
outside the old profile and every tag, set

```
r_L(z)=-B_total/(L² B_f) f(z/L),   h_new=h+r_L.
```

Then the global first moment is exactly zero; h and all tag moments
are unchanged on their original support. Disjoint supports give the
explicit action-norm increment

```
integral r_L² dz = B_total²/(L³ B_f²) integral f² ds.
```

Now `hhat_new=O(k³)` and `hhat_new/k²` is square integrable. For the
actual bounded rotating-Euler matrix propagator and Leray-prepared data,
its derivatives at nonzero k obey `|partial_epsilon^j| <= C_T |k|^-j`.
Splitting the k integral into `|k|<=|epsilon|` and its complement proves
an L2 second-jet remainder `O_T(|epsilon|^(7/2))` for even-symbol entries.
The mixed odd entries have an ordinary third-order term in general;
their second-jet remainder is `O_T(|epsilon|³)`. Thus the full state has
a genuine `o(epsilon²)` finite-time expansion. Smooth compact h controls
the high-frequency tail. No numerical floor or fitted threshold enters.

This extension costs a finite new geometric scale L and a known positive
action increment. Restoring an already exact action/spin match requires
an explicit profile retuning; unchanged tagged spin does not imply
unchanged canonical action. The parent must retain this cost and return
support in the declared geometry. The argument is a positive regularity
repair for uniform-rotation propagation, not an EPS-cell spectral theorem.

## Evidence and route decisions

`verify.py` checks the exact physical moments, fourth-order closure,
polynomial Euler/Lin mode, measured spin, pressure projector derivatives,
two-row control interface, return moment/action scalings and ordered jet
coupling. `first-run.txt` preserves a representation-only failed equality
check: two equal SymPy expressions had different factorizations. The
repair tests their expanded difference; `repaired-run.txt` and
`final-run.txt` preserve the subsequent receipts.

Candidate A is established as a physical two-branch closure and refutes
the scalar zero-splitting shortcut by its explicit mechanism. Candidate B
establishes the local visible polynomial repair; its finite-energy,
finite-tag same-EPS upgrade is active in 0128. The exterior return
establishes a full-state infrared regularity repair on the comparison
domain. Candidate C's exact periodic response/memory is recorded in
`periodic-second-jet.md`. None of these child verdicts asserts parent
completion or unrestricted invariant rotor dynamics.
