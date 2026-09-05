# One fixed physical tag matches three carrier jets

Parent proposal validation (263 claims / 12 proposals) preceded opening
0135/0137. This is an additive construction using the ordinary-column
surface branch of 0135 at one fixed nonzero carrier q0. It does not
retune the tag separately for each macro sideband.

## 1. Exact four-row independence

Use the pressure normalization directly from the Bessel equation:

```
P(q,r)=A(q) r^m F_m(t(q) r²),
F_m(z)=sum_(n>=0) c_n z^n,
c_n=(-1)^n/[4^n n! (m+1)_n],
t(q)=(lambda(q)/a)^2.
```

At the fixed small q0>0, `A(q0)!=0`, `t(q0)>0`, and
`t'(q0)!=0`. For the ordinary surface branch
`t(q)=3q²+O(q⁴)` in the stated units, so its last condition holds for
sufficiently small nonzero q0. The series coefficient is
`-t r²/[4(m+1)]`; the factor `1/[2(m+1)]` belongs to the logarithmic
Bessel derivative, not the normalized pressure itself. Only nonzero
coefficients and rank are used here.

The four analytic radial functions

```
r^m, P(q0,r), partial_q P(q0,r), partial_q² P(q0,r)
```

are linearly independent on every nonempty open radial interval.
Indeed the q derivatives are an invertible triangular change from
`P,partial_t P,partial_t² P`, with diagonal entries
`A,A t',A(t')²`. After division by r^m, the n>=1 coefficients of a
putative relation give a quadratic polynomial in n which vanishes at
every positive integer. Its coefficients vanish, and the n=0 row then
fixes the coefficient of r^m to zero.

An exposing finite determinant proves the same statement. Using the
first four powers `1,r²,r⁴,r⁶`, the coefficient determinant for the
normalized t-derivative rows is

```
2 c_1 c_2 c_3 t0³
=t0³/[24576 (m+1)³(m+2)²(m+3)] !=0.
```

The actual q-derivative determinant is multiplied by
`A(q0)³ t'(q0)³`. Near zero this behaves as q0^9 for the normalized
branch; the construction is deliberately at a fixed nonzero q0.
No uniform endpoint conditioning or zero-carrier rotor is asserted.

## 2. Four actual smooth radial controls

On an interval where chi is positive, independence gives four points
whose evaluation matrix is nonsingular. Otherwise all evaluation
vectors would span a proper subspace and contradict the analytic
independence. Four sufficiently narrow disjoint smooth radial bumps
retain that determinant by continuity. Integrating with the common
positive radial factor gives four independent moment rows

```
Q=integral chi b r^(m+1) dr,
T^(j)(q0)=integral chi b r partial_q^j P(q0,r) dr, j=0,1,2.
```

The signed b describes radial lobe orientation; the mass fraction
`mu chi chi_z[1+epsilon b cos(m theta)]` stays nonnegative by scaling
all four bump coefficients together. Fix mu, for example mu<=1/2,
first, then scale until `|epsilon b|<1`. All target ratios are
unchanged by this final scaling, and Q remains nonzero.

## 3. Match the physical spin/action jet with the same tag

For the ordinary-column conventions of 0135 write

```
R(q)=-beta(q) Z0 sigma(q)²[2 Omega+sigma(q)]
      /[rho mu m pi Zc(q)² G(q)].
```

These are actual pressure, KKS and transported-tag integrals. The
physical equality is `T(q)/Q=R(q)`. Choose the four moment controls to
prescribe any nonzero Q and

```
T(q0)=Q R(q0),
T'(q0)=Q R'(q0),
T''(q0)=Q R''(q0).
```

Then the *same frozen tag* has

```
j_tag(q)-M_can(q)=O((q-q0)³).
```

The proportionality factor converting `T-Q R` to this difference is
smooth and finite because q0 and the actual angle chart are fixed
before the macro expansion. The tag itself is not differentiated or
changed with q. This construction preserves the three mechanical/
canonical jets needed by a macro wavevector expansion about q0.

## 4. Interface to the generalized force-free construction

If the convention instead uses `T(q)=integral chi b r P(q,r)/s(q)`,
where s(q) is radially constant and nonzero, the pressure and its first
two derivatives undergo another invertible triangular row change.
Thus the same rank proof applies exactly to the ordinary-column
`P/s` convention.

For the generalized force-free case in 0140,
`s(r,q)=nu(q)-m Omega-q h(r)` is generally r-dependent. Division by
that function is **not** the preceding scalar row change. The proper
transfer is continuity of the particular fixed-q0 four-bump determinant
under the controlled pressure/Doppler perturbation, including two q
derivatives. If the new 4x4 matrix is M0+Delta and
`||M0^-1 Delta||<1`, it remains invertible. The smoothing/force-free
construction supplies the actual row bound; no new rank hypothesis or
zero-carrier limit is substituted for it.

`radial_rank_verify.py` checks the nonzero coefficient determinant,
the Bessel recurrence, the complete derivative row transformations,
and the mechanical-jet cancellation. `radial-rank-first-run.txt` is
its first execution receipt. This artifact is frozen for reuse by
0140; the parent spatial closure is not asserted here.
