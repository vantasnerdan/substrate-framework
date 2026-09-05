# An exact positive Euler material-stiffness witness

## Primary evidence and what transfers

East, Zrake, Yuan and Blandford's [primary paper, arXiv:1503.04793v3](https://arxiv.org/pdf/1503.04793v3)
gives a finite Fourier trial in Appendix B.2, Table II (PDF page 9), labeled
incompressible. Its background is the two-wave case of equation (10).
Table I and Tables III–IV instead label their displayed trials compressible;
they are not used here. Equation (6) defines the magnetic potential functional.
The article's dynamical growth rates and kinetic denominator are not Euler
imports. Only its trial field is used to construct the following independently
computed Euler quadratic functional. The periodic boundary choice avoids
importing a cylindrical wall contribution from a different source geometry.

## Fully specified rational field

On the torus `[0,2pi]^3`, use the real smooth field

    u=(-sin(2y),cos(2x),cos(2y)-sin(2x)),
    curl u=2u,    div u=0,    p0=-rho|u|²/2.

`verify.py` records all eight positive-lattice rows of Table II as integer
real and imaginary numerators divided by 10000. It includes their complex
conjugates at negative wavevectors. If `c_k` denotes that printed rational
coefficient, the actual trial used here is

    Xi(x)=sum_k P_k c_k exp(i k.x),
    P_k=I-k k^T/|k|²,
    k=(1,+/-1,0),(1,+/-3,0),(3,+/-1,0),(3,+/-3,0)
      and their negatives.

This exactly defines a mean-zero smooth real divergence-free displacement.
The projection is necessary: rounded rows at `(1,-3,0)` and `(3,-1,0)` have
small nonzero printed divergence. Nothing is inferred from a tolerance;
the projected coefficients and all arithmetic are rational complex numbers.
The finite witness is a newly specified exact trial derived from the table,
not the assertion that the rounded table itself is exactly solenoidal.

## Exact evaluation and independent form

Let brackets denote the volume average, `F=Xi cross u`, and `w=curl F`.
The exact convolution yields

    <F.w>      = 153402581/500000000,
    <|w|²>     = 190501151/312500000,
    <|Xi|²>    = 14496029/125000000,
    K/(rho V)  = 2<F.w>-<|w|²>
               = 5008301/1250000000 > 0.

Here `V=(2pi)^3`. Thus the material stiffness-to-metric ratio is exactly
`5008301/144960290`. The denominator is the material metric `rho integral
|Xi|²`, not the source's magnetic-field-weighted dynamical denominator.

The verifier independently forms `p0=-rho|u|²/2`, differentiates its Fourier
coefficients twice, and evaluates

    <Xi.Hess(p0).Xi/rho>-<|(u.grad)Xi|²>
       =5008301/1250000000.

This second evaluation uses the actual pressure and advective terms of the
Euler Jacobi action, rather than assuming the source's potential sign.
The 0084 exact identity relates the two expressions. Reversing the helicity
sign or omitting it produces a negative value and is explicitly exposed.
In the source convention the magnetic potential per volume is precisely
`-5008301/2500000000`, confirming the opposite functional sign without
importing any MHD evolution equation.

## Verdict and consumer boundary

`route_verdict: established` — an exact smooth periodic divergence-free
material displacement has strictly positive Euler Jacobi stiffness on the
specified stationary Beltrami field. `evidence_scope: symbolic_verified`
for this fixed finite Fourier functional, with no discretization error or
unresolved small eigenvalue. The source numerical search selected a trial;
the sign here is established by independent exact evaluation, not that search.

This result supplies a concrete positive margin for main's 0090 compact
localization and same-EPS-field transfer. It is not a proof of an invariant
Euler mode, fixed-Kelvin reduced inertia, Euler spectral stability, a core
rotation jet, or the parent Cosserat closure. Those physical constructions
retain their own action and material-observable conditions.

## Frozen provenance

`verify.py` reuses 0040's exact Fourier convolution, curl and projection
implementation and the canonical `CheckLedger`; it creates no duplicate
Fourier engine. Its first execution passed all eight checks, exit zero,
4.512290452 seconds. `stdout.txt` preserves that execution. Ruff and scoped
diff checks pass. No source-pressure sign repair or solver rerun was needed.

SHA-256:

```
cb9224b5951534f8170eaa17a8015b13c5c3cb15aee76250f529cfee7d81f9c4  1503.04793v3.pdf
8fab3eea98d2d96bc42c43650baa840ce1e537413198c3358ee560c73a8ddf40  ../0040/fourier_orbit.py
73a0e0a54c0aeb565bb334f97122143ed3e837288d8f11e25eb7e952378057ad  verify.py
```
