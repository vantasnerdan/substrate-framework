# Attempt 0009 — N2 block B3 (bend channel), outer energy via the K-identity

## Established

1. Exterior gradient-energy identity (modified-Bessel ODE + integration by
   parts): int_{x0}^oo x [K1'^2 + (1 + 1/x^2) K1^2] dx = -x0 K1(x0) K1'(x0).
   Verified numerically to machine precision at x0 = 0.3, 0.7, 1.5 and by
   symbolic small-x0 expansion:

     -x0 K1(x0) K1'(x0) = x0^{-2} - 1/2 + O(x0^2 log x0).

2. Outer m=1 wavenumber-k bending-channel energy (displacement amplitude
   xi, wavenumber k, x0 = k a):

     E_out/L = (rho Gam^2 xi^2 / 4 pi) · [ (k a)^{-2} - 1/2 + O((k a)^2 log) ].

## Matching structure (declared, per the analytic-receipt discipline)

- The (k a)^{-2} piece is the core-translation sector: the swirl field is
  translation-invariant to this order, so it cancels against the inner-core
  re-computation at the displaced position (inner swirl energy from block
  B1: rho Gam^2/(16 pi), translation-invariant). The cancellation is the
  standard matched-asymptotics bookkeeping and is checkable exactly.
- The -1/2 piece supplies the (ka)^0 constant of the bend modulus:

    B(k) = (rho Gam^2 / 4 pi) · [ ln(2/(k a)) - gamma + c_inner ],

  with gamma the Euler constant (from K1 ~ 1/x + (x/2)ln(x)-structure via
  K0 ~ -ln(x/2) - gamma) and c_inner the exact core-interior constant
  (top-hat mollifier), owed to the matching assembly (attempt 0010).

- The Krein finding of attempt 0008 applies ONLY to the m=2 polarization
  channel; the m=1 bending channel is a standard positive-energy wave
  (translation-sector cancellation + log), which is why B(k) is extractable
  without the canonical-receipt machinery.

## Status

K-identity + outer energy: established (machine-verified identity + exact
small-x0 expansion). Bend modulus assembly (c_inner via the core-interior
integral, Euler-constant bookkeeping, and the mutation battery): owed to
attempt 0010, which completes C-CST-002's verifier (B1 two-route, B2 m=2
energies + dispersion probes, B3 bend modulus with mutations).
