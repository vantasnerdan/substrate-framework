# Derivative-loss repair

This receipt records the scientific repair requested during the author-stage
audit of P253/0075. It is not an independent review.

## Exposing mechanism

The discarded velocity route placed `T=u tensor u` in a weighted space with
`s` derivatives, applied the order-zero pressure Riesz transforms, and then
used `grad p` in an `s`-times differentiated velocity estimate. The displayed
hypotheses supply only `s-1` weighted derivatives of `grad p`. Extra global
Sobolev regularity does not supply the missing far-field weight, so that
estimate was invalid.

## Positive replacement

The final derivation fixes integer `s>=4`, transports
`q=curl(u-U_*)` in `Z_(3+gamma)^(s-1,alpha) intersection L^1`, and uses
`integral q=0` in the cancellation formula

    Bq(x)=integral [L(x-y)-L(x)]q(y)dy.

This proves `Bq in Z_(2+gamma)^(s,alpha)`. The exact vorticity equation has no
pressure, closes through order `s-1`, and the reconstruction restores the one
velocity derivative. The pressure multipole remains a separately proved
low-order consequence.

## Provenance

Before reconciliation:

- `result.yaml`: `aab660920b47bcea61bef7c98de6157ca6146a8a2bf308f583c719c96164ebc4`
- `source-audit.md`: `89b969d3ff5e327eb8b3f74535512ef0dd8eb06b0abd1ad2a385a67b735e7f94`
- `validation.md`: `c0edb8bc45c55e851841a9828b75e9b10b72a2fc803b7a4338cc1520ab738d7a`

After the vorticity repair and reconciliation:

- `derivation.md`: `55b51d2769b821848ddbbd64be70df0a4b85320cda2e808d6e98c0c82fe4bea5`
- `result.yaml`: `f63e8797bcc74f9d32a55689eacee4389c4a0b36188c2ba845d17778bf4104c6`
- `source-audit.md`: `8bec4f6e251f85e61c30286062832e3dcfe2380e1efdce8a201890a0e3a49eaf`
- `validation.md`: `f98a7a28a1753ef25af35515898761fb420f98c590e2eb19860fc1b9dc78be0d`

The first execution attempted the nonexistent worktree-local interpreter and
is preserved as `weighted-tail.{command.txt,stderr,exit}` with exit `127`.
The corrected repository-interpreter run is append-only:

- command: `feed74fe616631f5f60ad0085ba7a86f4ef66c749931aa3d332152cbb05bd9f4`
- stdout: `4077a8b10d7a52caa21bd12014e6f11ea497c2f96fd2bbda48e7d5b0cb868d11`
- empty stderr: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- exit `0`: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`

The focused public-test receipt reports `12 passed` with exit `0`; it checks
the algebraic APIs, not the written Banach estimate.
