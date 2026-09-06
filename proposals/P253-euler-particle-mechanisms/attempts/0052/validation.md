# Validation and oracle boundary

## Exact execution

The first execution used the repository interpreter:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0052/verify_fixed_domain.py

The exact command, complete stdout, empty stderr, and exit `0` are preserved in
the `fixed-domain.*` quartet. Seventeen exact predicates pass. No production
numerics or scripts-pane job was used.

## Exposing checks

The verifier checks:

1. both distinct corrections in
   `n_*=A/(c delta L)-Ad/(c^2 delta L^2)+(B/A-beta_0)+...`;
2. preservation of the `O(delta^3L)` coupling,
   `Theta(delta^3L^2)` spacing, and `O(1/L)` ratio;
3. the constant physical Jacobian relative to the fixed polar volume in the
   action coordinate `I=I_a s^2`, without asserting equal physical volumes
   across `delta`;
4. the determinant in the contravariant Piola divergence identity;
5. the canonical `-i` KKS sign and exact Hessian/frequency identity;
6. both sides of the inverse formula for the triangular local
   transport/shear block;
7. the two modewise `HJ2` remainder-to-gap ratios `epsilon/L` and
   `delta epsilon`; and
8. the leading, logarithmic, and finite `rho_G` coefficients obtained from
   the exact elliptic-integral representation of the Cao Green kernel.

The check exposes four false shortcuts directly: the old `O(L^-1)` crossing
remainder, a determinant-free Piola map, the positive-`i` 0048 scratch sign,
and an unstructured `o(delta^2)` graph remainder compared directly with the
smaller resonant spacing.

## Analytic boundary

The action--angle/Hanzawa--Moser/Piola construction, physical energy metric,
all-sector triangular/Grushin reduction, and conditional Riesz implication are
analytic arguments in `derivation.md`. The verifier does not replace them by
a sampled matrix.

The full physical graph-Riesz theorem is conditional on the explicitly stated
`HJ2` common-domain/modewise global-Euler jet. Cao's published moving-support
elliptic remainder does not prove `HJ2`. The 0048 fixed-order residual series
remains approximate and is not validated as an exact branch.
