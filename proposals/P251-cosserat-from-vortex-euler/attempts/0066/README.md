# 0066 — complete joint-action jet and physical-angle pullback

Parent P251 / issue #198; assigned owner `/root/construction_review`; this
directory only. Frozen task: derive the complete time-even, parity-even,
isotropic translation/axial-angle mass and stiffness jets, and their common
normal-form field map. Inspect whether the separable `(U,q)` action of 0059
and exact determinant identity of 0061 supply the original N2/N3 joint
Cosserat action at the retained second-gradient degree. Preserve all terms
at that degree, including mixed stiffness and terms induced by the physical
map `Phi=q+curl U/2`.

Fixed theorem route: exact polynomial matrix congruence, Hermitian symmetry,
objectivity and determinant expansion. No empirical comparator, fitted
coefficient, numerical eigenvalue, all-wavelength invariant-subspace demand,
or new parent acceptance condition. The concrete protected failure is
manufacturing coupling by dropping terms of the same retained derivative
order after a field transformation. A genuine missing joint term will be
identified and its full coefficient map repaired here; canonical files stay
unchanged. 0059's owner authors this continuation, not its independent review.

Oracle: exact symbolic congruence and acoustic/optical root expansions,
with mutations that omit an induced term or decouple a genuinely mixed input.
Original N2/N3 premises and 0061 are source inputs. The local proof decides
the coefficient identities; microscopic realization of any newly required
coefficient remains a named parent construction, never silently supplied.

Completed receipt: `joint-jet.md` establishes the exact complete jet and
physical map. The induced relative-curvature pullbacks start at k³/k⁴, so
they do not invalidate the formal second-gradient Cosserat normal form.
However the physical centroid optical transfer is controlled by
`l=g-kappa*b/j`, which vanishes for the separable 0059 action. This is the
specific remaining joint response, not a missing same-degree curvature term.

`verify.py` imports the new canonical `relative_angle_field_map`; first
execution preserved in `stdout.txt` passes 20/20, exit zero, 1.765 seconds.
Ruff passes. Exact mutations expose dropping the kinetic cross and claiming
an exact determinant from only its second-order jet. Route verdict:
established as stated; evidence scope: exact conditional joint-jet theorem,
not a positive claim that the zero-transfer specialization meets the parent's
requested physical coupled response.

Failure-generated continuation 0070, assigned by the parent: use the eleven
dual responses to retain and prescribe rotational affine KKS pairings instead
of annihilating them. This constructs and computes a nonzero relative kinetic
connection from the same Euler action; its action and locality are the next
obligation, rather than assigning an otherwise forbidden potential cross.
