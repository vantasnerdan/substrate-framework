# First-run implementation diagnosis

The first execution passed the three Cartesian divergence/acceleration/
pressure checks, then failed the specialized curl assertion. Substitution
of a composite applied function f(x²+y²) in an already differentiated
SymPy expression left its Subs(Derivative(...)) terms unevaluated; this
was not differentiation of the chosen rational velocity field.

The repaired oracle differentiates that velocity directly, independently
of the general-profile curl expression, and additionally checks all
general radial curl components. No field or scientific equation changed.
The original stdout is preserved as first.stdout; repaired.stdout is the
next execution, not a replacement receipt.
