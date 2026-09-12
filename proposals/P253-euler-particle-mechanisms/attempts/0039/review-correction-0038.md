# Bounded torus-frequency correction requested by 0038

Independent review `0038` observed that `exp(i phi/h)` is not globally
single-valued for an arbitrary real `h` when `phi` is an action angle.  The
packet sequence now uses

    exp(i N phi),  N in the positive integers,  h=N^-1.

For each fixed circuit count `j`, the construction first chooses the tube
width and then an integer `N_j` above the same finite WKB threshold.  The
integers are unbounded, so every estimate and the weak-null diagonal sequence
remain available.  No density argument in a continuous semiclassical
parameter is used.

Before correction:

- `README.md`: `be8db67ef25fc542ef5721f87958386041ad1f70c8fe7ca95157076c54943fe4`;
- `derivation.md`: `0a47bea9732796586208d31f38490eee16e2d72c8a4a315ffee4e7a7bdea4fbb`.

After correction:

- `README.md`: `c4f51f772e5eb17c134df23a47304560944693e828594ab761a1e9b654374b94`;
- `derivation.md`: `4f78aa4a3f2057c2c883c17a885c13d76a6b95d64f44f4807e99af471cf90845`.

The physical packet, accessibility, relative remainder, per-circuit
quantifier, essential-norm bound, and linear-only verdict are unchanged.
