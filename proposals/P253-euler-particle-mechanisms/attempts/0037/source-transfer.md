# Source boundary for the whole-space long-wave reduction

Two primary historical results were checked after activation to test whether
the threshold calculation was already available as a transferable theorem.

S. Leibovich, *Weakly non-linear waves in rotating fluids*, Journal of Fluid
Mechanics 42 (1970), 803--822,
DOI `10.1017/S0022112070001611`, derives KdV for a rotating fluid in a tube.
Its abstract explicitly says that the method fails as the wall moves to
infinity and that a singular perturbation produces an integro-differential
equation. It therefore cannot be imported as a whole-space Euler
justification. The exact Bessel logarithm and scale (4) in this attempt are
the mechanism that must handle that singular limit.

F. J. Higuera and J. Jiménez, *Solitary waves on a vorticity layer*, Journal
of Fluid Mechanics 264 (1994), 303--319,
DOI `10.1017/S0022112094000674`, constructs a different contour-dynamics
axisymmetric family and reports a KdV small-amplitude law with logarithmic
corrections. Its vorticity class and dynamics differ from the smooth
pure-swirl column fixed here. It corroborates that logarithmic whole-space
corrections are physical, but supplies neither the `0027` field nor the
uniform full-Euler complement estimate.

The source verdict is therefore **blocked as a direct transfer**. The exact
leading equation in `threshold-reduction.md` is derived from this campaign's
own Bessel exterior, eigenvalue slope, nonlinear projection and exact branch.
The still-active achievement is a uniform evolution theorem on the physical
time `L_mu/mu`; neither historical source closes it.
