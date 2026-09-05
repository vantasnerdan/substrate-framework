# 0051 — exact material-cell pressure and spin balance

Main owns this attempt. Parent N4 needs coupled translation and angular
reaction from Euler, rather than merely variation of a prescribed micropolar
energy. Freeze a material partition D_a(t), its fluid mass centers X_a and
their exact linear/angular momenta. Shared-face pressure tractions obey
action/reaction. Derive distributional center-based stress and couple stress
by a bond-line identity, including the moment of noncentral pressure forces.

The microscopic Euler stress remains symmetric. Center-based force stress
need not be symmetric because its branch vector joins material cell centers,
while the pressure resultant acts across the actual common face. The missing
angular momentum is explicitly the material-cell intrinsic spin. No body
torque or fitted couple-stress coefficient is postulated.

This is a fixed exact weak-balance theorem, no empirical comparator and no
numerical remainder. Its oracle is integration against arbitrary smooth test
functions, checked by exact polynomial probes and wrong-sign/transpose
mutations. The subsequent constitutive identification of cell angular
momentum with the EPS orbit's impulse, including boundary terms, is kept
separate; this balance identity does not assume that identification.

The partition/center-based coarse-graining is explicit and differs from a
single Eulerian convolution of velocity: that latter filter can slave its
first moment to derivatives of the filtered velocity and does not by itself
introduce an independent material spin coordinate.

Completed balance receipt: pressure-balance.md derives exact center-based
force/couple stresses from shared-face Euler pressure. verify.py passes 9/9
checks, including arbitrary degree-five weak probes, angular stress signs,
total angular conservation and a solid-rotation counterexample to omitting
the impulse/spin boundary term. Ruff and diff checks pass. This child balance
theorem is established at its material-partition scope. Attempt 0052 now
executes the actual canonical-impulse/material-spin and joint-translation
matching; that missing identification has not been assumed here.
