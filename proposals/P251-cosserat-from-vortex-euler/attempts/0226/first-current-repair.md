# Explicit first toroidal displacement/spin return

The full 0222 calculation, not an identification of spin with a canonical
coordinate, gives the additional direct toroidal rows

    G_tor=-i rho pi^2 lambda R integral chi(s) A0(s)s^2 ds,
    S_vel= rho pi^2 lambda R integral chi(s)Omega(s)A0(s)s^2 ds.
                                                               (1)

These rows were retained there as finite-R errors. They are different
linear functionals from the leading physical tilt
B=integral s A0(s)chi'(s)ds. The new candidate executes their initial
and finite-clock-order removal rather than forgetting them during a
G/spin comparison.

For a chosen finite order q impose the additional actual profile rows

    integral chi(s)s^2[Omega(s)-nu]^j A0(s)ds=0,
                          0<=j<=q+1.                    (2)

Then both rows (1) vanish at the reference time, and their transported
versions have the Taylor remainder of order q+1 in the clock spread.
For G this follows directly from expanding exp(-i[Omega-nu]t); for S
write Omega=nu+(Omega-nu) before the same expansion. No stationary
geometric shape is inferred from a material label alone.

These rows can coexist with a nonzero leading tilt and the compact
pressure moments. Here is a concrete choice that verifies their
functional independence. Take the positive radial tag to be constant
near the axis, strictly decreasing on a selected annulus, and smoothly
flat at its outer support boundary, with
chi(s)=exp(-1/(b-s)) on a one-sided neighbourhood of that boundary,
joined smoothly and monotonically to the inner constant. All supports
are strictly inside the literal-curl region.

If the leading weight s chi' were a finite combination of the weights
in (2), it would obey chi'/chi=s P(Omega(s)) on that neighbourhood,
with a polynomial P. The right side stays bounded as s approaches b,
whereas chi'/chi=-(b-s)^-2 is unbounded. This is impossible. Polynomial
exterior-pressure rows cannot restore such a dependence: an open
off-tag part of the preparation support forces their polynomial
combination to vanish identically first. Consequently finite smooth
bumps solve (2), the pressure moments, and the nonzero leading tilt.

Including the frequency-flat leading tag rows only replaces the putative
identity by chi' P1(Omega)=s chi P2(Omega). Choose the centre frequency
and tag boundary so the required leading polynomial is not identically
zero there; alternatively, inspect its finite-order zero. The flat-tag
logarithmic derivative has a second-order pole while Omega has a
nonzero first derivative on the boundary annulus. To avoid coincident
polynomial cancellations altogether one can use

    chi(s)=exp[-exp(1/(b-s))]

near the boundary instead. Its logarithmic derivative grows faster
than every inverse polynomial. Then no two finite polynomials in the
analytic monotone Omega can satisfy that identity unless both vanish.
Use this latter positive smooth flat tag for the combined construction.
It is a fixed tag choice, not a parameter-dependent fitted weight.

The homogeneous kernel of all these finite linear rows still contains
nonzero compact profiles supported in either of the two frequency
intervals used by density-normalization.md. Therefore the disjoint
positive/negative quadratic-return construction survives unchanged:
H=nu beta is solved independently by its actual quadratic root, while
(2) removes the stated direct G and spin corrections.

This is an explicit physical-current repair of the two direct rows in
(1). Other first and higher curved Euler/Lin corrections belong to the
compact pressure recursion and its actual tag observations. They are
not silently equated with (1), and their full cancellation is not
claimed by this attachment.
