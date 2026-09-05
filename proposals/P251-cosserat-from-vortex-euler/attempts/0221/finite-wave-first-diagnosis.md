# First finite-wave oracle support correction

The first script source hash was
`9f8e180da0026a38ec4967c7d45b18603cbbfd79a8589dc02dbdd53fca0af1f0`.
Its full-pressure energy calculation passed eight exact checks and
printed both complete rational H2 coefficients. The subsequent phase
nonzero-sensitivity assertion failed: the chosen a-sideband difference
does not couple to the background omega_a=-sin(b) row. Its exact phase
is zero by Fourier selection, not because the Kelvin phase is absent.

The repaired phase probe uses the b-sideband, preserving the actual
divergence-free Bloch generator and original energy calculations. The
new --phase-only selection executes only the previously uncompleted
phase/sensitivity/scale rows. The eight unchanged energy checks reuse
their first captured evidence. A normal invocation still executes the
whole final script. This is an oracle-support repair, not a change of
the physical construction, symbol order, threshold or campaign scope.

The first corrected b-sideband run also failed the nonzero-mixed-phase
premise. Its preserved phase-wave-repaired.stdout contains that result.
Here the support was correct, but the opposite c=-1,c=0 sectors have
an ACTUAL cancellation between the axial-vorticity and planar-vorticity
phase rows at K=0. Directly expanding their cross product gives equal
opposite contributions. Thus a nonzero mixed-sector phase is not the
right sensitivity predicate. The final probe keeps the second axial
fraction symbolic, tests the nonzero same-positive-sector phase and
the exact opposite-sector cancellation separately, and derives the
second spatial coefficient of the entire fraction-dependent phase.
This additional null is compatible with the zero-wave normalization;
it does not invalidate either sector's own nonzero phase. The actual
cross rows remain retained at finite K.
