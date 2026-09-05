# First execution and exact implementation repair

The first invocation used
`env PYTHONPATH=src .venv/bin/python proposals/P251-cosserat-from-vortex-euler/attempts/0217/verify.py`.
It exited 1 after nine passing physical identities because the canonical
physical_scalar_chart API requires the separate measured-spin row. The
argument was omitted, not defaulted by the API. The native first output
is preserved in first.stdout (2.039 seconds).

The repair supplies the explicitly leading matched spin row
beta/(nu c^2) times the actual angle-rate row. This is the limiting algebra
already derived from the physical tag fraction in the proof, not a claim
that the API itself proves microscopic matching. The final energy check
was also strengthened from a scalar factor inequality to the full actual
coordinate pullback of the derived conserved energy matrix. The repaired
run exits 0, all 14 checks pass, in 2.976 seconds; its first passing output
is preserved in repaired.stdout. No numerical tolerance is used.

The source prose also makes the global/local meridional orientation
explicit: y_local=-z requires S=-R A_0 exp(-i theta-i varphi) in the
global curl convention, so the actual local field is the column used
by the exact verifier. Both quadratures receive this same sign; phase
and energy are unchanged, while the signed moment coefficient is then
the one displayed in the proof.
