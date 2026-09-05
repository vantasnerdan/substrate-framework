# Exact temporal stress receipt

First command: PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python
proposals/P251-cosserat-from-vortex-euler/attempts/0175/verify.py.
Native first.stdout:14/14 exact checks, exit0.

Continuation command uses higher_time_jets.py in the same environment.
Native higher-first.stdout:9/9 exact checks, exit0. This run changes the
analytic observable to the fourth time derivative and independently
contracts the isotropic tensor; it does not repeat the prior green tally.

The actual two-wave response has R_D''''(0)=-2/25 at unit amplitudes.
The analytic amplitude and symmetry argument extends this to
-A²*B²*(A²+B²)/25. All Fourier products are complete, all numbers exact;
there is no numerical spectral truncation or timestep.

SHA256:

- verify.py: 553ffdfb60d28209265fed5e3e47461fe02afcad9fa517779e7688ba1aa4a1dc
- first.stdout: ac34e21b9e86e056c1bba55ab62e0857ce2388e1bd48aa0776e2b75866b3cc53
- higher_time_jets.py: 597685a78dd3c18a53d9688160b2312f548162311a8d89199a1c5df4a2f3e855
- higher-first.stdout: 62168001cc4e795ef9ae478c757ecea690e5c71811d881d59f2e4481f710f66d
- temporal-stress.md: 52cdc3db9703c874a404d4a193e636e720b28ad5afaa6e789895f9d1c6e2b7ee

Targeted Ruff passes. Route verdict: refuted for automatic exact temporal
cancellation under stationary constant-curl/isotropic structure with this
bare preparation. Evidence scope: exact initial-time physical response,
not a no-go for corrected preparations or the full campaign. The next
actual stationary/current-corrector route is registered in0179 and executed
alongside0176's stationary optical observer construction.
