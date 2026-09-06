# Validation and evidence boundary

The construction analytically derives the Koopman commutator, KKS integral
class and polarized dimension, exact Euler action scaling, both `Z2` exchange
characters, the conditional Hopf/FR parity, and the bare uniform-Euler
convective frequency. The new API makes the selection inputs explicit rather
than hard-coding a fermionic sign or `N=1`.

Focused repository-interpreter command:

    /home/dan/substrate-framework/.venv/bin/python -m pytest -q \
      tests/test_euler_quantum_bridge.py tests/test_euler_field_maps.py

Result: `11 passed in 2.72s`, exit zero. The nonintegral-Chern and invalid-deck
mutations pass by rejecting the false inputs. The scaling test exposes the
continuous action family, and the exchange test proves that both characters
remain available.

Pinned SHA-256:

- README: `3e4a8dd16a28a51f40676e148c3a92231338412f6abb30bca4fb20a680e53675`
- source audit: `4d4a7ab43b057c487d0da788b5badf9ee33b66fb3207b736ca386070eb1e1415`
- construction: `828b5d26b12dad5b63956181f608253b4c6c0c933ea0cc112ab962a029fa0e64`
- API: `cc5035538ba233e255fe73b495b6b2ef53b924e69aadc9102c3f451c69708f8e`
- tests: `1e43e97a1b3fdcca8d3d25dc94382f38e9dd24713f02283b3027602682a6f6b6`
- focused stdout: `1b564d56b190fb6a550bb038ea4ced8533b5574ea4b3a63e6c41a456b4a2630b`
- focused exit: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`
- Krusch--Speight PDF: `b92709003a959ded272285525e35a56d53f1c297d9ad342607c37dfcdc30750d`

The tests verify the finite algebra only. They do not prove a persistent
carrier, an Euler-to-Hopf configuration-space injection, a Born rule, or an
emergent Lorentz sector. Independent scientific review and downstream replay
remain separate achievements.
