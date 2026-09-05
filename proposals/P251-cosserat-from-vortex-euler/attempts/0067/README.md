# 0067 — exact physical relative-angle field map

Main owns this bounded canonicalization of 0061/0066 algebra. Add the
importable map from physical (U,Phi) to (U,q=Phi-curl U/2), without asserting
that a coordinate change constructs a new microscopic force. Both action
forms use the same congruence. The exposing test retains the third/fourth
gradient entries and verifies exact spectral factorization in the separable
case, both helicities, and the mutation deleting relative kinetic coupling.
No discretization or numerical small-ratio sign inference is used.

Impact: additive unpromoted helper in micropolar.py, direct new test in
test_micropolar.py and forthcoming 0066 consumer. Existing consumers are
CST004 (energy), CST005 (Fourier stiffness), and test_micropolar.py, located
by direct rg. GitNexus at c9e1c0f locates the Fourier function but lacks its
test/attempt call edges; that is not evidence of no consumers. No existing
function, accepted claim, formal proof or generated documentation changes.
