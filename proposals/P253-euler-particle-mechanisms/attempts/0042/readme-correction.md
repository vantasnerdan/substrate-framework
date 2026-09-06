# Activated-target corrections

The first activated README had SHA-256
`2002c78268aa436824d5cf1beb7ae805a662da3c09ca86d42893df8301c47569`.
Direct derivation exposed two convention/domain defects before either was
used as a result:

1. With `(-Delta)^-1` denoting the positive inverse, the pressure is
   `p=rho*(-Delta)^-1 partial_i partial_j(u_i u_j)`.  The original displayed
   minus sign was wrong.
2. Projecting the entire Mori--Zwanzig equation by `P` kills the orthogonal
   initial term.  The corrected target is the full-observable Dyson identity.
   Its infinite-dimensional `exp(t QL)` is not assumed to exist for generic
   Euler.  The positive full-Euler construction instead uses the actual
   conditional unresolved solution `W[v,w_0]` for a smooth finite-rank
   solenoidal projection.

The corrected README SHA-256 is
`9c37f046edde581edd3e03e0c69fadaffb5d5ef91c7b1001c87bc1797c013a87`.
The corrected schema receipt has exit zero.  These corrections strengthen the
honest retained-memory statement and change no predecessor artifact.
