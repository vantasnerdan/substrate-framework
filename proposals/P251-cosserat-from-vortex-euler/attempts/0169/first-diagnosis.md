# First exact-verifier representation failure

The first run passed the actual nonadvected-pressure evaluation, its
independence of the Bernoulli-lift parameter and the divergence chain
rule. Direct trigsimp of a high-degree multiplied Euler residual then
left an unresolved symbolic zero, so the verifier correctly did not pass.

diagnose.py, captured on first execution in diagnosis.stdout, expands
the exact product rule before separately simplifying the lower-degree
stationary Euler residual. Both differences are identically zero.
The repair uses that analytic factorization. No field, coefficient,
selection criterion or tolerance changes; no numerical near-zero
decision is introduced. first.stdout remains the native failure record.
