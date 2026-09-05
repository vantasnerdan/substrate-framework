# 0132 — actual Euler versus magnetic elastic-wave mechanism

Owner /root; exact source/derivation route only. Parent scope and accepted
claims unchanged. Positive target: find a valid large-scale restoring
mechanism for the actual Euler mean/rotation join, informed by the
nearest primary magnetic/Euler stability comparison.0129's simplest
Beltrami mean-translation cell gives negative stiffness, while the
accepted conditional Euler action does not claim free invariance.

Candidate approaches: actual periodic Euler pressure/cell response;
magnetostatic analogy with its physical field/action map derived; and
finite vortex-structure interactions within stationary Euler. Selection
criteria remain actual Euler generator, physical mass/observables,
retained pressure and reaction, correct signs and limiting behavior.
No magnetic force, external rotation, viscosity or fitted modulus may
silently enter the incompressible-Euler object.

Source inventory before opening body: H.K. Moffatt, Magnetostatic
equilibria and analogous Euler flows of arbitrarily complex topology,
Part2: Stability considerations, JFM166 (1986),359 onward. Primary
author archive: https://www.damtp.cam.ac.uk/user/hkm2/PDFs/Moffatt_1986_JFM_166_359.pdf
Search snippets suggest separate magnetic wave and ABC Euler stability
sections; those snippets are not theorem evidence. Exact source
hypotheses and the licensed dynamic map will be checked at the body.
No empirical comparator or numerical representation is selected.

## Executed result and continuation

`euler-magnetic.md` and `verify.py` independently derive the actual ABC
Euler coadjoint Hessian and complete instantaneous slow stress response,
alongside the distinct magnetic functional. Twelve exact checks pass
in `final.stdout`, exit0; Ruff passes. The counterhelical displacement
plane has positive energy but zero KKS rank, so it does not by itself
provide the missing positive-mass acoustic oscillator.

`first.stdout` preserves four passes before an incorrectly anticipated
anisotropic force failed; its actual full-pressure result is retained.
The repair derives that force independently in Cartesian coordinates.
The unused magnetic cross-term sign in that first script was also
corrected from the defining variation before its first successful check.
`repaired.stdout` has nine passes; the final extension adds independent
Poisson/stress checks and the magnetic-force mutation, for twelve.
These are exact algebraic executions, not numerical stability evidence.

Proof SHA256: `28028b362387ded1049bf4d1be5f050bf1a746c144005a612f4383bc9f1e5f1b`.
Verifier SHA256: `c2d11bc11f45fb0c34480fda0725803ace9b4ee2458f4969bb965da2b05862ec`.

Route verdict: the exact Euler response construction is established;
the distinct direct magnetic-to-Euler constitutive-transfer candidate
is refuted by the field/action mismatch. Evidence scope: exact
coadjoint energy and initial acceleration, not autonomous dynamics.
The next candidates change the actual microscopic response, in active
0133/0134/0135. No parent no-go or completion is claimed.

Coordination checkpoint posted without refreshing the terminal PR:
https://github.com/vantasnerdan/substrate-framework/issues/198#issuecomment-5549741146

The attempt-only checkpoint passes `scripts/validate.sh --fixed-only`,
captured in `fixed.stdout`:263 accepted claims,1048 valid memory records
and all fixed repository checks. The43 existing memory warnings remain
warnings. The canonical implementation/test boundary is unchanged from
dbc2bf3 and its2580-test full receipt; no duplicate pytest run is needed.
