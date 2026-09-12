# Single bounded analytic-domain correction for 0045

The independent `0045` review found the pressure, split, history, Dyson, and
finite algebra correct and requested one bounded analytic-domain correction.
The correction now:

1. restricts every material observable to a transported bounded material
   region or compactly supported/integrable tag with the finite weighted
   moments and boundary regularity required by `0001`;
2. takes `v in C^1_t Ran Pi`, `w_0 in Q H^s_sigma`, `s>5/2`, with `Pi` a
   finite-rank `L2`-orthogonal smooth-range solenoidal projection bounded on
   `H^s_sigma` and smoothing `H^(s-1)_sigma` to `H^s_sigma`;
3. constructs `W[v,w_0]` by the forced-Euler transport/energy method on a
   common interval controlled by the initial and prescribed-history norms,
   and proves `Pi(u-v)` remains zero;
4. distinguishes the bounded observable idempotent `P_obs,Q_obs` from the
   velocity projection and states the common invariant-domain/generator
   hypotheses for the Dyson formula; and
5. retains the full ambient field generally while treating its integrated
   momentum as a scalar only under an additional `L1`/decay hypothesis.

The generic infinite-dimensional `Q_obs L` semigroup remains explicitly
unproved. The exact retained-state and causal-history verdicts are unchanged.
The API and tests did not change, so the existing 12-test receipt applies and
was not rerun.

Pre-correction SHA-256:

- README: `9c37f046edde581edd3e03e0c69fadaffb5d5ef91c7b1001c87bc1797c013a87`
- construction: `6e9993547fea6617150e4a7f8f72cefa75dc443026b4215f1e173e2767d1d277`
- result: `d23ebb23b0dc755ecb83926580ab49a932d689bfafa6e78f8037e80e6e8d6424`
- validation: `e8efa13afae3422c4e4526394c9dfb35da08fed4709fe66233c97a612758de28`

Post-correction SHA-256:

- README: `416423a720c879f9fd98666292305e2d9492b685fad424a4f53e437f7310fa84`
- construction: `b573d04f0d3a46819e569568c73c69c4d79a70cabb845f35e15dc56f93fa7061`
- result: `319b248ee41e84a9753198b57de016dcf98003c83744bb80d575ef8416a7ebd7`
- validation: `226bfd5cdf9330a3e2adc55780d5b52adbb2b8e091757153f13a437db7ee7934`

The corrected result YAML parses and the scoped whitespace check exits zero.
