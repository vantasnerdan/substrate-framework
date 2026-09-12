# Single bounded precision correction for independent review 0047

The non-author review preserved every route verdict and requested one bounded
precision package. The corrected artifacts now:

1. make the Koopman representation conditional on a specified global
   invariant phase space, measure-preserving flow, invariant probability and
   normalized state, with bounded multiplication observables or explicit
   common domains for unbounded ones;
2. use the inherited `i_X Omega=dH` convention explicitly as
   `Omega(Ax,y)=K(x,y)`, define `J=A(-A^2)^(-1/2)` on a declared skew-adjoint
   real Hilbert domain, and use the positive metric
   `G_J(x,y)=Omega(Jx,y)`;
3. scope the Krusch--Speight result to the proved spatial-rotation class
   `Q mod 2`, requiring a separate physical exchange-loop transfer before the
   same phase is assigned to exchange;
4. type nonzero uniform whole-space flow as periodic or relative-background,
   leaving rest as the absolute finite-energy symbol; and
5. reject `j=0` in the KKS-sphere API because the rotation orbit and two-form
   collapse to a point.

The API behavior changed only at the previously false zero-`j` edge. The same
focused test inventory was rerun with an exposing zero-`j` rejection and
reports `11 passed in 2.70s`, exit zero.

Pre-correction SHA-256:

- construction: `828b5d26b12dad5b63956181f608253b4c6c0c933ea0cc112ab962a029fa0e64`
- result: `5f4cfe8d6b140accb12ad2baa895404c67e8167c3e8bf57ce4b5f6b0c491593f`
- source audit: `4d4a7ab43b057c487d0da788b5badf9ee33b66fb3207b736ca386070eb1e1415`
- validation: `88733482406946fee6613fa574b7fe8391f8637613f6e8500cee2ce7f89b7ead`
- API: `cc5035538ba233e255fe73b495b6b2ef53b924e69aadc9102c3f451c69708f8e`
- tests: `1e43e97a1b3fdcca8d3d25dc94382f38e9dd24713f02283b3027602682a6f6b6`

Post-correction SHA-256:

- construction: `0558d180f013c2c15f0223fb65a5fae60cc4ae7427132814862ec46fff0c4ea1`
- result: `a0f3fcb7b627de1434e3761b4d22b29a7bd2892f778f67f5f4e2ff9ae0e17efb`
- source audit: `7acfe4f8454e616efa55df3a68cc4caae067553e9184967382374f4cf4e119de`
- validation: `bb97402f9418831da7c51f87cd2853e670b9526c9e71aa8bd9273648ef5b8a61`
- API: `b773d255b4adbf1419a3c8a94cbcbadeb29386aedd4edb0bff33c75a408c19f4`
- tests: `979218834d05afc39dfc1178221418cf1aa8395daf4f1db65d0f54e87890d90c`

Focused correction receipt SHA-256:

- command: `dd0baad1c919492417ce6ed66564aade29c992fe7e430db97bd2f015d13b2f18`
- stdout: `a2d60f5309db0bc78c093a7442afc289bccb1d2826f0e993fc4d0896fc659aa6`
- empty stderr: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- exit: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`

The README remains the frozen preregistration at SHA-256
`3e4a8dd16a28a51f40676e148c3a92231338412f6abb30bca4fb20a680e53675`.
The corrected result YAML parses and the scoped whitespace check passes.
