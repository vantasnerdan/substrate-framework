# Physical energy-action normalization correction

Sol-High strategy audit exposed that KKS normalization alone does not select
the complex coordinate used in the Schwinger map: reciprocal symplectic
rescalings of `q,p` preserve `Omega`. Root corrected the author attempt before
independent review.

The corrected construction requires the positive Hessian/compatible complex
structure. In a diagonal normal form it uses
`Q_a=sqrt(nu_a)q_a`, `P_a=p_a/sqrt(nu_a)` before
`z_a=sqrt(B/2)(Q_a+iP_a)`. It also records the exact
complex-linear/antilinear decomposition
`V_C=(V-JVJ)/2`, `V_A=(V+JVJ)/2`: a physical analyzer needs two
number-preserving noncollinear Pauli axes and a gate-time bound on squeezing,
action drift and complement leakage.

| Artifact | Before SHA-256 | After SHA-256 |
|---|---|---|
| `0051/README.md` | `3d31ab6820254ba09c60ab1fde20fa7501bc732df5a5a890fd33768869d46d4d` | `321f8429919a3b6c201693ff2220ef00114fcf44958dc910d0a00be6420c448b` |
| `0051/construction.md` | `78802a7206506a7daa8ee012a9362cb746695adcebc9dac098c94976f995879c` | `8df69d4c88f675af7fae4446656856d1b5f4b5069c544bff64bc3bdf8024d25d` |
| `0051/source-audit.md` | `dc3d53c8b90f5bc83079896e885d97767c62294a39096ad5983c1e1894b70fe3` | `e355dd7046e93595ff1f8062f391bd9f28380e11391ebe240c382408d11db50b` |
| `0051/result.yaml` | `4fdb8fe64ae135869fdc6f21553836a94bf1f04903537e2e681c69674949e65b` | `a889c0fb325bdb08db640dcc35dd5244c8a54be8d2db8fd471db48ce145dbe07` |
| `0051/validation.md` | `6bf6cb12f4807518754b1fd475df572046d21b272e5f3bb501e48b4c76a3bff1` | `3cd41bd877efcdc881be28eca92b1784c10caa53beda8729a181a7fe1be0fae8` |
| `euler_schwinger_hopf.py` | `7bc7ead5a9a8eb800028a17a408e08a7a042d23d2ccb4fe2c74afda1932bd3d9` | `3e3a9981e8be6ca7f21cc6e30087dca2924c624e40c4ac019f08d7dc5b6f9e4e` |
| `test_euler_schwinger_hopf.py` | `02ac131ea1015581badd0835636c3f18755997bf786683b384219faace97f402` | `9868f214c652d8629b27ea7fe524e6fdd344e9d7022cd31bb4137c3fac257ee6` |

The strengthened corrected focused run in the dedicated Herdr scripts pane
puts unequal positive frequencies directly into the symbolic Poisson oracle
and reports `14 passed in 2.00s`, empty stderr, exit `0`. Its stdout SHA-256 is
`ef3f62b56e136ea2c6fe5868700df626f91d87ca9ec74e69f125119686e4eae7`.
The exact abstract Hopf/Stokes algebra is unchanged; the correction closes its
physical normalization premise and sharpens the Euler analyzer dependency.
