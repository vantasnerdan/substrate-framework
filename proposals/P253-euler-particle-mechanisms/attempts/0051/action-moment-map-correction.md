# 0056 action/moment-map normalization correction

Independent review 0056 exposed one factor-of-two error in the otherwise
correct Schwinger--Hopf algebra.  With

    Omega=B sum dq_a wedge dp_a,
    z_a=sqrt(B/2)(sqrt(nu_a)q_a+i p_a/sqrt(nu_a)),

the diagonal physical oscillator has

    H_a=nu_a |z_a|^2,       J_a=H_a/nu_a=|z_a|^2.

Therefore `J=sum_a |z_a|^2` is the physical total action and common-phase
moment map, `{z_a,J}=-i z_a`.  The Stokes vector retains
`S_i=z^dagger sigma_i z/2`, hence `|S|=J/2`, the reduced KKS area is
`2*pi*J`, and prequantization is `N=J/hbar`.  The former convention called
`J/2` total action and would have doubled the area for a caller supplying the
physical action.

The public `total_action` and `reduced_kks_area` semantics, exact Hopf
residual, documentation, result ledger, and tests were corrected.  A new
`diagonal_mode_energy` API makes the unequal-frequency energy/action bridge
explicit.  The review's blocked Euler-supplier, mixing, leakage, action
selection, detector, exchange, and propagation verdicts do not change.

Frozen pre-correction SHA-256 values:

- API: `3e3a9981e8be6ca7f21cc6e30087dca2924c624e40c4ac019f08d7dc5b6f9e4e`;
- tests: `9868f214c652d8629b27ea7fe524e6fdd344e9d7022cd31bb4137c3fac257ee6`;
- README: `321f8429919a3b6c201693ff2220ef00114fcf44958dc910d0a00be6420c448b`;
- construction: `8df69d4c88f675af7fae4446656856d1b5f4b5069c544bff64bc3bdf8024d25d`;
- result: `a889c0fb325bdb08db640dcc35dd5244c8a54be8d2db8fd471db48ce145dbe07`;
- source audit: `e355dd7046e93595ff1f8062f391bd9f28380e11391ebe240c382408d11db50b`;
- validation: `3cd41bd877efcdc881be28eca92b1784c10caa53beda8729a181a7fe1be0fae8`.

First corrected post-run values before this receipt and validation update:

- API: `b91bbb6f39750cea2e846574c97d5c1f74ca86e91231b35db164869a25eea3c3`;
- tests: `185ecb8fe1cdcb99202e913d94ac28abb38d562677862ecf5cc7e5409575f65c`;
- README: `d8e7b92b388dc1a5f1efdbbf00913de32012b99ff58ef22a3fab0309f1cfdc50`;
- construction: `841ca0caf6516e5955b2de89630bcf05402cf54473773c70b679308a23cd0103`;
- result: `f47398cb566bd3d3b4edc6952b6f5682ae4ee8c1f3b4c2161fedc73471948931`;
- source audit: `82efd88350342fbcd726c0dfe262b35338ccd778b30fb16ad513ce6cd0f22da5`.

The focused repository-interpreter command reran the same 0051 plus direct
consumer inventory.  It reports `14 passed in 2.00s`, empty stderr, and exit
`0`.  The command, stdout, empty stderr, and exit hashes are respectively
`d8d75be024e6ecd2c6056ce7c1abcbd82e9fa2c7936b1c10172a25edde9a4f0c`,
`ef3f62b56e136ea2c6fe5868700df626f91d87ca9ec74e69f125119686e4eae7`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
