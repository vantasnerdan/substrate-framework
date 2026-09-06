# P253/0031 bounded Hodge-domain correction

## Finding and minimum repair

The original equation (18a) used compact support of `delta u` to identify the
vorticity impulse with the translation momentum pairing. For an actual compact
coadjoint generator `eta`, only
`delta omega=curl(eta cross omega)` is compact; the recovered velocity
`delta u=P(eta cross omega)` generally has a noncompact Hodge pressure tail.

The corrected derivation integrates by parts on the genuinely compact field
`eta cross omega` and proves

    delta I_z=rho_m integral (eta cross omega)_z
             =rho_m integral eta dot partial_z u
             =Omega(X_Z,X_eta).

For axial rotation the Leray gradient contributes
`integral partial_theta p=0` on each complete azimuthal cylinder. The finite
local KKS form, weak Hamiltonian identity, exact translation orbit, observable
values, asymptotics, and both route verdicts are preserved. The wording no
longer puts a general noncompact Euler generator in `Diff_vol,c`; construction
of an arbitrary-nearby asymptotic group remains conditional.

## Content-addressed boundary

| Artifact | Before SHA-256 | After SHA-256 |
| --- | --- | --- |
| `derivation.md` | `8a44074412f41b3fcee94dcd5154f87ab060adca51024af29fb38adc541002cd` | `f92079c746fd66a1b182c85fb5da3f1b10dbe9772e2ae22ff6d094f5e778d12a` |
| `result.yaml` | `e2bde6a06bca0dd968db1032e0d708e6d97aa0184e790d8650ffd4e8c437298e` | `9041e385d75de3fce6a71eeff64c9cd8c5449be5c04c704ac9cb7d3a7411fb42` |
| `validation.md` | `1918223d27c09888573369c84c60308de2d025f436fc921c08185aec6f1622f6` | `0e384df0f14019ab68d18df8683caa42bf1329975a696a158e7dc8d76735e916` |

The verifier remains
`08c6f7bd218d485cda047188648a11bb236ace2bdc1c81bdbbce53e0a1ba590a`
and its retained stdout remains
`c5b42dce0f83519f22ee4b9978698d56f4525a788dd08d1bade8c9ed107ddea2`.
No equation implemented by that verifier, API, test, or supplier changed, so
the exact oracle was not rerun.
