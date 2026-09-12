# One bounded 0050 target correction

The independent `0050` review requested one bounded precision package.  Root
changed only claim-bearing prose in the frozen `0046` and `0049` attempt
artifacts.  No API or test behavior changed, so the existing focused receipts
were not rerun.

## 0046

The KKS oscillator now states the complete positive-definite Hessian
hypothesis (`Hq>0` and `Hq P-N^2>0`).  The finite-dimensional trace argument
is scoped to an exact CCR truncation of the noncompact canonical plane and
explicitly leaves a separately constructed compact `su(2)` orbit/reduction
open.

| Artifact | Before SHA-256 | After SHA-256 |
|---|---|---|
| `0046/construction.md` | `744968cbedfaa2e1dfd424eec3e7d2ca8e98d303d481b76a88240bed66608b1f` | `7ee4b9d9679a9e2880b49a1efa4294a46f0c88be573124074dc39b55985a6d63` |
| `0046/result.yaml` | `7bc3acdf32661ea2dfc4eb8a23f10541bf7f3bffe647d8bb9c99724eb829d9ea` | `698841744a5ac7e1ed9a5e8a57fe3669cecbe31fe35b0338416293b303df5f69` |
| `0046/source-audit.md` | `77987e351b920ef74d0a06128f716a1546c07c00cc7c8c882a960ef4883a6147` | `74accea1e7b3a11410c8cf38c536936b36990ff517f61af47d57402cc3d671c2` |
| `0046/validation.md` | `bb85f0e02a4b5d05e1aaadb226505588f894b1f447a1fc1248548932a1e23582` | `c47ca3107b8fe6079f33758b0e3a62a2ddcab9ccb5a2c3619e2b237f15cac86b` |

## 0049

The deterministic threshold convention is now `lambda<=I` with
`F(x)=Prob(lambda<=x)`, so atoms obey the displayed formula.  The all-zero
clock boundary says that no finite event occurs and the conditional winner
law is undefined.  The projective reset is explicitly an open-system,
coarse-grained, or external repreparation hypothesis; a reversible autonomous
first-exit/scattering map must retain escaped degrees and does not itself
reprepare the sequential experiment.

| Artifact | Before SHA-256 | After SHA-256 |
|---|---|---|
| `0049/construction.md` | `d921309b780f44a8fee72ef723b6f11de9648a1806d71c34086edc68432a9897` | `b207eac0378998336d26af863b3467c6bbfc2f2cf9c7b6156427889c4b26a06c` |
| `0049/result.yaml` | `228cd30e6ad93eef33d647ef73c1a27498980623cb64d04b38a899dccc2478cc` | `b248708fa624b02bb478ae2ae02161ffc0f1b17a009120ac441dbaac40d3ca27` |
| `0049/source-audit.md` | `3eb37819c0c38767b46f5955a08f531fa19208e989d63cb7266f46f06bcdc970` | `f85f0dcd4aa831d010d05ce4ff0dbab9b7d50ae6974ba7369a9634430338c6f8` |
| `0049/validation.md` | `edc3cb99e471cb60e0af308817c2daaa8761e538f3457eae47619c131c82686a` | `282cc57296e77bbe2c3889a70bfc9311190d43090d1cc9c18a6d0ba32e096ae5` |

The stale-language scan found no old strict-threshold formula or unscoped
finite-CCR wording.  `git diff --check` passes on both target directories.
The classical KKS/Heisenberg boundary and conditional first-event theorem are
unchanged.
