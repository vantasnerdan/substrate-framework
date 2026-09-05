# 0268 independent finite-R inverse review

Reviewer: `herdr optical-review pane w3:p3`, a non-author process.  The
0263 geometry author and 0261/0256 authors are distinct; my 0250 work is not
used as independent evidence.

Pinned inputs: `trace-estimate.md`
`c13a1066d495e67df610ea0d64df0c2115808ab2ee42261a238a1dbd4a82b4e4`,
`perturbation-gluing.md`
`445f91702b63d0af92c4d91a7e5041578ac6b88b0a0cb1e064e540b40758a88e`,
`source-transfer.md`
`47fc7c0ccbbfc3f64eb2920deed3efe646e2bc99c0a2b25af446f39720d7b0f3`,
0261 radial construction
`55bb250a03963c097d8c918c854610591ff5e56e2f4d275ae561ba98d889803c`,
0261 inverse analysis
`faa86a06b3b1c57d76ad7400bd7fff3065358f497cd9f22cdf6fb82385040752`,
and 0256 balance
`b3b5d018cf844eceff9e64f1c3747236927c534bf8d8e92221ae2c0065ab604e`.

## Audit

The radial trace theorem retains the full normal second derivative.  The
Stratonovich-to-Itô correction is included, the defining-function drift is
uniformly inward, and the entrance/stopped representation uses only inner
exit data; it does not prescribe both Dirichlet and conormal traces.  The
Fourier barriers retain `d^2 v''` and give the stated exponential inner-data
factor and `|m|^{-2/3}` source gain, with low modes handled separately by the
recorded radial gap.  The moderate branch exclusion identifies this trace
with the physical form-domain branch.

The two-sided extension is a genuine construction: the half-collar is
invariant, so arbitrary negative-side data are not sampled, while the
extended solution is interior hypoelliptic.  The collar shrink makes remote
exit data enter only through a low norm.  The displayed intrinsic estimates
and tame product bound provide the parameter-uniform entrance resolvent
needed for a small Hanzawa perturbation.  The global form solution supplies
the matching traces automatically, avoiding an overdetermined interface.

Finally, the perturbative right inverse, translation projection, scalar
Hadamard division, and Grad--Shafranov reconstruction produce a finite-R
smooth compact Euler field with smooth zero extension and the stated fixed
annular twist.  The 0256 balance and 0261 nondegeneracy are consumed only at
their recorded scopes.

## Verdict

`route_verdict: established as stated` for finite-R compact Euler existence
in the intrinsic projective-limit scale, under the explicit radial spectral
gap, smooth collar/Hanzawa-smallness, entrance accessibility, and profile
nondegeneracy hypotheses.  The remaining mechanism is outside this review:
same-field response/action/current and hybrid observation transfer.  No
trace/domain/parameter or physical-gluing correction is required within the
finite-R boundary.
