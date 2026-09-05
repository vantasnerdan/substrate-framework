# Source-transfer receipt for the 0253 construction

This receipt records theorem scope, not a substitute for the construction in
`construction.md`.

## Smooth localizable fields

[Constantin--La--Vicol, arXiv:1903.11699](https://arxiv.org/abs/1903.11699)
starts from the Grad--Shafranov equation plus the additional speed first
integral and obtains a localizable smooth compactly supported Euler field by a
pressure cutoff.  Its source construction omits the singular central point
when the cutoff is chosen.  It licenses the exact cutoff identity used for
Candidate A, but not a collar matched to a prescribed nonlocalizable nonzero
core.  The archive records `0143` and `0169` check this same scope.

## Compact but nonlocalizable fields

[Dominguez-Vazquez--Enciso--Peralta-Salas,
arXiv:2005.04380](https://arxiv.org/abs/2005.04380) constructs compactly
supported, axisymmetric, nonlocalizable stationary Euler weak solutions using
an overdetermined Grad--Shafranov boundary problem.  The velocity is piecewise
smooth and discontinuous across a surface.  This is primary evidence that a
nonlocalizable toroidal matching problem is coherent, and that its boundary
data are load-bearing; it does not license the `C-infinity` velocity, pressure,
or flat free boundary required here.

## The logarithmic compact-support threshold

[Pucci--Serrin--Zou, *J. Math. Pures Appl.* 78 (1999), Theorem
2](https://pucci.sites.dmi.unipg.it/lavori/psz.pdf) proves a compact-support
principle for nonnegative solutions of

\[
 \mathop{\rm div}(A(|Du|)Du)-f(u)\geq0
\]

which tend to zero at infinity, under positivity near zero and the finite
integral condition involving the primitive of `f`.  For the Laplacian this is
`int_0 ds/sqrt(int_0^s f(t)dt)<infinity`.  Their Theorem 20 treats uniformly
elliptic variable coefficients under its stated comparison hypotheses.  For
`f(s)=s[log(1/s)]^p`, the primitive is asymptotic to
`s^2[log(1/s)]^p/2`, so the integral is finite exactly when `p>2`.

This theorem licenses the threshold and says that an already-existing
decaying solution of the applicable inequality becomes compactly supported.
It does **not** construct the present Grad--Shafranov solution, control the
topology of its positivity set, preserve the constant-curl core, or prove a
finite-radius toroidal inverse.

[Diaz--Hernandez--Ilyasov, DOI
10.1016/j.na.2014.11.019](https://doi.org/10.1016/j.na.2014.11.019) and their
[autonomous follow-up, arXiv:1808.03931](https://arxiv.org/abs/1808.03931)
provide variational existence and compact-support results for autonomous
strong-absorption spectral equations under their domain, dimension, exponent,
and stability hypotheses.  They support the variational/free-boundary route,
but their power nonlinearities do not by themselves give a physically
`C-infinity` logarithmic compacton, the weighted `Delta*` operator, or the
large-ring continuation required in 0253.

Finally, [Soave--Weth, DOI
10.1137/17M1144325](https://doi.org/10.1137/17M1144325) proves unique
continuation for a class including `-Delta u=|u|^(sigma-1)u`, `0<=sigma<1`.
That sign is not the 0253 flat-edge sign.  Near the proposed edge the equation
is `Delta psi=positive_coefficient*psi*log(1/psi)^p+lower_order`, precisely the
compact-support-principle orientation.  Non-Lipschitzness alone therefore
licenses neither compact support nor unique continuation; the equation sign
and the integral criterion must be retained.

## Net source license

The sources license (i) the exact pressure-cutoff operation after an actual
localizable collar has been built, (ii) the mathematical relevance of an
overdetermined nonlocalizable toroidal problem, and (iii) the sharp `p>2`
flat-edge threshold.  No located primary theorem supplies the missing smooth
finite-radius toroidal continuation with the prescribed nonzero Beltrami core.
