# Exact compact-tag coefficient and mode-cutoff clarification

This addendum completes an explicit cutoff step in the just-frozen
`fixed-cell-construction.md`; no prior receipt is overwritten. The
ideal Gaussian Laguerre coefficient5/21 does not follow from an arbitrary
fixed cutoff with merely a small nonzero tail. Either retain its exact
cutoff-dependent value q_chi>0, or make the following finite scalar repair.

For m=5,n=8 let f=exp(-x/2)L_8^4(x), E=21, q0=5/21. The calibrated
tag coefficient for its actual fixed radial weight chi is

    q_chi = [integral chi x^(m-1) (x/2-E) f dx]
                           /[integral chi x^(m-1) f dx].

Begin with chi0=x times a smooth cutoff that equals1 on0<=x<=R and
vanishes beyond2R. The full chi=x integral gives q0 exactly. Its
truncated residual

    F(chi0)=integral chi0 x^(m-1)(x/2-E-q0)f dx

is exponentially small as R grows. Choose one fixed smooth bump psi
strictly inside a region with chi0>0 and

    F(psi)!=0.

Such a region exists because the displayed analytic integrand is not
identically zero. Define

    chi=chi0-[F(chi0)/F(psi)]psi.                          (1)

For one sufficiently large but FINITE R this is smooth, compact,
nonnegative and has nonzero angle denominator. It satisfies F(chi)=0
EXACTLY, so q_chi=q0. A common positive rescaling can make chi bounded
by1 without changing the ratio. This uses an explicit material-tag
shape integral, not a fitted Euler frequency or an empirical comparator.
The finite bump rank remains available on open intervals where chi>0.
Thus(1) supplies the fixed physical chi assumed in the exact5/21 and
5/7 coefficients of the construction. Without(1), the valid statement
uses q_chi and3q_chi, with positivity guaranteed by the same small-tail
choice rather than claiming the ideal coefficient is exact.

The perturbation MODE cutoff is a different object. Its Gaussian
polynomial is cut off at a fixed small PHYSICAL invariant-tube radius,
not at the tag's fixed scaled R. Its scaled radius tends to infinity
like1/ell, and the resulting Euler/Lin residual is polynomial times
exp(-c/ell²), including two relative carrier derivatives. It is smaller
than the retained algebraic remainders. The tag is allowed to observe
only part of this full-cell mode; KKS still includes the entire mode.
This distinction prevents a fixed omitted Gaussian tail from being
silently counted as an O(delta²) error as p grows.
