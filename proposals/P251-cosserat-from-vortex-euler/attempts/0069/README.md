# 0069 — second slow jet while retaining physical rotational momentum

Parent: the original Euler-to-Cosserat objective, unchanged. The 0066
physical-field audit finds that eliminating every affine KKS moment in 0059
also eliminates leading physical mean-displacement/spin mixing. Its exact
normal form is useful, but cannot manufacture that physical response.
Failure-generated replacement 0070 retains the rotational moment while
annihilating translations and symmetric tracefree affine moments. This
attempt supplies its spatial-locality license, not its physical action.

Route: split the first force moment into scalar and antisymmetric parts,
represent the latter by a compact curl, and reuse the explicit compact
double-divergence construction of 0057 for the remainder. The Leray
projection fixes the curl, so its mixed energy is local. No statistical
mixing, infrared regularization or isolated-patch energy assumption enters.
The frozen oracle is the exact moment and full kinetic-symbol algebra;
the regularity argument is analytic. There is no numerical remainder.

## Compact decomposition

Let F=xi cross omega on one bounded good patch, in centered coordinates y.
The three translation constraints give integral F=0. Orthogonality to the
five symmetric tracefree affine generators gives

    M_ij=integral y_j F_i=c delta_ij+A_ij,  A^T=-A.

Every antisymmetric A has a unique vector a with A_ij=-epsilon_ijm a_m.
Choose a fixed compact smooth unit-integral cutoff chi. Set f=-c chi and
V=a chi. Integration by parts gives, with all boundary terms zero,

    integral y_j partial_i f=c delta_ij,
    integral y_j (curl V)_i=-epsilon_ijm a_m=A_ij.

Therefore G=F-grad f-curl V has zero zeroth and all first moments.
Apply the explicit averaged-Taylor-center formula in
0057/slow-locality.md to G. It supplies a smooth compact B with uniform
controlled Sobolev norms and G_i=partial_j partial_l B_ijl. Thus

    F=grad f+curl V+div div B

is an exact compact decomposition of the ACTUAL selected force. Neither
the displacement nor its angular momentum is changed by this calculation.
In particular this representation retains A, rather than setting it to
zero to make the projection differentiable.

## Full stationary kinetic energy and its differentiability

Modulate each patch by exp(i k.X_a), constant on that patch. All three
compact potentials have the same coefficient. As in 0057, their stationary
Bloch amplitudes contain exp(-i k.y); their k derivatives are uniformly
bounded by the finite patch radius. Define at total spectral frequency p

    C(p)=i[p cross],  (D(p)B)_i=-p_j p_l B_ijl,
    P(p)=I-pp^T/|p|^2.

The zero-mode convention is that the internal velocity has zero mean.
Since both C and D vanish at zero, it causes no discontinuity here. The
projected velocity is exactly C V+P D B. Its complete quadratic symbol is

    [ C* C       C* D   ]
    [ D* C       D* P D ].

Indeed P C=C and C* P=C*, so the off-diagonal entry is a degree-three
POLYNOMIAL, not a singular projection derivative. C*C is a degree-two
polynomial. Only D*P D is nonpolynomial; it is homogeneous of degree four,
smooth away from zero and C^3 at zero, as proved in 0057. Consequently the
entire energy symbol is C^2 with its second derivatives bounded by a fixed
polynomial in |p|. Uniform good-patch Sobolev bounds license differentiation
under the stationary spectral integral in operator norm, including arbitrary
square-integrable reaction amplitudes and atoms of their spectral measure.

This retains all ambient velocity and inter-patch interactions. The mixed
term is local because one ACTUAL summand is a curl, not because cross terms
were dropped or parceled into isolated energies. The local helicity and KKS
terms have smooth finite-support Bloch factors. A uniformly coercive full
momentum block therefore has a C^2 inverse; the complete second-order Schur
jet of 0062 applies unchanged. This licenses the same second-gradient
conditional action, not an all-wave-number local Euler PDE.

## Meaning and continuation

Established as stated by the compact moment decomposition and complete
symbol argument. Eight rather than eleven canceled affine moments suffice
for this locality theorem. Retaining physical angular momentum is therefore
compatible with that analytic license; 0070 still has to derive its actual
mean/angle symplectic and kinetic action. Positive gradient coefficients are
the separate construction in 0065. None of these route verdicts by itself
closes the original N1–N7 objective.
