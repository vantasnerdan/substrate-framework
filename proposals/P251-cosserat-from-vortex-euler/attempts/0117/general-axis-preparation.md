# Extension of the actual isotropic preparation to a general mark axis

The application of material-coarse-response.md does not need a transverse
spin endpoint.0114 proves a nonzero full three-component endpoint map
D from the two-dimensional actual Euler phase plane; it does not assert
that its transverse part is nonzero. Keep that exact source statement.

Let A be the physical core-rotation observation, a rank-two 3-by-2 matrix
with image in the normal plane to u. The set of vectors n with n^T A=0
is one proper subspace, and n^T D=0 is another. Their union is not R^3.
Choose a unit n outside both. It may have both axial and transverse
components relative to the local core; it is a declared observation mark,
not a claim that the physical tilt has gained a longitudinal component.

Put a=n^T A and d=n^T D. Both are nonzero row vectors. There is z0 with
a z0=1 and d z0!=0. To construct it, choose z1=a^T/(a a^T). If d z1!=0
use z1. Otherwise choose any z2 in ker(a) with d z2!=0 and set z0=z1+z2.
Such z2 exists in the remaining case: if d annihilated ker(a), it would
be proportional to a, contradicting d z1=0 and d!=0.

The covariant preparation z0(n.q), followed by Haar rotation of its FULL
state, gives exactly

    average_R R(A z0)n^T R^T q = q/3,
    average_R R(D z0)n^T R^T q = (d z0)q/3 !=0.

These follow from the invariant-tensor identity average(R F R^T)=tr(F)I/3.
Thus0114's actual nonconstant spin supplies the nonzero coherent endpoint
response used in0117 without strengthening its lemma or assuming a
canonical-spin identification. This is a constructive extension of the
earlier transverse sufficient condition, not a change to that condition's
proof or to an accepted claim. Ensemble/Bloch response regularity and the
remaining constitutive identification are still separate obligations.
