# Independent review: reusable Euler-orbit reduction

Reviewer: `/root/smooth_core_review`, distinct from the author of 0049 and
the API. Date: 2026-09-05. One bounded scientific/code review under AGENTS.md
and the physics skill; no earlier microscopic claim is reopened.

## Decision

Established as stated: `euler_orbit.py` implements the exact conditional
two-momentum reduction and the declared physical affine-cage angle map.
The API retains the full mixed Hessian, the gyroscopic connection and the
singular-map boundary. Its paired-ensemble tests independently vary all
reaction momenta before averaging. No load-bearing correction is requested.

These are reusable algebraic definitions, not a construction of the Euler
orbit or a proof that the constrained model is invariant under the full
Euler PDE. Symbolic positivity, reality, finiteness and nonzero conditions
that cannot be decided remain caller hypotheses. This conditional status
is explicit in the API and is appropriate for these inputs.

## Exact reduction and independent action checks

I read the attempt README, the complete module and four tests, and the
first, failed finite-input and repaired finite-input receipts. The repaired
receipt passes all four tests. The historical diagnostic mismatch is
preserved, not hidden. No redundant replay was required for this review.

For momentum p=(r,s), configuration q and velocity V=(Bdot,qdot), write
the complete action as

```
L=p^T D V-p^T P p/2-q g^T p-h q²/2.
```

Variation gives `p=P^-1(DV-gq)`. Substitution gives exactly
`M=D P^-1 D`, gyro `-D P^-1 g`, and
`K=h-g^T P^-1 g`. This is what the module computes from the (r,q,s)
input ordering. Positive full H implies P>0 and K>0 by its Schur
complement; nonzero pairings make D invertible and hence M>0.

The first test reconstructs the action after solving its momentum
Euler--Lagrange equations with a genuinely mixed positive matrix. It
therefore checks the source ordering and gyro sign against an independent
variation, rather than repeating the implementation formula. The q qdot
total derivative remains in the returned gyro, as documented.

For the paired ensemble, reversing both pairings changes only the sign
of gyro. The second test varies four distinct momenta in the unreduced
averaged action. Its positive time-even action agrees with the reduction.
The tied-momentum mutation explicitly loses every velocity-dependent
term. Thus the physical independent-reaction premise is exposed rather
than silently replaced by averaging symplectic forms first.

## Physical map and input boundary

Let `c=M_BB+M_Bq`, `d=M_BB+2 M_Bq+M_qq`, and beta=B-q.
Then the original kinetic form in (beta,q) has cross coefficient c and
q coefficient d. The map `Psi=beta+(d/c)q` eliminates that cross term,
giving `J_Psi=c²/d`, `J_beta=det(M)/d`, and restoring coefficient
`K_Psi=K c²/d²`. These are precisely the returned values and are positive
under the declared assumptions. Since Psi changes by the same constant
as beta under a common frame rotation, it is an absolute angle variable.

Positive M does not imply c is nonzero. The API correctly detects known
c=0, and its test supplies an explicitly positive matrix with that
singularity. Unknown symbolic c remains the documented nonzero premise;
the function does not invent a different physical cage to remove it.
The map diagonalizes only the time-even kinetic part. The docstring
properly requires the gyroscopic term to be separately retained or
canceled by the paired ensemble.

Known nonsymmetry, nonpositivity, zero pairings and nonfinite inputs are
rejected. The repaired validation checks finiteness before subtracting
matrices, avoiding the earlier infinity-induced misleading symmetry
diagnostic. Symbolic signs use exact leading principal minors and are not
guessed from numerical samples. This proof has no small-ratio numerical
remainder; exact action variation is the claim-appropriate oracle.

## Frozen hashes and disposition

- `src/substrate_framework/euler_orbit.py`: `98b41a31f0076d42213e94a28103021a96f70519f3c4a326eff172647f4597cb`
- `tests/test_euler_orbit.py`: `069c34d2975f720e35fc7a11896245635556e81a793b1e1b92a712c739c5d726`
- `finite-input-repaired.stdout.txt`: `24d9695b084f0b4f86323ff6fc6b47bb86cc59912b520d06328583c3208b8ad2`

Acceptance as conditional importable infrastructure is recommended.
This attachment neither promotes a microscopic existence claim nor
certifies completion of the parent continuum construction.
