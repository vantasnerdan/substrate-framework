# 0080 — invariant EPS material tubes with a continuous ambient phase

Owner `/root/orientation_construction`; this directory only. Parent P251 /
issue #198 remains active with its original full smooth-Euler physical
micropolar objective and permitted slow-affine constrained closure.

Frozen continuation from 0075: replace a time-dependent arbitrary cubic
material partition by the actual invariant EPS solid-torus domains as
finite material parcels. Their complement is a continuous ambient material
phase, not an infinite parcel with a fictitious centroid. Both phases are
the same Euler fluid; density, pressure matching, energy and momentum
remain those of the single material action.

Candidate A: exact tube-centroid cotangent reduction plus ambient material
coordinates, with a divergence-free affine macro lift and all induced
ambient motion retained in its full Schur operator. Candidate B: retain
additional phase-relative translation and shape coordinates before
elimination if their omission obstructs the physical spin connection.

Selection is by actual material observables, complete momentum and kinetic
bookkeeping, stationary reference interfaces, nonzero physical spin and
exact Euler/Kelvin admissibility. A frame convention is not assigned rotor
mass; an ambient Schur correction is not assumed to vanish. No empirical
comparator or numerical soft-mode gate is used. The first achievement is
the exact cotangent/kinetic decomposition and its physical spin row.

Earned scope: `tube-moment-construction.md` proves stationary physical
tags and the actual nineteen-row tube/global response construction,
including its nondegenerate trefoil/ABC prototype, disjoint dual fields,
exact physical tube spin, and explicit affine lift/current identity.
It retains all exterior velocity and boundary terms. The material tag
transport/reconstruction is a distinct obligation, not inferred from
reference tangent matching; independent 0082/0084 work addresses it.

The first exact run of `verify.py` passed 17/17 checks. Its unchanged
output is `first.stdout`. No numerical existence or soft eigenvalue
claim is made by this narrow algebraic verifier. This attempt is ready
for independent review; it is not a parent-completion declaration.
