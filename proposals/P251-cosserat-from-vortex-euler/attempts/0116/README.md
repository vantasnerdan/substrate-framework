# 0116 — periodic Beltrami coherence and actual Fourier transfer

Owner `/root/construction_review`; this directory only. Parent authorizes
a materially different candidate for the STRONGER actual full-Euler
finite-wave-number realization. Accepted conditional results through
v0.175.0 are inputs and are not reopened. This is neither promotion work
nor a second audit of the old Cauchy--Born construction.

Frozen question: can identical periodic cells carrying actual knotted
Beltrami tubes prevent local-frequency dephasing and provide nonzero
COMPLETE-fluid coarse Fourier translation from an optical angular packet,
with an explicit scale-controlled physical-observable remainder?

Candidates: A primary periodic inverse-localization theorem plus exact
Bloch fiber and complete packet-moment coupling; B a weak intercell
pressure/return interaction calculated from the exact periodic Leray
operator. A packet's subparcel angular moment alone is not the complete
fluid Fourier momentum. Selection is exact periodic background existence,
same-cell coherence, physically defined Fourier transfer, and a controlled
linear/full-Euler remainder. Neither an arbitrary static CB profile nor a
positive local fixed-frame action is treated as that stronger object.

Imports to read: primary Enciso/Peralta-Salas/Torres de Lizaur periodic
Beltrami inverse-localization results;0112/0114 only for their established
actual elliptic-packet scope, without duplicating that action calculation.
No empirical comparator or numerical frequency-splitting selection. Use
exact Bloch/projection identities and analytic bounds first; any later
soft-splitting numerical design would load the small-ratio prescriptions
before its representation is frozen.

## Frozen result and receipts

Established bounded construction: exact same-core periodic background;
identical-cell coherence; complete-fluid Bloch stress transfer; compact
Kelvin-generator preparation with exactly zero initial Fourier mean and
strictly nonzero stress; finite-time actual Euler transfer; and a small
transported-seed attachment preserving the positive physical-frame
angular action with controlled core-observation error. Full proof and
scale order: `periodic-bloch-transfer.md`. Primary imports and digests:
`source-receipt.md`.

The autonomous optical-band candidate has one named missing construction:
a Bloch spectral/resolvent residue connecting the physical angular packet
to the measured optical-frequency coarse stress. This is distinct from
the established finite-time response and from 0117's actual hybrid
centroid-plus-ambient observable. No unrestricted constitutive closure,
spectral exhaustion, or parent completion is inferred.

Executed oracle: `verify.py`, 18/18 exact checks, exit 0 in
`corrected-run.txt`. The initial missing-PYTHONPATH invocation is preserved
in `first-run.txt`; it ran no scientific checks. Exact fixture outputs:
positive Gram residual `mu=1`; complete transverse acceleration
`-7*I/50` at rational Bloch wave number `1/7`, after exact zero-mean
Kelvin preparation. The fixture checks algebra and projection, not
knotted-tube existence; the source theorem and analytic compact Gram
argument supply that part of the construction.

`ruff check` passes. PDF and text digests were replayed and match the
source receipt. New-file whitespace check reports no whitespace errors
(the no-index diff exit 1 denotes a new file, not a check failure).
All work is confined to this attempt directory; existing modules and
accepted records are unchanged by this child. Frozen for the parent's
bounded integration review; any optical-residue construction belongs to
a new attempt.
