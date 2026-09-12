# Validation and evidence boundary

The proof in `construction.md` derives the survival law, first-event density,
winner probabilities, zero-channel and all-zero boundaries, finite-window
conditioning, and post-reset transition matrix. The importable API exposes the projective
intensities and declared-clock calculation without claiming that Euler creates
either one.

The focused tests include exact symbolic analyzer rotations, arbitrary
positive intensities, a zero-intensity channel, repeated and rotated
measurements after an assumed reset, and invalid normalization/unitarity/rate
mutations. They validate probability algebra only. They do not validate an
Euler analyzer, mixing theorem, event clock, exclusive capture,
open-system/external reset,
action-scale selection, or particle interpretation.

No production numerics or empirical comparator enters. Independent review is
required before the conditional probability theorem is treated as reviewed
campaign evidence.

The repository-interpreter focused run reports `5 passed in 1.59s`, exit
zero. Pinned SHA-256 values are:

- source audit: `3eb37819c0c38767b46f5955a08f531fa19208e989d63cb7266f46f06bcdc970`;
- construction: `d921309b780f44a8fee72ef723b6f11de9648a1806d71c34086edc68432a9897`;
- result: `228cd30e6ad93eef33d647ef73c1a27498980623cb64d04b38a899dccc2478cc`;
- API: `89d459c0121459584a8cab39f0db8e1bec2c20a6ef161c98551e4d51f34f4180`;
- tests: `19d3a4723b8bbfb8e8f4150118a25444b225df53d496a8f1f5c052d7dff00b5d`;
- stdout: `d5579e263d3cf70acb717b91007b4e8f4a22f00eeaa6558ee55446d3ad9d35d1`;
- exit: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
