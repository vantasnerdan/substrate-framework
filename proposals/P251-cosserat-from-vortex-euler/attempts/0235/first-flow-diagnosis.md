# Exactness repair of the independent flow oracle

The first flow-trace execution failed its first identity. Its source hash
was `46e66b6c8ae84339008b9c5cd03f7028017886214a722796422e70f5d6259a0d`;
first-flow.stdout preserves the traceback, hash
`460cb900685ec446a07fc5c9a1cbca90e7e06cbc7f221942d473f2e65faa77f8`.

The Gaussian moment routine used Python `/2` on even integer powers,
introducing Float exponents such as width**1.0 into an intended exact
symbolic polynomial. The correction uses integer `//2`, justified by
the already enforced even-power branch. It also prints the directly
computed trace and its unsimplified scientific comparison before the
predicate. No material-flow equation, moment, or predicted correlation
identity changed. The exact repaired output has zero defining defect
and all three checks pass. Both executions remain recorded.
