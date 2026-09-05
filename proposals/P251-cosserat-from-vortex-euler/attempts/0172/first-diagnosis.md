# First execution: reporting API repair

All23 exact mathematical predicates passed. The final attempted
CheckLedger.to_json call does not exist. This is a reporting implementation
failure, not a failed scientific predicate. The output in first.stdout
preserves stdout and traceback (the terminal interleaved the first two
traceback lines ahead of buffered stdout). Replace only that summary call
with the canonical finish method and retain its exit-status semantics.
