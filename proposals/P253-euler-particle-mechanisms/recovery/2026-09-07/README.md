# Temporary-file harvest at the user-requested pause

The coordinator and both workers shared `/tmp/substrate-particles`. Its entire
uncommitted P253 surface was harvested into this checkpoint. The three P253
detached validation worktrees (`substrate-particles-validate-bd42096`,
`substrate-validate-6269e87`, `substrate-validate-efd222b`) were all clean.
Ignored files under campaign/source/test/memory paths were generated Python
bytecode, not additional source or proof drafts.

`tmp-scripts/` preserves 46 historical author scripts and coordination text
from `/tmp`; `tmp-receipts/` preserves 39 scratch check/validation outputs.
Three GitNexus stdout copies omit terminal blank lines for repository
whitespace hygiene; the durable local backup preserves their exact bytes.
They have `.archived.txt` suffixes because they are historical evidence, not
current entry points. Some scripts overwrite old receipts or address obsolete
Herdr panes: inspect them and adapt to a new append-only destination before
any future execution. No archived script was executed during harvesting.

A separate durable local copy at `/home/dan/p253-recovery/2026-09-07/` preserves
the identified P253 temporary files and primary-source caches (about 41 MB),
including raw earlier transcripts and third-party PDFs/source extracts.
Those raw caches/transcripts are not published as scientific source files.
`local-backup.sha256` records their content hashes and durable paths. The
canonical source URLs/hashes remain in the campaign source audits. All author
work needed to resume 0107 is committed; the local backup supplements it.

p4 had not yet created the Hill proof body when stopped. Its reported
calculation and remaining proof obligations are transcribed in
`../../attempts/0107/pause-state.md` (the direct path is
`proposals/P253-euler-particle-mechanisms/attempts/0107/pause-state.md`).
The campaign-level `PAUSED.md` is the resume entry point. No agents were
restarted and no mathematical route was continued during this harvest.
