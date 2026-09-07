# First exact-oracle failure

The first execution exited `1` at the Maxwell momentum identity.  This was an
exposing implementation omission, not a physical sign reversal: the script
used arbitrary symbolic `B` but silently tested the constrained identity as if
`div B=0` had already been imposed.  For arbitrary `B`, the exact identity is

    partial_t(epsilon_EM E cross B)
      =div sigma_EM-f-(B/mu_EM)div B.

The repaired oracle must first expose this magnetic-monopole defect and then
parameterize `B=curl A` before asserting the constrained momentum law.  The
failed command, empty stdout, traceback, and exit value are preserved beside
this note.  It has no scientific route verdict.
