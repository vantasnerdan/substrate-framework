#!/usr/bin/env bash
# herd/board.sh — render herd/BOARD.md from STATUS v1 lines + INDEX + HEAD.
# NEVER hand-edit BOARD.md: regenerate anytime with this script (idempotent).
# Usage: bash proposals/P253-euler-particle-mechanisms/herd/board.sh
set -u
HERD="proposals/P253-euler-particle-mechanisms/herd"
OUT="$HERD/BOARD.md"
HEAD="$(git log -1 --format=%h 2>/dev/null || echo nogit)"
NOW="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
STATUS="$HERD/STATUS.md"
INDEX="$HERD/INDEX.md"

sig_of()  { echo "$1" | sed -n 's/.*\[\([A-Z]*\)\] \[.*/\1/p'; }
obl_of()  { echo "$1" | sed -n 's/.*\] \[\([^]]*\)\].*/\1/p'; }
fld_of()  { echo "$1" | sed -n "s/.*$2:\\([^ ]*\\).*/\\1/p"; }
msg_of()  { echo "$1" | sed -n 's/.*:: //p' | cut -c1-140; }

{
echo "# Herd board — GENERATED, do not hand-edit"
echo ""
echo "Regenerate: \`bash $HERD/board.sh\` (anyone, anytime; idempotent)."
echo "Generated: $NOW at HEAD \`$HEAD\`. If \`git log -1\` shows a newer commit, re-run — this board predates the branch."
echo "Sources: STATUS v1 \`[SIGNAL] [OBL]\` lines + INDEX + HEAD. Freeform lines are debt, not state."
echo ""
echo "## PR-readiness (latest v1 signal per agent)"
echo ""
echo "| agent | signal | obligation | frontier | blocked-on | one-liner |"
echo "|-------|--------|------------|----------|------------|-----------|"
for n in shepherd atlas beacon cipher drift; do
  last="$(grep " $n \[" "$STATUS" 2>/dev/null | tail -n 1 || true)"
  if [ -z "$last" ]; then
    echo "| $n | (no v1 signal) | - | - | - | freeform STATUS only |"
    continue
  fi
  echo "| $n | $(sig_of "$last") | $(obl_of "$last") | $(fld_of "$last" frontier) | $(fld_of "$last" blocked-on) | $(msg_of "$last") |"
done
echo ""
echo "## Open handoffs (latest line per agent, blocked-on is not -)"
echo ""
found=0
for n in shepherd atlas beacon cipher drift; do
  last="$(grep " $n \[" "$STATUS" 2>/dev/null | tail -n 1 || true)"
  [ -z "$last" ] && continue
  case "$last" in *'blocked-on:-'*) continue;; esac
  echo "- $(echo "$last" | cut -c1-220)"
  found=1
done
[ "$found" -eq 0 ] && echo "- none"
echo ""
echo "## Latest landings (INDEX tail)"
echo ""
tail -n 8 "$INDEX" 2>/dev/null || echo "- INDEX missing"
} > "$OUT"
echo "board rendered: $OUT @ $HEAD"
