#!/usr/bin/env bash
# herd/health.sh — v0+v1 comms health. Exit 0 = healthy, 1 = broken v0.
# Usage: bash proposals/P253-euler-particle-mechanisms/herd/health.sh
set -u
HERD="proposals/P253-euler-particle-mechanisms/herd"
fail=0; warn=0
say() { printf '%s\n' "$*"; }

say "== herd health $(date -u +%Y-%m-%dT%H:%M:%SZ) =="
for n in shepherd atlas beacon cipher drift; do
  if [ -f "$HERD/inbox/$n.md" ]; then say "inbox $n: ok"; else say "inbox $n: MISSING"; fail=1; fi
done
[ -f "$HERD/STATUS.md" ] && say "STATUS.md: ok" || { say "STATUS.md: MISSING"; fail=1; }
[ -f "$HERD/README.md" ] && say "README.md: ok" || { say "README.md: MISSING"; fail=1; }
for f in protocol-v1.md INDEX.md health.sh; do
  [ -f "$HERD/$f" ] && say "v1 $f: ok" || { say "v1 $f: missing (warn)"; warn=1; }
done
if command -v herdr >/dev/null 2>&1; then
  list="$(herdr agent list 2>&1)"
  for n in shepherd atlas beacon cipher drift; do
    echo "$list" | grep -q "\"name\":\"$n\"" && say "agent $n: present" || { say "agent $n: ABSENT"; fail=1; }
  done
else
  say "herdr CLI: not found"; fail=1
fi
blocked="$(grep -c '\[BLOCKED\]' "$HERD/STATUS.md" 2>/dev/null || true)"
say "STATUS [BLOCKED] lines: $blocked"
if [ "$fail" -eq 0 ]; then
  [ "$warn" -eq 0 ] && say "HEALTHY: v0+v1" || say "HEALTHY: v0 (v1 partial)"
  exit 0
else
  say "BROKEN: v0 needs repair"
  exit 1
fi
