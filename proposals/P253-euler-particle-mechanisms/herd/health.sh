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
say "--- freshness (staleness trap probes) ---"
# Probe 1: every worktree attempt dir must have an INDEX row.
for d in proposals/P253-euler-particle-mechanisms/attempts/01[01][0-9]-*/; do
  [ -d "$d" ] || continue
  base="$(basename "$d")"
  if grep -q "$base" "$HERD/INDEX.md" 2>/dev/null; then
    say "index covers $base: ok"
  else
    say "index gap: $base in worktree, no INDEX row (STALE)"; fail=1
  fi
done
# Probe 2: uncommitted attempt paths are invisible to index readers.
untracked="$(git status --short -- proposals/P253-euler-particle-mechanisms/attempts 2>/dev/null | grep '??' || true)"
if [ -n "$untracked" ]; then say "uncommitted attempt paths (warn):"; echo "$untracked" | head -n 10; warn=1; else say "attempts tree committed: ok"; fi
# Probe 3: HEAD's attempt files must be INDEX-covered (catches landed-but-unindexed).
if git rev-parse --verify HEAD >/dev/null 2>&1; then
  headgaps=0
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    case "$f" in proposals/P253-euler-particle-mechanisms/attempts/01[01][0-9]-*/*)
      dir="$(echo "$f" | cut -d/ -f3)"
      grep -q "$dir" "$HERD/INDEX.md" 2>/dev/null || { say "HEAD gap: $f landed, $dir unindexed (STALE)"; headgaps=1; } ;;
    esac
  done <<EOF
$(git show --name-only --pretty=format: HEAD -- proposals/P253-euler-particle-mechanisms/attempts 2>/dev/null)
EOF
  [ "$headgaps" -eq 0 ] && say "HEAD attempt files indexed: ok" || fail=1
fi
say "--- per-agent latest v1 signal (PR-readiness view) ---"
for n in shepherd atlas beacon cipher drift; do
  last="$(grep " $n \[" "$HERD/STATUS.md" 2>/dev/null | tail -n 1 || true)"
  if [ -n "$last" ]; then say "$n: $last"; else say "$n: no v1 signal line (debt)"; warn=1; fi
done
debt="$(grep -cE '^- [^-[]*[^]\[]$' "$HERD/STATUS.md" 2>/dev/null || true)"
say "non-signal STATUS lines (debt): $debt"
if [ "$fail" -eq 0 ]; then
  [ "$warn" -eq 0 ] && say "HEALTHY: v0+v1" || say "HEALTHY: v0 (v1 partial)"
  exit 0
else
  say "BROKEN: v0 needs repair"
  exit 1
fi
