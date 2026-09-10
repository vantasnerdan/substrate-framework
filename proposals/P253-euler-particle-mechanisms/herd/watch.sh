#!/usr/bin/env bash
# herd/watch.sh — digest notifier (default: atlas self-watch).
# Watches: new commits on the branch (landings), excluding TARGET's own.
# One digest prompt per interval max. Uncommitted edits NOT watched (v1 limit).
# Usage: bash proposals/P253-euler-particle-mechanisms/herd/watch.sh [interval_s] [target] [statefile]
set -u
export LC_ALL=C.UTF-8
INTERVAL="${1:-60}"
TARGET="${2:-atlas}"
STATE="${3:-/tmp/${TARGET}-watch.head}"
[ -f "$STATE" ] && SEEN="$(cat "$STATE" 2>/dev/null || true)"
if [ -z "$SEEN" ]; then SEEN="$(git log -1 --format=%H 2>/dev/null || true)"; echo "$SEEN" > "$STATE"; fi
echo "watch armed @ ${SEEN:0:8}, every ${INTERVAL}s (peer commits only)"
while true; do
  sleep "$INTERVAL"
  HEAD="$(git log -1 --format=%H 2>/dev/null || true)"
  [ -z "$HEAD" ] && continue
  if [ "$HEAD" != "$(cat "$STATE" 2>/dev/null || true)" ]; then
    RANGE="$(git log --format='%h %s' "$SEEN..$HEAD" 2>/dev/null | grep -vE "^[0-9a-f]+ $TARGET" | head -n 8 || true)"
    if [ -z "$RANGE" ]; then RANGE="$(git log -1 --format='%h %s' "$HEAD")"; fi
    echo "$HEAD" > "$STATE"; SEEN="$HEAD"
    N="$(echo "$RANGE" | wc -l)"
    herdr agent prompt "$TARGET" "[WATCH] $N peer landing(s): $(echo "$RANGE" | tr '\n' ';' | cut -c1-380)" 2>&1 | head -n 2
  fi
done
