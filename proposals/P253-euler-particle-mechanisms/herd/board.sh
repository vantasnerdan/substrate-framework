#!/usr/bin/env bash
# herd/board.sh — render herd/BOARD.md from STATUS v1 lines + INDEX + HEAD.
# NEVER hand-edit BOARD.md: regenerate anytime with this script (idempotent).
# Usage: bash proposals/P253-euler-particle-mechanisms/herd/board.sh
set -u
export LC_ALL=C.UTF-8
HERD="proposals/P253-euler-particle-mechanisms/herd"
OUT="$HERD/BOARD.md"
HEAD="$(git log -1 --format=%h 2>/dev/null || echo nogit)"
NOW="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
STATUS="$HERD/STATUS.md"
INDEX="$HERD/INDEX.md"

sig_of()  { echo "$1" | sed -n 's/.*\[\([A-Z]*\)\] \[.*/\1/p'; }
obl_of()  { echo "$1" | sed -n 's/.*\] \[\([^]]*\)\].*/\1/p'; }
fld_of()  { echo "$1" | sed -n "s/.*$2:\\([^ ]*\\).*/\\1/p"; }
kind_of() { k="$(fld_of "$1" bkind)"; case "$k" in physics|ack|decision) echo "$k";; *) echo "ack";; esac; }
msg_of()  { echo "$1" | sed -n 's/.*:: //p' | cut -c1-140; }
age_min_of() { ts="$(echo "$1" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}Z' | head -n 1)"; [ -z "$ts" ] && { echo "?"; return 0; }; age=$(( ($(date -u +%s) - $(date -u -d "$ts" +%s)) / 60 )); [ "$age" -lt 0 ] && age=0; echo "$age"; }

{
echo "# Herd board — GENERATED, do not hand-edit"
echo ""
echo "Regenerate: \`bash $HERD/board.sh\` (anyone, anytime; idempotent)."
echo "Generated: $NOW at HEAD \`$HEAD\`. If \`git log -1\` shows a newer commit, re-run — this board predates the branch."
echo "Sources: STATUS v1 \`[SIGNAL] [OBL]\` lines + INDEX + HEAD. Freeform lines are debt, not state."
echo ""
echo "## TL;DR (10-second scan)"
echo ""
for n in shepherd atlas beacon cipher drift; do
  last="$(grep " $n \[" "$STATUS" 2>/dev/null | tail -n 1 || true)"
  if [ -z "$last" ]; then echo "- $n: no v1 signal"; continue; fi
  blo="$(fld_of "$last" blocked-on)"
  if [ "$blo" = "-" ] || [ -z "$blo" ]; then echo "- $n: $(sig_of "$last") [$(obl_of "$last")], clear"; else echo "- $n: $(sig_of "$last") [$(obl_of "$last")], waiting on $blo [$(kind_of "$last")] ($(age_min_of "$last")m)"; fi
done
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
attn=""
for n in shepherd atlas beacon cipher drift; do
  hit="$(grep -n " $n \[" "$STATUS" 2>/dev/null | tail -n 1 || true)"
  [ -z "$hit" ] && continue
  waitln="${hit%%:*}"; last="${hit#*:}"
  case "$last" in *'blocked-on:-'*) continue;; esac
  blo="$(fld_of "$last" blocked-on)"
  kind="$(kind_of "$last")"
  ack="dep"
  if [ "$kind" = "ack" ]; then
    ackhit="$(grep -nF "ack:$blo" "$STATUS" 2>/dev/null | awk -F: -v s="$waitln" '$1>s' | tail -n 1 || true)"
    if [ -n "$ackhit" ]; then
      ackts="$(echo "${ackhit#*:}" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}Z' | head -n 1)"
      ack="ACKED${ackts:+ @$ackts}"
    else
      ack="no-ack"
    fi
  fi
  [ "$kind" = "decision" ] && ack="owner-decision"
  wage="$(age_min_of "$last")"
  echo "- [waiting ${wage}m, $kind, $ack] $(echo "$last" | cut -c1-180)"
  [ "$wage" = "?" ] && attn="${attn}- UNDATED: $n re-post this wait with a full YYYY-MM-DDTHH:MMZ timestamp (§8), then it can age
"
  if [ "$kind" = "ack" ] && [ "$ack" = "no-ack" ]; then
    esc=""; [ "$wage" != "?" ] && [ "$wage" -gt 60 ] && esc=" — STALE past 60m, ESCALATE-TO-SHEPHERD"
    attn="${attn}- UNACKED: whoever starts on $n's block, post a STATUS line containing \`ack:$blo\` (§9)$esc
"
  fi
  found=1
done
[ "$found" -eq 0 ] && echo "- none"
echo ""
echo "## Needs attention (do these, oldest strain first)"
echo ""
if [ -z "$attn" ]; then echo "- nothing outstanding"; else printf '%s' "$attn"; fi
echo ""
echo "## Latest landings (INDEX tail)"
echo ""
tail -n 8 "$INDEX" 2>/dev/null || echo "- INDEX missing"
echo ""
echo "## Validation receipts (INDEX rows carrying vrfy:cmd:scope:exitN)"
echo ""
vrows="$(grep 'vrfy:' "$INDEX" 2>/dev/null | tail -n 6 || true)"
if [ -z "$vrows" ]; then echo "- none yet — add \`vrfy:<cmd>:<scope>:exit<N>\` to your landing INDEX row (§11)"; else echo "$vrows"; fi
echo ""
echo "## Gap closure (herd/GAPS.md: gap → next artifact → owner/class)"
echo ""
grows="$(grep '^|' "$HERD/GAPS.md" 2>/dev/null | grep -v '^|---' || true)"
if [ -z "$grows" ]; then echo "- GAPS.md missing"; else echo "$grows"; fi
} | sed 's/[[:space:]]*$//' > "$OUT"
echo "board rendered: $OUT @ $HEAD"
