# drift firewall review — beacon D-D1 pair table (c44572ed, #107): PASS

Reran dd1_pair_table.py from checkout (HEAD bytes identical to banked
commit, tracked path, check-ignore negative): exit 0, numbers
digit-exact (Q0 −0.01872, climb 0.01592). Freeze c23e0eb7 (09:57:22)
precedes compute c44572ed (09:58:28) ✓.

## (1) Repulsion exactness: genuine import, honestly receipted

F = μb₁b₂/2πd is the textbook Volterra/PK screw result, transcribed
under the frozen "no new formalism" license — NOT field-derived from
the F-A medium, and the paper claims no more than that ("exact PK").
The differentiation check guards transcription, not physics; the
physics content is the sign (like-sign → F>0) on an analytic form
with ZERO fitted constants. No fit knobs exist anywhere in the
script. Tier for the fold: imported-analytic + regression-guarded,
not derived-from-medium. Genuine as labeled.

## (2) Q0 tolerance: honest guard, reported-not-gated

|0.01872−0.01881|/0.01881 = 0.46% ≈ 0.5% claimed ✓ arithmetic exact.
The 25% bar is the gate (bug-catcher); the 0.5% is reported headroom,
never promoted to a precision claim. Single-point comparison stays a
guard — the receipts don't round it up. Honest. (The 0.5% is partly
fortuitous given B1's own ~15% analytic gap; as a guard with a 25%
bar this is fine, as metrology it would not be — and it isn't.)

## (3) Exclusions honored, order held

Screw+tilt null: no tilt/J(eps) terms anywhere in D-D1 ✓ (direct
elastic only, as frozen). Fluctuation mediation absent, explicitly
out per P3-C ✓. Bare bound stays dead: Q2 monotonic-repulsion →
EMPTY, the pre-registered kill honored with no bound-pair reading
attempted ✓. Edge fence intact: edge forms banked as D-D2 INPUTS,
no ALIVE language on edges ✓. D-D2/D-D3 not run, stated in receipts
✓ — ordered rounds held.

Slips disclosed in-tree (Q2 reword, NU NameError, orphan tail) are
drafting-level with no physics iteration — correctly banked as slips,
not silently fixed. Verdict: PASS. D-D2 may proceed on charter.
