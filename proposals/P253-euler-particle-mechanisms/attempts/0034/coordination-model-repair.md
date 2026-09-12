# Operational model correction during the live campaign

Root inspected actual Codex turn_context records, not only Herdr process argv.
The particle-foundations session01a07656-64cf-74e2-ac94-9ba86652daef began as
gpt-5.6-sol/high at2026-09-06T10:50:13.607Z, then recorded
gpt-5.6-luna/low from2026-09-06T13:20:24.656Z through13:47:00.603Z.
Its launch argv still named Sol High, so the earlier process-only check did
not establish the model used for those later turns. The cause of the switch
was not established. The user's Sol High scientific-work instruction remains
controlling.

Root stopped only that worker, preserved its files/session, and resumed the
same session with explicit `--dangerously-bypass-approvals-and-sandbox
-m gpt-5.6-sol -c model_reasoning_effort="high"`. The actual next turn_context
at2026-09-06T13:49:29.271Z recordsgpt-5.6-sol/high, approval_policynever.
The other two active scientific workers' current turn_context records remain
Sol High. No unrelated pane or user worktree was changed.

The affected0032 additions remain unaccepted campaign evidence. The resumed
worker is auditing their actual source-derived C0/C1 coefficients, singular
resonance scaling, missing-mu correction, claimed error bound and cubic
recurrence before consuming those results. Model identity alone neither
proves nor refutes a mathematical claim; concrete source-equation checks do.
In particular, printing an unsolved recurrence is not the requested analytic
coefficient construction. The resumed task explicitly continues that
construction in-run.
