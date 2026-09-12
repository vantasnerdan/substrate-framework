# 0073 shared-worktree collision repair

The content-blind review README was frozen and reported at SHA-256
`43653f748a86d09c9ab0f7c286fa52fc470e36082d5ddb9a9e8e8cf31ec02277`.
A concurrent foundations draft then overwrote that path at SHA-256
`557c33042feacf21b11b5b95fde496cb07665de84c0d8b00f3e2fe79c070071a`.
The constructive draft is preserved at `0074/README.md`; it does not belong to
0073.

The first activation executed while the wrong README occupied 0073 is
preserved as `collision-invalid-activation-schema.*` and has no activation
authority. A concurrent first restoration then raced with the reviewer and
left the README absent; that replay is preserved as
`collision-missing-readme-activation-schema.*`, and the interrupted first
receipt is preserved as `collision-repair-first-race.md`. Neither carries
activation authority.

The reviewer had frozen hash `43653f74...` content-blind before the collision.
Because the wrong constructive README was present when the first schema run
returned zero, the reviewer opened `0066/source-audit.md`, the 0066 supervisor
precision note, `0066/README.md`, `0066/derivation.md`, `0066/result.yaml`, and
`0066/validation.md`, and inspected cited Cao sections under that apparent
activation. It opened no 0072 file and produced no review verdict. The hash
mismatch was detected at the beginning of that review interval, and the
reviewer paused when root reported the collision.

Root recovered the exact original review README from the reviewer's
append-only Codex session patch, verified the frozen SHA-256 above, and only
then replayed the standard `activation-schema.*` quartet. Its exit is exactly
`0`. The pre-exposure review contract and all frozen target hashes therefore
remain independent of the later source inspection. This repair changes no
0066 target file, target hash, review criterion, or scientific claim.
