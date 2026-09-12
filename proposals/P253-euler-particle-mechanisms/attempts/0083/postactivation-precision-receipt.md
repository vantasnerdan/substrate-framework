# Postactivation convention-precision receipt

The central activation completed before the supervisor's bounded convention
audit arrived.  No 0083 scientific body had been opened.

## Chronology and immutable inputs

- originally activated README SHA-256:
  `327ec0ed473129a01059df0f7da47863d57f98c4d583a440caecbb747bd94683`;
- activation command SHA-256:
  `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`;
- activation stdout SHA-256:
  `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`;
- activation stderr SHA-256 (empty file):
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
- activation exit SHA-256:
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`;
- activation exit content: exactly `0`;
- activation stdout: `WORKFLOW VALID`.

## Bounded repair

The original README conflated the cylindrical base measure with the
vorticity-weighted measure and mixed two scale-coefficient conventions.  The
corrected README now freezes

    dnu=r dr dz,
    dm_epsilon=zeta_epsilon dnu,
    kappa=integral dm_epsilon,
    Rbar=kappa^(-1) integral r dm_epsilon,

and

    s_epsilon=C_s(kappa,R) epsilon [1+eta(epsilon)],
    delta=(C_s/R) epsilon [1+eta(epsilon)],
    delta'=(C_s/R)[1+o(1)].

No route, selection criterion, or scientific body changed.  Corrected README
SHA-256:
`12e5f0cd72aa392d2959d2608be98fd49a3f3bda699bf885271ece4de32ce10d`.
The corrected README awaits central schema replay before bodies open.

## Corrected-contract schema replay

Root centrally replayed the schema after receiving corrected frozen README
SHA-256
`12e5f0cd72aa392d2959d2608be98fd49a3f3bda699bf885271ece4de32ce10d`.
The replay completed with exit content exactly `0`, stdout `WORKFLOW VALID`,
and empty stderr.  Its receipt hashes are:

- command:
  `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`;
- stdout:
  `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`;
- stderr:
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
- exit:
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The original activation receipts were preserved under
`preprecision-activation-schema.*`; they have the same byte hashes because
the repository-schema command and repository state output were unchanged.
The README status was then changed from “awaiting replay” to “replayed”; that
provenance-only status version has SHA-256
`17ad166ba1bdeed34432032239c198c232ae16118d229510da5f6dd283919d95`.
No frozen route, equation, criterion, or scientific conclusion changed in the
status edit.

## Completed-review authority refresh

Before author freeze, the README dependency prose was refreshed from the
historical “P253/0082 pending” state to the completed independent-review
boundary.  The pinned review hashes are `0082/review.md`
`75e55672e359fdf91b91017d6098e266eb605f4ac7e289132153e22a4f7cf334`
and `0082/verdicts.yaml`
`69eee8688f986e7ee0463d92d9a115f6e397a74f2bbc13370d8d0ba414e175a4`.
This is an authority/provenance refresh only: 0082 still leaves the continuous
fixed-mean Cao path, physical KKS normalization, and controls open, exactly as
the activated candidate contract required.

After that provenance-only refresh, the final author-freeze README SHA-256 is
`30d5ef578b3474b9af8afcacb059e8e08c2c88d55673728366256f7c46aad56a`.
