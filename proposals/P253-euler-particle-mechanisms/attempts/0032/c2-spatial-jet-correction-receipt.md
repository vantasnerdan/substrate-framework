# C2 spatial-jet correction provenance

This append-only receipt pins the bounded 0038 correction without replacing
the original successful source-recurrence evidence.

| Object | SHA-256 |
| --- | --- |
| pre-correction `verify_c2_source_recurrence.py` | `194fa4b83f08689bd77dae3cc0da6f0a26b12ffe72dd570e716367754a5d95c8` |
| pre-correction `c2-first-success.stdout.txt` | `6fb56ad2baf1f7ebdbf2f67e88a000e80bca753c11c0a7859edfd0949c71e878` |
| post-correction `verify_c2_source_recurrence.py` | `74535827251d9cb3edc23e2a9c016984f7cc0e7b2df70f16fecb0d16d218ddc2` |
| first corrected `c2-corrected-spatial-jet.command.txt` | `f3829407139d2e944e2a7e936478fb067812e7b7b01c82999cf2baf536c76d9c` |
| first corrected `c2-corrected-spatial-jet.stdout.txt` | `712618c03a43aa2d931608120af61e7b184a9b9a7102ad604e7fc61a89e0a247` |
| first corrected `c2-corrected-spatial-jet.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |
| repository-interpreter `c2-repository-interpreter.command.txt` | `7b47e0c64dd8a4b5d54b2fbbada624c65508c54dc628fac42947ceb0635112ad` |
| repository-interpreter `c2-repository-interpreter.stdout.txt` | `712618c03a43aa2d931608120af61e7b184a9b9a7102ad604e7fc61a89e0a247` |
| repository-interpreter `c2-repository-interpreter.stderr.txt` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| repository-interpreter `c2-repository-interpreter.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

The repository-interpreter command was

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0032/verify_c2_source_recurrence.py

It exited `0`, emitted no stderr, and reproduced the corrected stdout byte for
byte.  The new eighth assertion checks

    C2_new[1,0]-C2_spatial_degree_at_most_2[1,0]
      =sqrt(2)*(3*cos(q)^4/2-13*cos(q)^2/8-9/16).

The restored entry is outside the single-`C2` trace functional.  The exact
`9*pi^2+2*pi^2` trace coefficient and `88*pi^2` discriminant coefficient are
therefore unchanged.
