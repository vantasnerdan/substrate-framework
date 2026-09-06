# 0069 bounded correction receipt

Date: 2026-09-06

This is the single correction package requested by independent review0069.
It preserves every compact-vorticity, pressure, flux, and scalar-to-transverse
formula and changes their claim boundary as follows:

1. The `r^-3` impulse, `d^-3` kinetic cross term, and `r^-3` pressure
   quadrupole are the first permitted coefficients. They lead when their
   displayed contractions are nonzero; cancellation gives faster decay.
2. The old broad direct-smooth-localized headline is restricted to the proved
   compact-vorticity/weighted-moment/radial-flux and local-scalar transverse
   classes.
3. The exact smooth finite-energy escape field
   `u_a=(a cross x)/(1+|x|^2)^(3/2)` is recorded. It has an oriented `r^-2`
   tail and anisotropic `d^-1` translated cross energy, while supplying no
   scalar charge, Gauss law, current, or persistent carrier.

No implementation predicate or prior oracle changed, so neither 0064 nor 0065
was rerun merely to create a new tally.

## Pre/post hashes

| Artifact | Pre-correction SHA-256 | Post-correction SHA-256 |
| --- | --- | --- |
| `0064/README.md` | `84b86e0eae440081f7db84c36e6e6db799e091a8645c1b78c2032cb8aaad4808` | `3ea0f3b8f0b1697032f240b36f48e94b3f278a0160780e616e1c841199e34e7d` |
| `0064/derivation.md` | `39afd045652ca89b3e888e5394857f9d54969a4493b81eb69d781e4be0ff9fd7` | `4e2babcb76e06a0b9d53e67697dd3c861ddeadf36fabbf88ab23fb82a3386add` |
| `0064/result.yaml` | `b205029a3e1d28e8cc646c3b22906f197ceaccc03cf444bc4c4d750100082fb6` | `a118834a38c970581531d54b500866aa4bcf1fc3a7294f107b91068d69cb9232` |
| `0065/derivation.md` | `2b1a34608665741b41a0add3d0a941ce68d5b0df8bee8a0f14a886a1d97555ab` | `71d137d5aa14737609d3484781455623e937f59c4ef4a0c61dd0b1ea88fbeca0` |
| `0065/result.yaml` | `d86bd7589121f82d67b24635764bffef906f756ecd1d8b4851d4f8d130703807` | `9740bb5b86c17bcbaaeddf94a144e084ee2b8e8181e329b99138abd15194ff56` |
