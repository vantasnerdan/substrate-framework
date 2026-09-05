# One exposing strengthening after the first execution

The first execution passed 24 checks. Its all-angular ground-state
identity was fully symbolic in m, but the additional inequality check
used only m=2,3,4. The proof itself already used the exact factor m^2-1.
The supplementary inequality check was strengthened to all integers
m>=2, with m=h+2, h>=0 and gap=h^2+4h+3>0. This is a verifier
coverage repair, not a changed scientific theorem or tolerance.

To reconstruct the first source from final verify.py, replace exactly:

    excess = s.symbols("m_minus_two", integer=True, nonnegative=True)
    angular_gap = s.expand((excess+2)**2-1)
    checks.check("higher-angular correction is strictly positive for every m at least two",
                 angular_gap == excess**2+4*excess+3 and angular_gap.is_positive is True)

with:

    checks.check("higher-angular correction is strictly positive for m at least two",
                 all((value**2-1) > 0 for value in (2, 3, 4)))

No other source change separates the first execution and final source.
