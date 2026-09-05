# Bounded integration correction: signed Kelvin reconstruction

The parent identified a concrete convention mismatch: the initial API
returned the positive stationary circulation-constraint offset, while the
campaign's reconstruction generator A has its NEGATIVE sign. The norm
identity H-K=rho Gram(A) could not distinguish these signs. The original
implementation and pytest receipts above remain historical evidence, not
the final signed reconstruction convention.

`material_kelvin_operator` now returns

    A xi=-P[(u.grad)xi+(D xi)^T u]
        =P(xi cross omega)-curl(xi cross u),
    xi_t=A xi

for fixed Kelvin data. The earlier0100 README equation `A=curl F-v`
describes the original offset and is superseded by this correction.

The existing nonzero pointwise Fourier assertion is reversed consistently.
An independent signed assertion now also checks the mean of A on the
explicit phase fixture: it is lambda*(2/3)e_x, since its induced velocity
has that mean and a periodic curl has zero mean. This check uses positive
AND negative lambda and fails for the old API sign, whereas its squared
norm does not. All pressure/transport and H-K identities remain intact.

Density validation additionally rejects explicit complex and nonfinite
values in BOTH coadjoint and material APIs. Undetermined symbolic density
remains allowed under the documented positive-real-finite hypothesis.

One bounded correction replay: `python -m pytest tests/test_euler_fourier.py
-q` returns16 passed in3.54s, exit0; captured in `correction-pytest.txt`.
Ruff and scoped diff checks pass. No unrelated scientific replay or graph
reindex was run. Prior logs are preserved unchanged.

Corrected SHA256:

- src/substrate_framework/euler_fourier.py:
  4e00da99e211f340e729e362c59c7c524eba24b9cae98d5dfbd266a67f4e94fb
- tests/test_euler_fourier.py:
  2aab041922e7d7983e20ea6898fbaf9bd83c63d789425d251cdb5c29ba1519f3

Correction verdict: established at the requested signed-reconstruction and
density-domain boundary. No parent scientific claim is promoted here.
