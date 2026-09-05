# Exact input-domain repair

The first10 direct scientific/API tests passed. An additional explicit
NaN domain probe then showed that SymPy reports neither is_real nor
is_finite as False for NaN; the constructor admitted invalid data.
domain-probe.stdout preserves that failing test. The repair explicitly
rejects non-expression, NaN and infinite amplitude inputs. All field,
moment, energy and current formulas stay unchanged. The directly changed
API/test pair is replayed in repaired-pytest.stdout, without a full suite.
