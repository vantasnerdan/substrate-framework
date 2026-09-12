# Tool receipts - 0144 (beacon)

- `run_m1.py`: exit 0 (31.0 s). 6 rows printed (d/base/dlike/
  dopp/noise); fit A=-1.3609e-3, like-res 0.545, opp-res 0.899,
  significance 0.1x → KILL (a).
- A3 run_a3.py (ring_state/shoot backbone): READ ONLY. Custom
  per-ring-aa rhs + Neumann integral local to run_m1.py.
- Bug trail: window-block splice + ds-line drop (both caught by
  parse/run before compute counted); separation coordinate
  switched R2 → measured passage distance pre-run (frozen-form
  change documented here, numbers unaffected in kind).
