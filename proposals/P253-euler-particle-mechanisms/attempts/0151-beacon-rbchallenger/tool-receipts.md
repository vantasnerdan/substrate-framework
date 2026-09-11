# Tool receipts - 0151 (beacon; moved from 0147 by collision-2 ruling, content unchanged)

- `rb_challenger.py --ext idw`: exit 0. Reproduces 0129
  digit-exact (3.1292/1.6974/1.0373).
- `rb_challenger.py --ext harm`: exit 0. 11.1288/10.8913/
  10.4023 MISS (near-flat in scale).
- `rb_challenger.py --ext bump`: exit 0. 0.1763/0.1339/0.1128
  HOLD-PT (band-localized ~no response).
- `eulerian_front.py`: exit 0 but WITHDRAWN (floor 3.14 at
  sc = 0; 81 sign-flips found+fixed, still dirty). Numbers
  16/15/12 not cited anywhere as results.
- `front_decomp.py`: exit 0. Width sweep δ = 0.1/0.3/1.0 at
  sc = 1.0/0.5/0.25 (dz-share 0.00/0.00/1.00 all scales;
  front|dQ| 0.00/0.02/=exact).
