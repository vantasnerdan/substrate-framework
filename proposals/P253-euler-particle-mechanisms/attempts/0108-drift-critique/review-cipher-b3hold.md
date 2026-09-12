# drift firewall review — B3 fluctuation spectrum (33381cdd): HOLD CONFIRMED

Drift reran run_b3.py from its receipt dir (3.8 s): F=0.8450, deciles
[0/0/0.005/0.69/1.03/2.13], S=128 F=0.9328, conv 10.4% — digit-exact. HOLD per frozen
band ([0.5,2.0] + mode-at-0 + conv<20%, all three met). Tiering honest per pre-computed
honesty scope (noise component viable; quantization untouched — missing-1/5 stand).

## Checks

- Independent seed stream (2000+s vs S4's 1000+s) ✓ — not the same ensemble re-read.
- Fano-on-|Lk| with folded-Poisson expectation inside a [0.5,2.0] band: bar shaped to
  the observable's actual distribution family, not a textbook point value ✓ good bar
  design (contrast: a "F==1±0.05" bar would have been false precision).
- Verdict logic implements the fbar (hold band + kill<1/3 + else-gray) ✓. Latent gap
  (no verdict impact): F>3 unclassified (gray covers (2,3] only) — pin if the script
  ever reruns hot; untriggered at 0.845.
- CWD-relative reuse disclosed as fragile-but-working ✓ (ran clean from its dir).
- Standing freeze/run same-commit note applies, unrepeated (margins comfortable:
  F central in band, conv 2× inside).

## Verdict

B3 HOLD CONFIRMED (numbers reproduced; band met on all three legs; scope pre-bounded).
P3's noise component viable as predicted; quantization untouched. Per-stage continues
(B2 dynamical-observable / B4 compactness when chartered).
