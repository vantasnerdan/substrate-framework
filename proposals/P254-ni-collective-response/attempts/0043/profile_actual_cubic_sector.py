"""Profile exactly one missing actual sector, retaining its normal immutable output."""
import cProfile
from pathlib import Path
import pstats
import sys
import cubic_material_source as source
import run_material_cubic_kernel as cli

original=source._new_sector
calls=0

def measured(*args,**kwargs):
    global calls
    if calls:
        raise SystemExit('Intentional checkpoint after one newly committed actual sector; full source remains partial')
    calls+=1
    profile=cProfile.Profile()
    try:
        return profile.runcall(original,*args,**kwargs)
    finally:
        path=Path(__file__).with_name('actual-sector-cost-001.prof')
        profile.dump_stats(path)
        pstats.Stats(profile,stream=sys.stdout).sort_stats('cumulative').print_stats(35)
        sys.stdout.flush()

source._new_sector=measured
cli.main()
