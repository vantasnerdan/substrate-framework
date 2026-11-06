"""Original Coulomb inputs and public A63 saved electronic legs only.

No FrozenReference, DressedHistory or mechanical/Gaussian bath is constructed.
A63 SavedHistory still authenticates/maps its complete common-history provider;
only its electron retarded/lesser banks enter the unchanged EE contraction.
"""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / 'src')] + [str(HERE.parent / n) for n in (
    '0063', '0062', '0061', '0060', '0056', '0055', '0054', '0053',
    '0052', '0051', '0050', '0048', '0046', '0045', '0044', '0043',
    '0042', '0041', '0030', '0029', '0028', '0027', '0025', '0024', '0064')]
import numpy as np
from history_provider import SavedHistory
from dressed_reference import greater
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors


class ElectronHistory:
    """Own A63 mappings; expose A56 greater/lesser over time-major views."""
    def __init__(self, provider):
        self.saved = SavedHistory(provider)
        self.ret = self.less = None
        try:
            self.times = self.saved.times
            if not np.array_equal(self.times, np.arange(31, dtype=float)):
                raise ValueError('exact complete saved31 time grid required')
            bank = self.saved.banks['electron']
            # Transpose returns views of the read-only slot-first NPY mappings.
            # A56 owns the spectral algebra, including the raw theta(0)=1 jump.
            self.ret = bank['ret'].transpose(1, 2, 0, 3, 4)
            self.less = bank['less'].transpose(1, 2, 0, 3, 4)
            if (self.ret.dtype != np.complex128 or self.less.dtype != np.complex128
                    or self.ret.flags.writeable or self.less.flags.writeable):
                raise ValueError('read-only complex128 saved electron banks required')
            self.bindings = self.saved.bindings
        except BaseException:
            self.close()
            raise

    def electron(self, a, b):
        return greater(self.ret, self.less, a, b), self.less[a, b]

    def close(self):
        self.ret = self.less = None
        self.saved.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def load_inputs(provider):
    """Return (BareCoulombFactors, ElectronHistory); caller owns history.close()."""
    factors = BareCoulombFactors(load_coulomb_inputs())
    return factors, ElectronHistory(provider)
