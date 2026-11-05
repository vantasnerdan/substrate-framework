"""Reuse original Gaussian Cartesian propagation without any electron/H replay."""
import numpy as np
from material_electron_ion_inputs import FrozenReference


class CartesianReference:
    """Original raw quadratic S, original W0, all coordinates including Gamma field.

    Only ionic readout is taken after FULL phase-space propagation. No new
    preparation, stable-mode clipping, stationarity or FDT is introduced.
    This fixed-reference evaluation is explicitly not a selfconsistent history.
    """
    label = 'original-Gaussian-raw-quadratic-Cartesian-fixed-reference'
    transitions = FrozenReference.transitions
    boson = FrozenReference.boson

    def __init__(self,coefficients,times=(0.,1.,10.,30.)):
        self.c = coefficients; self.times = np.asarray(times,float)
        self._q = None; self._transitions = None
        self.bindings = []

    def ordered_grid(self,q):
        greater = []; lesser = []
        for a in range(len(self.times)):
            grow = []; lrow = []
            for b in range(len(self.times)):
                g,l = self.boson(q,a,b)
                grow.append(1j*g); lrow.append(1j*l)
            greater.append(grow); lesser.append(lrow)
        return np.asarray(greater),np.asarray(lesser)
