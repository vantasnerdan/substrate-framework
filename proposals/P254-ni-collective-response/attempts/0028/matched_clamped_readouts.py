"""Bounded spatial/field Gaussian readouts, with explicit field multiplicity."""
import numpy as np
from matched_clamped_regions import IndependentRegions


def connected_readouts(model, left, right, anchor, coordinates):
    """Ionic R vs anchor; field R vs anchor-field for independent regions.

    Periodic field covariance has a singleton leading field index; its mixed
    covariance is identical at every cell because only Gamma sees the field.
    Independent regions instead have N local fields; correlations between
    different regions vanish exactly in the tensor-product Gaussian state.
    """
    ix = coordinates
    independent = isinstance(model, IndependentRegions)
    if independent:
        mask = np.all(model.indices == model.cell(anchor),axis=1).astype(float)[:,None,None]
        local = model.local
        equal = mask*local.ordered_pair(left,left,(0,0,0),(0,0,0),ix)
        unequal = mask*local.ordered_pair(left,right,(0,0,0),(0,0,0),ix)
        reverse = mask*local.ordered_pair(right,left,(0,0,0),(0,0,0),ix).T
        fequal = mask*local.field_ordered_pair(left,left)
        funequal = mask*local.field_ordered_pair(left,right)
        mixed = mask*local.ion_field_ordered_pair(left,right,ix)
    else:
        equal = np.asarray([model.ordered_pair(left,left,cell,anchor,ix) for cell in model.indices])
        unequal = np.asarray([model.ordered_pair(left,right,cell,anchor,ix) for cell in model.indices])
        reverse = np.asarray([model.ordered_pair(right,left,anchor,cell,ix).T for cell in model.indices])
        fequal = model.field_ordered_pair(left,left)[None,:,:]
        funequal = model.field_ordered_pair(left,right)[None,:,:]
        mixed = np.broadcast_to(model.ion_field_ordered_pair(left,right,ix),(model.n,len(ix),2*model.f))
    return dict(spatial_equal_time_ordered=equal, spatial_two_time_ordered=unequal,
        spatial_two_time_commutator=unequal-reverse,
        field_equal_time_ordered=fequal, field_two_time_ordered=funequal,
        ion_field_two_time_ordered=mixed)


def time_readout(model, left, right, anchor, coordinate, kicks, initial_energy):
    baseline = model.mean(left)
    independent = isinstance(model,IndependentRegions)
    field = lambda value: value if independent else value[None,:]
    arrays = dict(baseline_ions=baseline['ions'],baseline_field=field(baseline['field']))
    diagnostics = {}
    for name, kick in kicks.items():
        value = baseline if kick is None else model.mean(left,kick)
        if kick is not None:
            arrays[name+'_kick_delta_ions'] = value['ions']-baseline['ions']
            arrays[name+'_kick_delta_field'] = field(value['field']-baseline['field'])
        record = model.diagnostics(left,value['q'])
        record['energy_drift_ry'] = record['energy_ry']-initial_energy[name]
        record['energy_drift_relative'] = abs(record['energy_drift_ry'])/max(1.,abs(initial_energy[name]))
        diagnostics[name] = record
    arrays.update(connected_readouts(model,left,right,anchor,[coordinate,model.d+coordinate]))
    return arrays, diagnostics
