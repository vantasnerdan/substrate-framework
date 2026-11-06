"""Read complete saved joint Gamma numeric kernels; no source replay."""
import json
from pathlib import Path
import numpy as np
from joint_material_inputs import file_binding,authenticate
from represented_electron_channel import polarization_difference


def consume(record_path):
    record=json.loads(Path(record_path).read_bytes())
    outputs={}
    with np.load(authenticate(record['common_archive']),allow_pickle=False) as saved:
        for component in ('greater','lesser','retarded'):
            pi=saved['joint_pi_'+component]
            ionic=saved['ionic_pi_'+component]
            sigma=saved['joint_sigma_'+component]
            sigma0=saved['ionic_sigma_'+component]
            outputs[component]=dict(joint_pi_norm=float(np.linalg.norm(pi)),
                ionic_pi_norm=float(np.linalg.norm(ionic)),
                mixed_ion_field_norm=float(np.linalg.norm(pi[:,:,:111,111:])),
                mixed_field_ion_norm=float(np.linalg.norm(pi[:,:,111:,:111])),
                field_field_norm=float(np.linalg.norm(pi[:,:,111:,111:])),
                joint_minus_ionic_sigma_norm=float(np.linalg.norm(sigma-sigma0)),
                ionic_pi_retention_error=float(np.linalg.norm(pi[:,:,:111,:111]-ionic)))
        new={name:saved['joint_'+name].copy() for name in ('pi_retarded','pi_greater','pi_lesser')}
    with np.load(authenticate(record['represented_archive']),allow_pickle=False) as saved:
        represented={name:saved[name] for name in new}
        difference=polarization_difference(new,represented)
        delta={name:float(np.linalg.norm(value)) for name,value in difference.items()}
        sensitivity={name:float(np.linalg.norm(new[name]+represented[name]-difference[name])) for name in new}
        noise_norm=float(np.linalg.norm(saved['ordered_ab']))
        response_norm=float(np.linalg.norm(saved['pi_retarded']))
    return dict(schema='p254-a50-saved-numeric-consumer-v1',source=file_binding(record_path),
        times=record['times'],scope='all114 Gamma joint component and same128-active-sector represented channel; no closure',
        common=outputs,new_minus_represented_pi_norms=delta,
        stacked_instead_of_subtracted_differences=sensitivity,
        represented_response_norm=response_norm,represented_ordered_noise_norm=noise_norm)


if __name__=='__main__':
    import sys
    from moving_reference_inputs import atomic_json
    result=consume(sys.argv[1]);atomic_json(sys.argv[2],result)
    print(json.dumps(result,allow_nan=False),flush=True)
