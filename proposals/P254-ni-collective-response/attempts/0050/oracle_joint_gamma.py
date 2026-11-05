"""Only A50 new-boundary finite oracles; no old proof suite replay."""
import json
import sys
from pathlib import Path
from joint_gamma_adapter import HERE
from oracle_fixture import fixture
from oracle_joint_normalization import check as normalization
from oracle_represented_channel import check as represented
from moving_reference_inputs import atomic_json


def main(path):
    c,baseline=fixture()
    result=dict(schema='p254-a50-new-boundary-oracle-v1',normalization=normalization(c),
                represented=represented(c,baseline),
                scope='finite differentiated scalar + normalized functional + original coupled tangent Schur realization + explicit CAR Fock initial ordered noise',
                limits='finite numerical/SymPy evidence, not material calibration or rigorous enclosure')
    atomic_json(Path(path),result)
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__=='__main__':main(sys.argv[1])
