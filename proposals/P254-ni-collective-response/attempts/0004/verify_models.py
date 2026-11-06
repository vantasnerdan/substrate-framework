"""Check actual pseudovalence spaces against the completed native calculation."""

from collections import Counter
import gzip
import json
from pathlib import Path
import re
import sys
import traceback
import xml.etree.ElementTree as ET

import sympy as sp

from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
PROPOSAL = HERE.parents[1]
# Reuse the already executed evidence writer, rather than another receipt convention.
sys.path.insert(0,str(HERE.parent / "0003"))
from oracle_artifacts import sha256, save_receipt


def upf_valence(path):
    header = re.search(r"<PP_HEADER>(.*?)</PP_HEADER>",path.read_text(),re.S).group(1)
    value = re.search(r"(?m)^\s*(\S+)\s+Z valence\s*$",header).group(1)
    return sp.Rational(value)


def source_models(ledger, metadata, receipt):
    qe = {}
    paw = {}
    for element,record in metadata["qe_gbrv"].items():
        path = Path(record["path"])
        ledger.check(f"actual {element} UPF hash",sha256(path) == record["sha256"])
        qe[element] = upf_valence(path)
        ledger.check(f"actual {element} UPF metadata",qe[element] == record["valence"])
        receipt["source_sha256"][str(path)] = sha256(path)
    for name,record in metadata["gpaw"].items():
        path = Path(record["path"])
        ledger.check(f"actual {name} PAW hash",sha256(path) == record["sha256"])
        with gzip.open(path,"rt") as handle:
            atom = ET.fromstring(handle.read()).find("atom")
        paw[name] = sp.Rational(atom.attrib["valence"])
        ledger.check(f"actual {name} PAW metadata",paw[name] == record["valence"])
        receipt["source_sha256"][str(path)] = sha256(path)
    return qe,paw


def native_count_block(ledger, record, qe, paw, receipt):
    result_path = Path(record["inputs"]["completed_scf"])
    xml_path = Path(record["inputs"]["native_xml"])
    result = json.loads(result_path.read_text())
    ledger.check("native source SCF is complete",result["status"] == "complete")
    stdout_record = result["data"]["artifacts"]["pw.stdout"]
    stdout_path = Path(stdout_record["path"])
    ledger.check("published native stdout hash",sha256(stdout_path) == stdout_record["sha256"])
    counts = re.findall(r"(?m)^\s*number of electrons\s*=\s*(\S+)\s*$",stdout_path.read_text())
    ledger.check("one native electron-count field",len(counts) == 1)
    published_count = sp.Rational(counts[0])
    ledger.check("supplemental XML pinned bytes",sha256(xml_path) == record["inputs"]["native_xml_sha256"])
    ledger.check("supplemental XML remains in completed attempt",
                 xml_path.resolve().is_relative_to(Path(result["data"]["attempt_directory"]).resolve()))
    native = ET.parse(xml_path).getroot()
    composition = Counter(atom.attrib["name"] for atom in native.findall("input/atomic_structure/atomic_positions/atom"))
    ledger.check("XML composition matches published SCF symbols",composition == Counter(result["data"]["symbols"]))
    charge = sp.Rational(native.findtext("input/bands/tot_charge"))
    expected = sum(count*qe[element] for element,count in composition.items())-charge
    reported = published_count
    ledger.check("XML electron count matches published native stdout",
                 sp.Rational(native.findtext("output/band_structure/nelec")) == reported)
    ledger.check("exact UPF composition matches native electron count",expected == reported)
    ledger.check("spin2 fixture is the actual collinear LSDA source",native.findtext("input/spin/lsda") == "true")
    wrong_spin = 2*(expected+charge)-charge
    ledger.check("extra spin degeneracy mutation exposed",wrong_spin != reported)
    mutants = {"extra_spin_degeneracy":str(wrong_spin)}
    for name,valence in paw.items():
        wrong = sum(count*(valence if element == "Ni" else qe[element])
                    for element,count in composition.items())-charge
        ledger.check(f"{name} substitution into QE space exposed",wrong != reported)
        mutants[name] = str(wrong)
    receipt["source_sha256"].update({str(result_path):sha256(result_path),str(xml_path):sha256(xml_path),
                                    str(stdout_path):sha256(stdout_path)})
    return {"composition":dict(composition),"net_charge":str(charge),"qe_valence_by_element":{k:str(v) for k,v in qe.items()},
            "paw_Ni_valence_by_dataset":{k:str(v) for k,v in paw.items()},
            "native_reported_electrons":str(reported),"independent_UPF_composition_electrons":str(expected),
            "exposed_wrong_model_counts":mutants}


def main():
    record_path = HERE / "record.json"
    record = json.loads(record_path.read_text())
    metadata_path = PROPOSAL / record["inputs"]["model_metadata"]
    metadata = json.loads(metadata_path.read_text())
    receipt = {"status":"running","source_sha256":{str(Path(__file__)):sha256(__file__),
                str(record_path):sha256(record_path),str(metadata_path):sha256(metadata_path),
                str(HERE.parent/"0003/oracle_artifacts.py"):sha256(HERE.parent/"0003/oracle_artifacts.py")},
               "scope":"Source-bound pseudovalence identity only, not a material response correspondence or a stationary structure."}
    ledger = CheckLedger("P254-electronic-model-identity")
    try:
        qe,paw = source_models(ledger,metadata,receipt)
        receipt["model_counts"] = native_count_block(ledger,record,qe,paw,receipt)
        ledger.finish()
        receipt["status"] = "complete"
    except Exception:
        receipt["status"] = "failed"
        receipt["raw_exception"] = traceback.format_exc()
        raise
    finally:
        receipt["passed_checks"] = ledger.passed
        save_receipt(HERE/"receipt.json",receipt)
    print("P254_ELECTRONIC_MODEL_IDENTITY_COMPLETE",flush=True)


if __name__ == "__main__":
    main()
