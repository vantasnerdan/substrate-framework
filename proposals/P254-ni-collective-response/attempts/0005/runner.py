"""Compose checked finite-matrix blocks and atomically preserve source-bound receipts."""
import argparse
import importlib
import json
import os
from pathlib import Path
import signal
import time
import traceback
from checkpoints import atomic, recover, save
from provenance import manifest, loaded_runtime

BASE = Path(__file__).resolve().parent
BLOCKS = ("native", "transformations", "embedding", "numerical")



def interrupted(signum, frame):
    raise InterruptedError("Explicit service interruption: signal " + str(signum))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--block", choices=BLOCKS)
    parser.add_argument("--resume-from", help="Explicit prior run ID; never automatically retry")
    args = parser.parse_args()
    directory = BASE / "runs" / args.run_id
    directory.mkdir(parents=True, exist_ok=False)
    receipt = {"status": "running", "started_unix": time.time()}
    atomic(directory / "started.json", receipt)
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        receipt["provenance"] = manifest()
        atomic(directory / "started.json", receipt)
        assert receipt["provenance"]["affinity"] == [5], "CPU5 isolation is required"
        receipt["blocks"] = {}
        receipt["recovered"] = {}
        for name in (args.block,) if args.block else BLOCKS:
            restored = recover(BASE / "runs" / args.resume_from, name, receipt["provenance"]) if args.resume_from else None
            if restored:
                receipt["blocks"][name] = restored["result"]
                receipt["recovered"][name] = {k: v for k, v in restored.items() if k != "result"}
                print("RECOVERED", name, restored["source"], flush=True)
            else:
                print("RUNNING", name, flush=True)
                receipt["blocks"][name] = importlib.import_module("blocks." + name).run()
                print("VERIFIED", name, json.dumps(receipt["blocks"][name], default=str), flush=True)
            save(directory, name, receipt["blocks"][name], receipt["provenance"])
        receipt["status"] = "complete"
    except BaseException:
        receipt["status"] = "failed"
        receipt["traceback"] = traceback.format_exc()
        (directory / "traceback.txt").write_text(receipt["traceback"])
        raise
    finally:
        receipt["ended_unix"] = time.time()
        receipt["elapsed_seconds"] = receipt["ended_unix"] - receipt["started_unix"]
        receipt["loaded_runtime_sha256"] = loaded_runtime()
        atomic(directory / "receipt.json", receipt)
        print("RECEIPT", directory / "receipt.json", receipt["status"], flush=True)


if __name__ == "__main__":
    main()
