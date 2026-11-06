"""Own one bounded native construction and publish all-exit receipts/notification."""
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
WORK = Path('/home/administrator/projects/catalyst-nc-operator/runs/2026-10-19_collective-response/coherence-action-0014')


def exit_receipt(unit):
    WORK.mkdir(parents=True, exist_ok=True)
    outcome = {'unit': unit, 'time': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'invocation_id': os.environ.get('INVOCATION_ID'),
        'service_result': os.environ.get('SERVICE_RESULT'),
        'exit_code': os.environ.get('EXIT_CODE'), 'exit_status': os.environ.get('EXIT_STATUS'),
        'meaning': 'Bounded NiH action construction exit; not coherence-domain or campaign completion.'}
    path = WORK / ('outcome-' + str(outcome['invocation_id']) + '.json')
    path.write_text(json.dumps(outcome, indent=2) + '\n')
    tag = 'dft-done' if outcome['service_result'] == 'success' else 'dft-fail'
    notify = '/home/administrator/projects/catalyst-dft/runs/campaign3/notify.py'
    result = subprocess.run(['/usr/bin/python3', notify, tag, unit, str(path)],
                            capture_output=True, timeout=30)
    (WORK / 'notification.stdout').write_bytes(result.stdout)
    (WORK / 'notification.stderr').write_bytes(result.stderr)
    (WORK / 'notification-exit.json').write_text(json.dumps({'exit': result.returncode})+'\n')


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--notify-exit':
        exit_receipt(sys.argv[2])
    else:
        from run_native import main
        main()
