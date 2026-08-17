#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/validate.sh
  scripts/validate.sh --full
  scripts/validate.sh --pytest-scope SELECTOR [SELECTOR ...]

Every mode runs the repository, generated-state, memory, skill, import, and
compile checks. With no arguments or --full, pytest runs the complete suite.
Arguments after --pytest-scope must be repository test files, directories, or
node IDs. Pytest options are rejected so collection-only or similar flags cannot
be mistaken for executed validation.
EOF
}

pytest_mode="full"
pytest_args=()
case "${1:-}" in
  "")
    ;;
  --full)
    if [ "$#" -ne 1 ]; then
      echo "ERROR: --full does not accept additional arguments" >&2
      usage >&2
      exit 2
    fi
    ;;
  --pytest-scope)
    shift
    if [ "$#" -eq 0 ]; then
      echo "ERROR: --pytest-scope requires at least one pytest selector" >&2
      usage >&2
      exit 2
    fi
    pytest_mode="scoped"
    pytest_args=("$@")
    for selector in "${pytest_args[@]}"; do
      case "$selector" in
        -* )
          echo "ERROR: pytest options are not allowed in --pytest-scope: $selector" >&2
          exit 2
          ;;
        tests|tests/*|tools/agent-memory/tests|tools/agent-memory/tests/*)
          ;;
        *)
          echo "ERROR: pytest scope must select repository tests: $selector" >&2
          exit 2
          ;;
      esac
    done
    ;;
  -h|--help)
    usage
    exit 0
    ;;
  *)
    echo "ERROR: unknown validation option: $1" >&2
    usage >&2
    exit 2
    ;;
esac

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -n "${PYTHON:-}" ]; then
  python_bin="$PYTHON"
elif [ -x "$repo_root/.venv/bin/python" ]; then
  python_bin="$repo_root/.venv/bin/python"
else
  python_bin="python3"
fi

PYTHONPATH="$repo_root/src" "$python_bin" "$repo_root/scripts/validate_repository.py"
PYTHONPATH="$repo_root/src" "$python_bin" "$repo_root/scripts/render_docs.py" --check
PYTHONPATH="$repo_root/src" "$python_bin" "$repo_root/scripts/render_memory.py" --check
"$python_bin" -c 'import numpy, scipy, sympy, yaml'
if ! command -v memory >/dev/null 2>&1; then
  echo "ERROR: memory is not on PATH; run scripts/bootstrap.sh" >&2
  exit 1
fi
bundled_memory_version="$(PYTHONPATH="$repo_root/tools/agent-memory/src" "$python_bin" -c 'from agent_memory import __version__; print(__version__)')"
installed_memory_version="$(memory --version)"
case "$installed_memory_version" in
  *"$bundled_memory_version"*) ;;
  *)
    echo "ERROR: installed memory CLI does not match bundled version $bundled_memory_version: $installed_memory_version" >&2
    exit 1
    ;;
esac
memory --help >/dev/null
memory validate "$repo_root/memory"
"$python_bin" "$repo_root/.agents/skills/physics-erdos-loop/scripts/validate_skill.py" \
  "$repo_root/.agents/skills/physics-erdos-loop"
"$python_bin" -m compileall -q "$repo_root/src" "$repo_root/tools/agent-memory/src"

# Pin the pytest collection count so every release manifest and campaign
# adjudication quotes one source of truth for the test metric. `pytest
# --collect-only -q` ends with "<N> tests collected in <T>s"; extract the
# numeric prefix on the matching line and emit it as "Tests: N" before the
# pytest run so the metric is reproducible regardless of selector.
tests_collected="$(
  PYTHONPATH="$repo_root/src" "$python_bin" -m pytest --collect-only -q 2>/dev/null \
    | awk '/^[0-9]+ tests collected/ {print $1; exit}'
)"
tests_collected="${tests_collected:-0}"
echo "Tests: ${tests_collected}"

if [ "$pytest_mode" = "full" ]; then
  PYTHONPATH="$repo_root/src" "$python_bin" -m pytest -q
  echo "ALL REPOSITORY WORKFLOW CHECKS PASS (full pytest suite, ${tests_collected} tests)"
else
  printf 'Running requested pytest scope:'
  printf ' %q' "${pytest_args[@]}"
  printf '\n'
  PYTHONPATH="$repo_root/src" "$python_bin" -m pytest -q -- "${pytest_args[@]}"
  echo "ALL FIXED REPOSITORY CHECKS AND REQUESTED PYTEST SCOPE PASS (collection: ${tests_collected} tests)"
fi
