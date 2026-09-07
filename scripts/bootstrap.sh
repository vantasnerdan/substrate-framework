#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if ! command -v pipx >/dev/null 2>&1; then
  echo "ERROR: pipx is required so the memory CLI is available outside a virtual environment." >&2
  echo "Install pipx with your operating-system package manager, then rerun this script." >&2
  exit 1
fi

python3 -m venv "$repo_root/.venv"
"$repo_root/.venv/bin/python" -m pip install -e "$repo_root[dev]"

# The memory command deliberately lives in its own pipx-managed environment.
# It must remain callable without activating the framework development venv.
pipx install --force "$repo_root/tools/agent-memory"

if ! command -v memory >/dev/null 2>&1; then
  echo "ERROR: pipx installed agent-memory, but its binary directory is not on PATH." >&2
  echo "Run 'pipx ensurepath', start a fresh shell, and rerun this script." >&2
  exit 1
fi

memory --version

# ripwire is the structural view named in AGENTS_START_HERE.md: a single
# deterministic binary, pinned like the Lean toolchain so every machine maps the
# tree identically. The upstream installer also activates its agent skills into
# ~/.claude/skills unconditionally, so this fetches the checksummed release
# asset directly and installs only the binary.
ripwire_version="0.4.0"
ripwire_bin="${HOME}/.local/bin/ripwire"
if [ "$("$ripwire_bin" --version 2>/dev/null | awk '{print $2}')" != "$ripwire_version" ]; then
  case "$(uname -s)/$(uname -m)" in
    Linux/x86_64) ripwire_asset="ripwire-${ripwire_version}-linux-x64.tar.gz" ;;
    Linux/aarch64|Linux/arm64) ripwire_asset="ripwire-${ripwire_version}-linux-arm64.tar.gz" ;;
    Darwin/x86_64) ripwire_asset="ripwire-${ripwire_version}-macos-x64.tar.gz" ;;
    Darwin/arm64) ripwire_asset="ripwire-${ripwire_version}-macos-arm64.tar.gz" ;;
    *) echo "ERROR: no prebuilt ripwire for $(uname -s)/$(uname -m); build it from https://github.com/redhat-et/ripwire" >&2; exit 1 ;;
  esac
  ripwire_url="https://github.com/redhat-et/ripwire/releases/download/v${ripwire_version}/${ripwire_asset}"
  ripwire_tmp="$(mktemp -d)"
  trap 'rm -rf "$ripwire_tmp"' EXIT
  curl -fsSL -o "$ripwire_tmp/$ripwire_asset" "$ripwire_url"
  curl -fsSL -o "$ripwire_tmp/$ripwire_asset.sha256" "$ripwire_url.sha256"
  (cd "$ripwire_tmp" && sha256sum -c "$ripwire_asset.sha256")
  tar -xzf "$ripwire_tmp/$ripwire_asset" -C "$ripwire_tmp"
  mkdir -p "$(dirname "$ripwire_bin")"
  install -m 755 "$ripwire_tmp/${ripwire_asset%.tar.gz}/ripwire" "$ripwire_bin"
fi
if ! command -v ripwire >/dev/null 2>&1; then
  echo "ERROR: ripwire is installed at $ripwire_bin but that directory is not on PATH." >&2
  exit 1
fi
ripwire --version

"$repo_root/scripts/setup_lean.sh"

echo "Bootstrap complete. No virtual-environment activation is needed for memory or repository Lean commands."
echo "Run: $repo_root/scripts/validate.sh"
echo "Formal developments: $repo_root/scripts/check_lean.sh"
