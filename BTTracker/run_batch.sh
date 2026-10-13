#!/usr/bin/env bash
set -eo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Discovery and --help also work on machines without CERN ROOT.
needs_root=true
for arg in "$@"; do
    case "$arg" in --dry-run|--help|-h) needs_root=false;; esac
done
analysis_python="${BTTRACKER_PYTHON:-python3}"
if "$needs_root"; then
    if ! command -v root-config >/dev/null 2>&1 || \
       ! "$analysis_python" -c 'import ROOT' >/dev/null 2>&1; then
        source "${BTTRACKER_ROOT_SETUP:-/home/thoriatong/root/bin/thisroot.sh}"
    fi
    if ! "$analysis_python" -c 'import ROOT' >/dev/null 2>&1; then
        if [[ -x /usr/bin/python3 ]] && \
           /usr/bin/python3 -c 'import ROOT' >/dev/null 2>&1; then
            analysis_python=/usr/bin/python3
        else
            echo "ERROR: no Python interpreter with PyROOT is available." >&2
            exit 1
        fi
    fi
    export BTTRACKER_PYTHON="$analysis_python"
fi
exec "$analysis_python" "$script_dir/scripts/run_manager.py" batch "$@"
