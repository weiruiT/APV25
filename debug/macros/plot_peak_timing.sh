#!/usr/bin/env bash
set -eo pipefail
if [ "$#" -ne 1 ]; then
  echo "Usage: plot_peak_timing.sh /absolute/analysis-output-directory"
  exit 2
fi
macro_dir="$(cd "$(dirname "$0")" && pwd)"
root -l -b -q "$macro_dir/plot_peak_timing.C(\"$1\")"
