#!/bin/sh

# Use Mathematica's kernel directly.  The standalone macOS WolframScript
# launcher can spin before starting a kernel while opening its SharedMemory
# WSTP transport; in that state BuildRelease.wls has not begun executing.

set -eu

script_directory=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "$script_directory/../.." && pwd)
kernel=${WOLFRAM_KERNEL:-/Applications/Mathematica.app/Contents/MacOS/WolframKernel}

if [ ! -x "$kernel" ]
then
  printf '%s\n' "WolframKernel not found or not executable: $kernel" >&2
  exit 2
fi

cd "$project_root"

exec env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$kernel" -noinit -script \
  "$script_directory/BuildRelease.wls" "$project_root"
