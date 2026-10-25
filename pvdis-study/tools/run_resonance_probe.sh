#!/bin/sh
set -eu
study=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
case_dir=$1
export LD_LIBRARY_PATH="${PVDIS_RUNTIME_LIB:-$study/build/compiler/root/usr/lib/x86_64-linux-gnu}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd "$case_dir"
"$study/build/sf_probe" < points.dat > raw_sf.dat
"$study/build/sf_probe" < points.dat > raw_sf_repeat.dat
cmp raw_sf.dat raw_sf_repeat.dat
printf '%s\n' 'PASS: fresh resonance-study model executions are byte-identical'
