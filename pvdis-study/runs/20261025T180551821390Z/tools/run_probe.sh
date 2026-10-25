#!/bin/sh
set -eu
study=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
export LD_LIBRARY_PATH="$study/build/compiler/root/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd "$study/build"
./sf_probe < points.dat > sf_values.dat
./sf_probe < points.dat > sf_values_repeat.dat
cmp sf_values.dat sf_values_repeat.dat
printf '%s\n' 'PASS: fresh-process repeat is byte-identical'
if printf '1 4.867 -1 1.263\n' | ./sf_probe > invalid.log 2>&1; then
  echo 'FAIL: invalid Q2 accepted'; exit 1
fi
grep 'outside benchmark domain' invalid.log
printf '%s\n' 'PASS: invalid kinematics rejected'
mkdir -p missing_tables
if (cd missing_tables; printf '1 4.867 .95 1.263\n' | ../sf_probe) > missing_tables.log 2>&1; then
  echo 'FAIL: missing smearing tables accepted'; exit 1
fi
grep 'required smearing table missing' missing_tables.log
printf '%s\n' 'PASS: missing tables rejected before model evaluation'
