#!/bin/sh
set -eu
study=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
compiler="$study/build/compiler/root/usr/bin/x86_64-linux-gnu-gfortran-13"
export LD_LIBRARY_PATH="$study/build/compiler/root/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
mkdir -p "$study/build"
cd "$study/build"
"$compiler" --version > compiler_version.txt
"$compiler" -B/usr/lib/gcc/x86_64-linux-gnu/13/ -O1 -std=legacy -ffixed-line-length-132 \
  -fno-automatic -fdefault-real-8 -fdefault-double-8 -fcheck=bounds \
  -ffpe-trap=invalid,zero,overflow -ffunction-sections -fdata-sections \
  "$study/src/sf_probe.f90" "$study/vendor/F1F2IN21_v1.0.f" \
  -Wl,--gc-sections -L"$study/build/compiler/root/usr/lib/gcc/x86_64-linux-gnu/13" \
  -o "$study/build/sf_probe" > build.log 2>&1
