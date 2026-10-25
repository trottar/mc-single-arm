#!/bin/sh
set -eu
study=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mkdir -p "$study/build/compiler/packages" "$study/build/compiler/root"
cd "$study/build/compiler/packages"
apt-get download gfortran-13-x86-64-linux-gnu libgfortran-13-dev libgfortran5
for package in ./*.deb; do dpkg-deb -x "$package" ../root; done
../root/usr/bin/x86_64-linux-gnu-gfortran-13 --version
