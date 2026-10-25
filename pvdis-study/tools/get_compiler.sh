#!/bin/sh
set -eu
study=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mkdir -p "$study/build/compiler/packages" "$study/build/compiler/root"
cd "$study/build/compiler/packages"
apt-get download \
  gfortran-13-x86-64-linux-gnu=13.3.0-6ubuntu2~24.04.1 \
  libgfortran-13-dev=13.3.0-6ubuntu2~24.04.1 \
  libgfortran5=14.2.0-4ubuntu2~24.04.1
sha256sum -c "$study/data/compiler_packages.sha256"
while read -r digest package; do
  dpkg-deb -x "$package" ../root
done < "$study/data/compiler_packages.sha256"
../root/usr/bin/x86_64-linux-gnu-gfortran-13 --version
