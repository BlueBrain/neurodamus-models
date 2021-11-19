#!/usr/bin/env sh

set -euo pipefail

for d in */; do
    echo "==> Handling repo $d"
    pushd $d
    git fetch origin main
    git reset --hard FETCH_HEAD
    rm -rf common
    ln -s ../common ./
    popd
done

echo "==> Done"
