#!/usr/bin/env bash
# This script fetches common mods for local development as a submodule
# It always fetches the latest version in the 'main' branch
# 
# BBP builds always take common into account so dont build on a specific commit
# as this is ignored by BBP deployment
set -euo pipefail

# If we can access common/mod (whatever the nature) we are done]

if [ -d common/mod ]; then
    echo "Common mods found. Checking for update..."
    cd common
    if ! git pull 2> /dev/null; then
        echo "Warning: git pull failed. To update ensure common is in a branch"
    fi
    exit 0
fi

# Otherwise get a fresh clone
rm  -rf common
echo "Common mods not found. Getting a new fresh repo clone"
git clone --depth=2 git@bbpgitlab.epfl.ch:hpc/sim/models/common.git

