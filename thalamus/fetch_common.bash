#!/usr/bin/env bash
# This script fetches common mods for local development as a submodule
# It always fetches the latest version in the 'main' branch
#
# BBP builds always take common into account so dont build on a specific commit
# as this is ignored by BBP deployment
set -euo pipefail


pull_or_clone() {
    local MODEL_NAME=$1
    local MODEL_PATH=$2
    if [ -d $MODEL_PATH/mod ]; then
        echo "$MODEL_NAME mods found. Checking for update..."
        pushd $MODEL_PATH
        if ! git pull 2> /dev/null; then
            echo "Warning: git pull failed. To update ensure $MODEL_NAME is in a branch"
        fi
        popd > /dev/null
        return 0
    fi

    # Otherwise get a fresh clone
    rm  -rf $MODEL_PATH
    echo "$MODEL_NAME mods not found. Getting a new fresh repo clone"
    pushd $(dirname $MODEL_PATH)
    git clone --depth=2 git@bbpgitlab.epfl.ch:hpc/sim/models/$MODEL_NAME
    popd > /dev/null
}

pull_or_clone common common
pull_or_clone neocortex deps/neocortex
