#!/usr/bin/env bash

set -euo pipefail

declare -A MODEL_BRANCH
MODEL_BRANCH["common"]=${MODELS_COMMON_BRANCH:-main}
MODEL_BRANCH["neocortex"]=${NEOCORTEX_BRANCH:-main}
MODEL_BRANCH["hippocampus"]=${HIPPOCAMPUS_BRANCH:-main}
MODEL_BRANCH["thalamus"]=${THALAMUS_BRANCH:-main}
MODEL_BRANCH["mousify"]=${MOUSIFY_BRANCH:-main}


if [ -n "$(git ls-files -m)" ]; then
    echo "Error: Repo is not clean. Please commit or discard ongoing changes"
    echo " hint: You may need to 'git submodule update' to discard changing submodule commits"
    exit 1
fi

echo "==> Fetching latest models..."
updated_modules=""
update_details=""
nl="
"

for d in */; do
    model_name="${d%/}"
    if [ ${d::1} == "_" ]; then continue
    fi
    echo "Repo: $model_name"
    pushd $d > /dev/null
    printf " -> Fetching tip of branch ${MODEL_BRANCH[$model_name]}..."
    fetch_log=$(git fetch origin ${MODEL_BRANCH[$model_name]} 2>&1) || {
        echo "[Error] Log: $nl$fetch_log"
        exit 1
    }
    if [ $(git rev-list HEAD..FETCH_HEAD --count) -gt 0 ]; then
        commit_info="**$model_name**:
$(git log HEAD..FETCH_HEAD --oneline --no-decorate --ancestry-path | sed 's/^/  - /')"
        updated_modules="$updated_modules $model_name"
        update_details="$update_details$nl$commit_info$nl"
        git reset --hard FETCH_HEAD
    else
        echo " [ok] Already up to date."
    fi
    git checkout -B ${MODEL_BRANCH[$model_name]} >& /dev/null  # Move to a branch of that name
    popd > /dev/null
done

echo "==> Done. Cheching if there were updates..."

if [ "$updated_modules" ]; then
    echo " -> Modules updated. Checking commits are new...!"
    git add $updated_modules
    outp=$(git diff --cached --name-only)
    if [ "$outp" ]; then
        echo " -> New commits found. Creating new bump commit"
        git commit -m "Bumping models: $updated_modules$nl$update_details"
        echo " -> Modules were updated. Commit created! Done."
    else
        echo " -> Already up to date. Done"
    fi
else
    echo " -> All modules are up to date."
fi
