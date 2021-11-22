#!/usr/bin/env bash

set -euo pipefail

echo "==> Fetching latest models..."
updated_modules=""
update_details=""
nl="
"

echo " -> Sync (Submodule update)"
git submodule update

for d in */; do
    echo " -> Handling repo $d"
    pushd $d > /dev/null
    fetch_log=$(git fetch origin main 2>&1)
    if [ $? -ne 0 ]; then
        echo "Error fetching. Log: $nl $fetch_log"
        exit 1
    fi
    if [ $(git rev-list HEAD..FETCH_HEAD --count) -gt 0 ]; then
        commit_info="$d: $(git log -1 FETCH_HEAD --oneline --no-decorate)"
        updated_modules="$updated_modules $d"
        update_details="$update_details $nl $commit_info"
        git reset --hard FETCH_HEAD
    fi
    rm -rf common
    ln -s ../common ./
    popd > /dev/null
done

echo "==> Done. Cheching if there were updates..."

if [ "$updated_modules" ]; then
    echo " -> Modules updated. Checking commits are new...!"
    git add $updated_modules
    outp=$(git diff --cached --name-only)
    if [ "$outp" ]; then
        echo " -> New commits found. Creating new bump commit"
        git commit -m "Bump models: $updated_modules $nl $update_details"
        echo " -> Modules were updated. Commit created! Done."
    else
        echo " -> Already up to date. Done"
    fi
else
    echo " -> Modules are already up to date."
fi
