#!/bin/zsh
set -e
setopt NULL_GLOB

for taxa_list in ../trees/taxa_*.txt; do
    if [ -e "$taxa_list" ]; then
        echo "Found taxa list: $taxa_list"

        lineage=$(basename "$taxa_list")
        lineage="${lineage#taxa_}"
        lineage="${lineage%.*}"
        echo "Processing lineage: $lineage"
        echo "Using taxa list: $taxa_list"

        echo "Looking for tree files in ../trees/original/*.tre"
        tree_files=(../trees/original/*.tre)
        if [ ${#tree_files[@]} -eq 0 ]; then
            echo "No tree files found in ../trees/original/"
            continue
        fi

        echo "Tree files found: ${tree_files[@]}"

        python3 0b_prune-tree.py ../trees/original/*.tre "$taxa_list" "$lineage" || echo "Error running 0b_prune-tree.py"

        mkdir -p "../trees/$lineage"
        mv ../trees/original/*_"$lineage".tre "../trees/$lineage/" || echo "No files to move for lineage $lineage"
        echo "Moved pruned trees to ../trees/$lineage/"
    else
        echo "No taxa list files found matching pattern: $taxa_list"
    fi
done