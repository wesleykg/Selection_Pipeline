#!/bin/zsh

# Enable zsh options for compatibility
setopt extended_glob nullglob

# Define source and target base directories
source_base="../alignments"
target_base="../paml"
test_taxa_dir="../lists/test_taxa"
tree_base="../trees"

# Ensure the target base directory exists
mkdir -p "$target_base"

# Loop through all subdirectories in the source base
for source_dir in "$source_base"/*; do
    if [[ -d "$source_dir" ]]; then
        # Extract the subdirectory name
        subdir_name=$(basename "$source_dir")
        target_dir="$target_base/$subdir_name"
        
        # Ensure the target directory exists
        mkdir -p "$target_dir"
        
        # Loop through non-empty inner directories
        for inner_dir in "$source_dir"/*(/); do
            if [[ -n "$(ls -A "$inner_dir" 2>/dev/null)" ]]; then
                inner_dir_name=$(basename "$inner_dir")
                mkdir -p "$target_dir/$inner_dir_name"
                cp -r "$inner_dir/"* "$target_dir/$inner_dir_name/"
                echo "Copied $inner_dir to $target_dir/$inner_dir_name"
            fi
        done

        # Copy matching test_taxa files
        test_taxa_file="$test_taxa_dir/test_taxa_${subdir_name}.txt"
        if [[ -f "$test_taxa_file" ]]; then
            cp "$test_taxa_file" "$target_dir/"
            echo "Copied $test_taxa_file to $target_dir/"
        else
            echo "No matching test_taxa file for $subdir_name in $test_taxa_dir"
        fi

        # Copy matching tree files
        tree_dir="$tree_base/$subdir_name"
        if [[ -d "$tree_dir" ]]; then
            for tree_file in "$tree_dir"/*_"$subdir_name".tre; do
                if [[ -f "$tree_file" ]]; then
                    cp "$tree_file" "$target_dir/"
                    echo "Copied $tree_file to $target_dir/"
                fi
            done
        else
            echo "No matching tree directory for $subdir_name in $tree_base"
        fi
    fi
done