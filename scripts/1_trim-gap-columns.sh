#!/bin/zsh

set -euo pipefail

ALIGNMENTS_DIR="../alignments"

# Iterate over each subdirectory in the alignments directory
for dir in "$ALIGNMENTS_DIR"/*/; do
    if [ ! -d "$dir" ]; then
        echo "Skipping $dir: Not a directory" >&2
        continue
    fi
    pushd "$dir" > /dev/null || continue

    # Process each .fasta file
    for alignment in *.fasta; do
        if [ ! -f "$alignment" ]; then
            echo "No .fasta files found in $dir" >&2
            continue
        fi

        echo "Working on $alignment in $dir"

        # Run trimal
        trimal -in "$alignment" -out "$alignment.noallgaps" -noallgaps -keepseqs
        if [ $? -ne 0 ]; then
            echo "Error: trimal failed on $alignment" >&2
            continue
        fi

        # Process the output file
        if [ -e "$alignment.noallgaps" ]; then
            sed -E 's/ [0-9]+ bp//g' "$alignment.noallgaps" > "$alignment"
            rm "$alignment.noallgaps"
        else
            echo "Warning: $alignment.noallgaps not found after trimal" >&2
        fi
    done
    popd > /dev/null || exit
done