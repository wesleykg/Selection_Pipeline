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

        echo "Working on $alignment"

        # Run trimal with error handling
        if ! trimal -in "$alignment" -out "$alignment.noallgaps" -noallgaps -keepseqs; then
            echo "Error: trimal failed on $alignment" >&2
            continue
        fi

        # Process the trimal output file with error handling
        if [ -e "$alignment.noallgaps" ]; then
            if ! sed -E 's/ [0-9]+ bp//g' "$alignment.noallgaps" > "$alignment"; then
                echo "Error: sed failed on $alignment.noallgaps" >&2
                continue
            fi
            rm "$alignment.noallgaps"
        else
            echo "Warning: trimal couldn't find $alignment.noallgaps" >&2
        fi
    done
    popd > /dev/null || exit
done