#!/bin/zsh

setopt +o nomatch  # Adjust globbing behavior

# Remove old issues.txt files, and doesn't throw an error if it doesn't exist yet
rm -f ../alignments/*/issues.txt(N)

# Run the check_frame.py script in Parallel
parallel 'cd {//} && python3 ../../scripts/1a_check-frame.py {/} >> issues.txt' ::: ../alignments/*/*.fasta

# Sort and remove duplicates from the issues.txt files
parallel 'cd {} && if [ -f issues.txt ]; then sort -u issues.txt -o issues.txt; fi' ::: ../ ../alignments/*/