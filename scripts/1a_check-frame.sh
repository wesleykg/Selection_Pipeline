#!/bin/zsh

rm alignments/*/issues.txt
parallel "cd {//} && sed 's/-/N/g' {/} > temp_{/} && python3 ../../scripts/1a_check-frame.py temp_{/} >> issues.txt && rm temp_{/}" ::: alignments/*/*.fasta