#!/bin/zsh

parallel "cd {//} && python3 ../../scripts/0a_filter.py {/} taxa_*.txt" ::: alignments/*/*.fasta
