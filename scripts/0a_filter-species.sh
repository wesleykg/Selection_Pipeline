#!/bin/zsh

parallel 'cd {//} && python3 ../../scripts/0a_filter-species.py {/} taxa_*.txt' ::: ../alignments/*/*.fasta
