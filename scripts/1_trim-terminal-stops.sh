#!/bin/zsh

parallel 'cd {//} && python3 ..//../scripts/1_trim-terminal-stops.py {/}' ::: ../alignments/*/*.fasta