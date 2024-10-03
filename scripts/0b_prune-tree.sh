#!/bin/zsh

for taxa_list in ../trees/*/taxa_*.txt; do
	lineage=$(basename "$taxa_list")
	lineage=${lineage:5:r}
	python3 0b_prune-tree.py ../trees/original/*.tre $taxa_list $lineage
	mv ../trees/original/*_$lineage.tre ../trees/$lineage/
done
