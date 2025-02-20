from pathlib import Path
import glob

# Get all subdirectories inside alignments/
subdirs = [Path(d) for d in glob.glob("alignments/*") if Path(d).is_dir()]

# Find all `.fasta` files in those subdirectories and extract alignment names
alignment_files = [f for subdir in subdirs for f in subdir.glob("*.fasta")]
# Extract directory names and alignment names
alignments = [(f.parent.name, f.stem) for f in alignment_files]

rule filter_all:
    input:
        expand(
            "alignments/{subdir}/{alignment}_{subdir}_filtered.fasta", 
            subdir=[s[0] for s in alignments], 
            alignment=[s[1] for s in alignments]
        )

rule filter:
    input:
        "alignments/{subdir}/{alignment}.fasta",
        "lists/taxa/taxa_{subdir}.txt"
    output:
        "alignments/{subdir}/{alignment}_{subdir}_filtered.fasta"
    shell:
        "scripts/0a_filter-species.py {input} {output}"

# Detect all tree files in `trees/original/`
original_tree_file = glob.glob("trees/original/*.tre")

# Detect all taxa list files in `trees/taxa_*.txt`
taxa_lists = glob.glob("lists/taxa/taxa_*.txt")
# Extract the lineage name by removing 'taxa_' and file extension
lineages = [Path(t).stem.replace("taxa_", "") for t in taxa_lists]

rule prune_all:
    input:
        expand("trees/{lineage}", lineage=lineages)

rule prune_trees:
    input:
        original_tree_file=original_tree_file,
        taxa_list="lists/taxa/taxa_{lineage}.txt"
    output:
        directory("trees/{lineage}")
    shell:
        """
        python scripts/0b_prune-tree.py {input.original_tree_file} {input.taxa_list} {wildcards.lineage}

        mkdir -p trees/{wildcards.lineage}

        mv trees/original/*_{wildcards.lineage}.tre trees/{wildcards.lineage}/ 2>/dev/null || true
        """