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
tree_files = glob.glob("trees/original/*.tre")

# Detect all taxa list files in `trees/taxa_*.txt`
taxa_lists = glob.glob("lists/taxa/taxa_*.txt")
# Extract the lineage name by removing 'taxa_' and file extension
lineages = [Path(t).stem.replace("taxa_", "") for t in taxa_lists]

rule prune_all:
    """
    Aggregation rule that ensures all lineages have been pruned.
    We define the output as the set of all lineage directories,
    indicating that all pruned trees are in place.
    """
    input:
        expand("trees/{lineage}", lineage=lineages)

rule prune_trees:
    """
    For each lineage file, prune every .tre in trees/original/ 
    and move them into trees/{lineage}/
    """
    input:
        taxa_list = "lists/taxa/taxa_{lineage}.txt",
        # Provide the entire list of .tre files as input so the script can prune each
        tree_files = tree_files
    output:
        # We treat the entire lineage directory as the 'output'
        # so Snakemake knows something was produced there.
        directory("trees/{lineage}")
    shell:
        """
        # 1. Run the prune script on all original trees for this lineage.
        python scripts/0b_prune-tree.py {input.tree_files} {input.taxa_list} {wildcards.lineage}

        # 2. Create a directory for pruned trees, if it doesn't exist.
        mkdir -p trees/{wildcards.lineage}

        # 3. Move any newly created pruned *.tre files for this lineage
        #    into trees/{wildcards.lineage}/
        mv trees/original/*_{wildcards.lineage}.tre trees/{wildcards.lineage}/ 2>/dev/null || true
        """