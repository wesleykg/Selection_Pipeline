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
        expand("alignments/{subdir}/{alignment}_{subdir}_filtered.fasta", 
               subdir=[s[0] for s in alignments], 
               alignment=[s[1] for s in alignments])

rule filter:
    input:
        "alignments/{subdir}/{alignment}.fasta",
        "alignments/{subdir}/taxa_{subdir}.txt"
    output:
        "alignments/{subdir}/{alignment}_{subdir}_filtered.fasta"
    shell:
        "scripts/0a_filter-species.py {input} {output}"