import glob
import os

# Get all subdirectories inside alignments/
subdirs = [d for d in glob.glob("alignments/*") if os.path.isdir(d)]

# Find all `.fasta` files in those subdirectories and extract alignment names
alignment_files = []
for subdir in subdirs:
    alignment_files.extend(glob.glob(f"{subdir}/*.fasta"))

# Extract directory names and alignment names
alignments = [
    (os.path.basename(os.path.dirname(f)), os.path.basename(f).replace(".fasta", ""))
    for f in alignment_files
]

rule all:
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