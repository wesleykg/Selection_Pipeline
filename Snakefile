import glob
import os

# Find all `.fasta` files in alignments/burmPartials/ and extract alignment names
alignment_files = glob.glob("alignments/burmPartials/*.fasta")
alignments = [os.path.basename(f).replace(".fasta", "") for f in alignment_files]

rule all:
    input:
        expand("alignments/burmPartials/{alignment}_burmPartials_filtered.fasta", alignment=alignments)

rule filter:
    input:
        "alignments/burmPartials/{alignment}.fasta",
        "alignments/burmPartials/taxa_burmPartials.txt"
    output:
        "alignments/burmPartials/{alignment}_burmPartials_filtered.fasta"
    shell:
        "scripts/0a_filter-species.py {input} {output}"