from pathlib import Path
import glob

# Get all subdirectories inside alignments/
subdirs = [Path(d) for d in glob.glob("alignments/*") if Path(d).is_dir()]

# Each processing step deletes its input, so only one tier of fasta files
# exists at a time. Detect which stage we're at and derive base names accordingly.
# Stages: original → _filtered → _stops_trimmed → _gaps_trimmed
_orig_files = [
    f for subdir in subdirs for f in subdir.glob("*.fasta")
    if not any(kw in f.stem for kw in ("_filtered", "_stops_trimmed", "_gaps_trimmed"))
]
_filtered_files     = [f for subdir in subdirs for f in subdir.glob("*_filtered.fasta")]
_stops_trimmed_files = [f for subdir in subdirs for f in subdir.glob("*_stops_trimmed.fasta")]
_gaps_trimmed_files  = [f for subdir in subdirs for f in subdir.glob("*_gaps_trimmed.fasta")]

if _orig_files:
    alignments = [(f.parent.name, f.stem) for f in _orig_files]
elif _filtered_files:
    alignments = [(f.parent.name, f.stem[: -len("_filtered")]) for f in _filtered_files]
elif _stops_trimmed_files:
    alignments = [(f.parent.name, f.stem[: -len("_stops_trimmed")]) for f in _stops_trimmed_files]
elif _gaps_trimmed_files:
    alignments = [(f.parent.name, f.stem[: -len("_gaps_trimmed")]) for f in _gaps_trimmed_files]
else:
    alignments = []
filtered_subdirs = sorted({s[0] for s in alignments})

# Detect all tree files in `trees/original/`
original_tree_file = glob.glob("trees/original/*.tre")

# Detect all taxa list files in `lists/taxa/taxa_*.txt`
taxa_lists = glob.glob("lists/taxa/taxa_*.txt")
# Extract the lineage name by removing 'taxa_' and file extension
lineages = [Path(t).stem.replace("taxa_", "") for t in taxa_lists]


# ── Step 0: filter species & prune trees ──────────────────────────────

rule filter_all:
    input:
        expand(
            "alignments/{subdir}/{alignment}_filtered.fasta",
            zip,
            subdir=[s[0] for s in alignments],
            alignment=[s[1] for s in alignments]
        )

rule filter:
    input:
        "alignments/{subdir}/{alignment}.fasta",
        "lists/taxa/taxa_{subdir}.txt"
    output:
        "alignments/{subdir}/{alignment}_filtered.fasta"
    shell:
        "scripts/0a_filter-species.py {input} {output}"

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


# ── Step 1a: trim terminal stop codons ────────────────────────────────
# 1_trim-terminal-stops.py derives its own output name ({stem}_trimmed.fasta),
# so we rename it to _stops_trimmed.fasta to distinguish this step.

rule trim_stops_all:
    input:
        expand(
            "alignments/{subdir}/{alignment}_stops_trimmed.fasta",
            zip,
            subdir=[s[0] for s in alignments],
            alignment=[s[1] for s in alignments]
        )

rule trim_terminal_stops:
    input:
        "alignments/{subdir}/{alignment}_filtered.fasta"
    output:
        "alignments/{subdir}/{alignment}_stops_trimmed.fasta"
    shell:
        """
        python scripts/1_trim-terminal-stops.py {input}
        mv alignments/{wildcards.subdir}/{wildcards.alignment}_filtered_trimmed.fasta {output}
        """


# ── Step 1b: remove all-gap columns ───────────────────────────────────
# trimal writes to a fresh output file so no sentinel is needed.

rule trim_gaps_all:
    input:
        expand(
            "alignments/{subdir}/{alignment}_gaps_trimmed.fasta",
            zip,
            subdir=[s[0] for s in alignments],
            alignment=[s[1] for s in alignments]
        )

rule trim_gap_columns:
    input:
        "alignments/{subdir}/{alignment}_stops_trimmed.fasta"
    output:
        "alignments/{subdir}/{alignment}_gaps_trimmed.fasta"
    shell:
        """
        trimal -in {input} -out {output}.tmp -noallgaps -keepseqs
        sed -E 's/ [0-9]+ bp//g' {output}.tmp > {output}
        rm {output}.tmp {input}
        """


# ── Step 1c: check reading frame (diagnostic, independent branch) ─────
# Produces issues.txt per lineage; does not block the main pipeline.

rule check_frame_all:
    input:
        expand("alignments/{subdir}/issues.txt", subdir=filtered_subdirs)

rule check_frame:
    input:
        lambda wc: [
            f"alignments/{wc.subdir}/{stem}_gaps_trimmed.fasta"
            for subdir, stem in alignments
            if subdir == wc.subdir
        ]
    output:
        "alignments/{subdir}/issues.txt"
    shell:
        """
        rm -f {output}
        for fasta in alignments/{wildcards.subdir}/*_trimmed.fasta; do
            python scripts/1a_check-frame.py "$fasta" >> {output} 2>/dev/null || true
        done
        sort -u {output} -o {output} 2>/dev/null || true
        """


# ── Step 1d: concatenate alignments by gene group ─────────────────────
# 1b_alignment_concatenator.py creates gene-group subdirs inside the
# alignment dir, so outputs are dynamic — a sentinel records completion.

rule concat_all:
    input:
        expand("alignments/{subdir}/concat_done", subdir=filtered_subdirs)

rule concat_alignments:
    input:
        lambda wc: [
            f"alignments/{wc.subdir}/{stem}_gaps_trimmed.fasta"
            for subdir, stem in alignments
            if subdir == wc.subdir
        ]
    output:
        touch("alignments/{subdir}/concat_done")
    shell:
        """
        cd alignments/{wildcards.subdir}
        python ../../scripts/1b_alignment_concatenator.py
        """


# ── Step 1e: stage PAML working directories ───────────────────────────
# Declared as a checkpoint so Snakemake re-evaluates the downstream DAG
# after running it — the gene-group subdirs it creates aren't known until
# the rule actually executes.

rule paml_prep_all:
    input:
        expand("paml/{lineage}", lineage=lineages)

checkpoint paml_prep:
    input:
        concat_done="alignments/{lineage}/concat_done",
        taxa_file="lists/test_taxa/test_taxa_{lineage}.txt",
        tree_dir="trees/{lineage}"
    output:
        directory("paml/{lineage}")
    shell:
        """
        mkdir -p {output}

        for inner_dir in alignments/{wildcards.lineage}/*/; do
            if [ -d "$inner_dir" ] && [ -n "$(ls -A "$inner_dir" 2>/dev/null)" ]; then
                inner_name=$(basename "$inner_dir")
                mkdir -p {output}/$inner_name
                cp -r "$inner_dir/"* {output}/$inner_name/
            fi
        done

        cp {input.taxa_file} {output}/

        cp trees/{wildcards.lineage}/*_{wildcards.lineage}.tre {output}/ 2>/dev/null || true
        """


# ── Step 2: run PAML ──────────────────────────────────────────────────
# One job per gene-group directory so Snakemake can parallelise with -j N.
# paml_gene_group_outputs() uses the checkpoint to discover which
# gene-group subdirs were created for a given lineage.

def paml_gene_group_outputs(lineage):
    paml_dir = checkpoints.paml_prep.get(lineage=lineage).output[0]
    gene_groups = [d.name for d in Path(paml_dir).iterdir() if d.is_dir()]
    return expand(
        "paml/{lineage}/{gene_group}/paml_done",
        lineage=lineage,
        gene_group=gene_groups
    )

rule run_paml_all:
    input:
        lambda wc: [f for lineage in lineages for f in paml_gene_group_outputs(lineage)]

rule run_paml:
    input:
        fasta=lambda wc: glob.glob(f"paml/{wc.lineage}/{wc.gene_group}/*.fasta"),
        tree=lambda wc: glob.glob(f"paml/{wc.lineage}/*_{wc.lineage}.tre"),
        taxa=lambda wc: glob.glob(f"paml/{wc.lineage}/test_taxa_*.txt")
    output:
        touch("paml/{lineage}/{gene_group}/paml_done")
    shell:
        """
        cd paml/{wildcards.lineage}/{wildcards.gene_group}
        python3 ../../../scripts/2a_paml.py *.fasta ../*.tre branch ../test_taxa_*.txt
        """


# ── Step 3: collect and analyse results ───────────────────────────────

rule collect_results:
    input:
        lambda wc: [f for lineage in lineages for f in paml_gene_group_outputs(lineage)]
    output:
        combined="combined_results.csv",
        corrected="combined_results_corrected.csv"
    shell:
        """
        find paml/*/*/*.csv -exec cat {{}} \\; > {output.combined}
        python scripts/3b_paml_stats.py {output.combined}
        """