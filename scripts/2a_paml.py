#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from ete3 import EvolTree
from itertools import combinations
import typer

app = typer.Typer()

tests = {'branch': ['M0', 'b_free'],
         'bsA': ['bsA1', 'bsA'],
         'cmD': ['M3', 'bsD'],
         'cmC': ['XX', 'bsC']
         }  # NOTE XX refers to M2a_rel - this is just how to do a user-defined model

@app.command()
def main(
    alignment: str = typer.Argument(..., help="Path to the alignment file"),
    tree: str = typer.Argument(..., help="Path to the phylogenetic tree file"),
    test: str = typer.Argument(..., help="Test name to run (branch, bsA, cmD, cmC)"),
    test_taxa: str = typer.Argument(..., help="Path to the test_taxa file")
):
    alignment_file = alignment
    tree_file = tree
    test_models = tests[test]
    test_taxa_file = test_taxa

    alignment_name = Path(alignment_file).name
    gene_name = alignment_name.split('_')[0]
    clade_name = Path(alignment_file).stem.split('_')[1]
    alignment_format = Path(alignment_file).suffix[1:]  # Remove '.' from filetype


    # Check for empty sequences
    empty_seq_count = 0
    for record in SeqIO.parse(alignment_file, format=alignment_format):
        gapSeq = '-' * len(record.seq)
        if str(record.seq).upper().replace("N", "-") == gapSeq:
            empty_seq_count += 1

    # Process tree pruning if there are empty sequences
    if empty_seq_count >= 1:
        taxa_in_alignment = [
            record.id
            for record in SeqIO.parse(alignment_file, format=alignment_format)
            if str(record.seq).upper().replace("N", "-") != '-' * len(record.seq)
        ]

    tree = EvolTree(tree_file)
    out_tree_name = f"{Path(tree_file).stem}_{gene_name}.tre"

    if empty_seq_count >= 1 and len(taxa_in_alignment) >= 1:
        tree.prune(taxa_in_alignment, preserve_branch_length=True)
        tree.unroot()
        tree.write(outfile=out_tree_name, format=0)
        tree = EvolTree(out_tree_name)

    tree.link_to_alignment(alignment_file)
    tree.workdir = str(Path.cwd())

    # Collect node IDs for later use
    list_of_node_ids = [node.node_id for node in tree.traverse('postorder')]

    # Mark test taxa
    test_taxa = []
    with open(test_taxa_file, 'r') as test_taxa_list:
        test_taxa = [taxon.rstrip() for taxon in test_taxa_list]

    marked_taxon_ids = []
    for taxon in test_taxa:
        taxon_node = tree & taxon
        marked_taxon_id = taxon_node.node_id
        tree.mark_tree([marked_taxon_id])
        marked_taxon_ids.append(marked_taxon_id)

    # Find and mark internal nodes below test taxa
    for i in range(len(test_taxa), 1, -1):
        taxa_groups = combinations(test_taxa, i)
        for group in taxa_groups:
            common_node = tree.get_common_ancestor(*group)
            marked_taxon_id = common_node.node_id
            tree.mark_tree([marked_taxon_id])
            marked_taxon_ids.append(marked_taxon_id)

    # Best model and likelihood dictionaries
    best_model = {key: None for key in ['M0', 'b_free', 'bsA1', 'bsA', 'M3', 'bsD', 'XX', 'bsC']}
    best_lnL = {key: float('-inf') for key in ['M0', 'b_free', 'bsA1', 'bsA', 'M3', 'bsD', 'XX', 'bsC']}

    # Quicker version of running PAML for testing
    # for model in test_models:
    #     model_specifications = model + '.' + 'bl' + '_' + '0.7' + 'w'
    #     print(f"Testing model: {model} on: {alignment_name} with starting branch length option: bl and initial omega: 0.7w")
    #     tree.run_model(model_specifications, fix_blength=1, omega=0.7)
    #     current_model = tree.get_evol_model(model_specifications)
    #     print(f"""Model fitting of {alignment_name} complete, the likelihood was: {current_model.lnL} with these settings: 
    #           - Model: {model} 
    #           - Starting Branch Length Option: bl 
    #           - Initial Omega: 0.7w""")
    #     if current_model.lnL > best_lnL[model]:
    #         best_lnL[model] = current_model.lnL
    #         best_model[model] = current_model

    # Run each test
    for model in test_models:
        for starting_branch_length_option in [1, -1]:
            branch_estimation = 'bl' if starting_branch_length_option == 1 else 'random'
            for initial_omega in [0.2, 0.7, 1.2]:
                if model == 'bsA1':
                    initial_omega = 1.0
                model_specifications = f"{model}.{branch_estimation}_{initial_omega}w"
                print(f"Testing model: {model} on: {alignment_name} with starting branch length option: {branch_estimation} and initial omega: {initial_omega}w\n")
                
                if model == 'XX':
                    tree.run_model(model_specifications, fix_blength=starting_branch_length_option, omega=initial_omega, NSsites=22, ncatG=3)
                    tree.get_evol_model(model_specifications).properties['typ'] = 'branch-site'
                    tree.get_evol_model(model_specifications)._load(f"{model_specifications}/out")
                else:
                    tree.run_model(model_specifications, fix_blength=starting_branch_length_option, omega=initial_omega)
                
                current_model = tree.get_evol_model(model_specifications)
                print(f"""Model fitting of: {alignment_name} complete, the likelihood was: {current_model.lnL} with these settings:
                    - Model: {model}
                    - Starting Branch Length Option: {branch_estimation}
                    - Initial Omega: {initial_omega}w\n""")
                
                if current_model.lnL > best_lnL[model]:
                    best_lnL[model] = current_model.lnL
                    best_model[model] = current_model

                if model == 'bsA1':
                    break

    # Output results for all models
    for model in test_models:
        current_model = best_model[model]
        model_name = current_model.name
        lnL = current_model.lnL

        if model == 'M0':
            all_branch_stats = current_model.branches
            one_branch = all_branch_stats[1]
            omega = one_branch.get('w')
            results = f"{clade_name},{gene_name},{model_name},{lnL},{omega}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}\n")

        elif model == 'b_free':
            all_branch_stats = current_model.branches
            fg_branch = all_branch_stats[marked_taxon_id]
            fg_omega = fg_branch.get('w') if fg_branch.get('mark') == ' #1' else None

            bg_omega = None
            for node_id in list_of_node_ids:
                if node_id not in marked_taxon_ids:
                    bg_branch = current_model.branches[node_id]
                    if bg_branch.get('mark') == ' #0':
                        bg_omega = bg_branch.get('w')
                        break

            results = f"{clade_name},{gene_name},{model_name},{lnL},{bg_omega},{fg_omega}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}\n")

        elif model == 'M3':
            classes = current_model.classes
            proportions = classes.get('proportions')
            omegas = classes.get('w')
            results = f"{clade_name},{gene_name},{model_name},{lnL}," \
                    f"{proportions[0]},{omegas[0]}," \
                    f"{proportions[1]},{omegas[1]}," \
                    f"{proportions[2]},{omegas[2]}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}")

        elif model == 'bsD':
            classes = current_model.classes
            proportions = classes.get('proportions')
            background_omegas = classes.get('branch type 0')
            foreground_omegas = classes.get('branch type 1')
            results = f"{clade_name},{gene_name},{model_name},{lnL}," \
                    f"{proportions[0]},{background_omegas[0]},{foreground_omegas[0]}," \
                    f"{proportions[1]},{background_omegas[1]},{foreground_omegas[1]}," \
                    f"{proportions[2]},{background_omegas[2]},{foreground_omegas[2]}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}")

        elif model == 'XX':  # If model is M2a_rel
            classes = current_model.classes
            proportions = classes.get('proportions')
            omegas = classes.get('w')
            results = f"{clade_name},{gene_name},{model_name},{lnL}," \
                    f"{proportions[0]},{omegas[0]}," \
                    f"{proportions[1]},{omegas[1]}," \
                    f"{proportions[2]},{omegas[2]}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}")

        elif model == 'bsC':
            classes = current_model.classes
            proportions = classes.get('proportions')
            background_omegas = classes.get('branch type 0')
            foreground_omegas = classes.get('branch type 1')
            results = f"{clade_name},{gene_name},{model_name},{lnL}," \
                    f"{proportions[0]},{background_omegas[0]},{foreground_omegas[0]}," \
                    f"{proportions[1]},{background_omegas[1]},{foreground_omegas[1]}," \
                    f"{proportions[2]},{background_omegas[2]},{foreground_omegas[2]}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}")

        elif model == 'bsA':
            classes = current_model.classes
            proportions = classes.get('proportions')
            background_omegas = classes.get('background w')
            foreground_omegas = classes.get('foreground w')
            results = f"{clade_name},{gene_name},{model_name},{lnL}," \
                    f"{proportions[0]},{background_omegas[0]},{foreground_omegas[0]}," \
                    f"{proportions[1]},{background_omegas[1]},{foreground_omegas[1]}," \
                    f"{proportions[2]},{background_omegas[2]},{foreground_omegas[2]}," \
                    f"{proportions[3]},{background_omegas[3]},{foreground_omegas[3]}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}")

        elif model == 'bsA1':
            classes = current_model.classes
            proportions = classes.get('proportions')
            background_omegas = classes.get('background w')
            foreground_omegas = classes.get('foreground w')
            results = f"{clade_name},{gene_name},{model_name},{lnL}," \
                    f"{proportions[0]},{background_omegas[0]},{foreground_omegas[0]}," \
                    f"{proportions[1]},{background_omegas[1]},{foreground_omegas[1]}," \
                    f"{proportions[2]},{background_omegas[2]},{foreground_omegas[2]}," \
                    f"{proportions[3]},{background_omegas[3]},{foreground_omegas[3]}\n"
            out_filename = f"{clade_name}_{gene_name}_{model}.csv"
            with open(out_filename, 'w') as out_results:
                out_results.write(results)
            print(f"Results for the best model {model} written to {out_filename}")

if __name__ == "__main__":
    app()
