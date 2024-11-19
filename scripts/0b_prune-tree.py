from pathlib import Path
from ete3 import Tree  # Manipulating trees
import typer

def main(
    tree_file: Path = typer.Argument(
        ..., help="Path to the tree file in Newick format."
    ),
    taxa_to_keep_file: Path = typer.Argument(
        ..., help="Path to the file containing taxa to keep."
    ),
    clade: str = typer.Argument(
        ..., help="Name of the clade being pruned down to."
    ),
):
    # Construct the output tree filename
    tree_name = tree_file.stem
    out_tree_filename = f"{tree_name}_{clade}.tre"
    out_tree_path = tree_file.parent / out_tree_filename

    # Read the tree
    try:
        full_tree = Tree(str(tree_file))
    except Exception as e:
        typer.echo(f"Error reading tree file '{tree_file}': {e}")
        raise typer.Exit(code=1)

    # Read the list of taxa to keep
    try:
        with taxa_to_keep_file.open('r') as taxa_list:
            taxa_to_keep = [line.strip() for line in taxa_list if line.strip()]
    except Exception as e:
        typer.echo(f"Error reading taxa file '{taxa_to_keep_file}': {e}")
        raise typer.Exit(code=1)

    # Check for taxa not present in the tree
    tree_leaf_names = set(leaf.name for leaf in full_tree.get_leaves())
    missing_taxa = [taxon for taxon in taxa_to_keep if taxon not in tree_leaf_names]
    if missing_taxa:
        typer.echo(
            f"Warning: The following taxa were not found in the tree and will be ignored:\n{', '.join(missing_taxa)}"
        )

    # Keep only taxa present in the tree
    taxa_to_keep_in_tree = [taxon for taxon in taxa_to_keep if taxon in tree_leaf_names]
    if not taxa_to_keep_in_tree:
        typer.echo("Error: None of the taxa to keep are present in the tree.")
        raise typer.Exit(code=1)

    # Prune the tree
    full_tree.prune(taxa_to_keep_in_tree, preserve_branch_length=True)

    # Ensure the tree is unrooted
    full_tree.unroot()

    # Write the pruned tree to a file
    try:
        full_tree.write(outfile=str(out_tree_path))
        typer.echo(f"Pruned tree written to {out_tree_path}")
    except Exception as e:
        typer.echo(f"Error writing pruned tree: {e}")
        raise typer.Exit(code=1)

if __name__ == "__main__":
    typer.run(main)