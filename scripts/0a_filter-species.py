#!/usr/bin/env python3

from pathlib import Path  # Manipulating filenames
from Bio import SeqIO  # Reading DNA sequences
from Bio.SeqRecord import SeqRecord  # Create an empty sequence if the species is missing
from Bio.Seq import Seq  # Create an empty sequence if the species is missing
import sys # Retrieve output directory from Snakemake
import typer  # CLI argument handler

def main(
    alignment_file: Path = typer.Argument(
        ..., help="Path to the alignment file."
    ),
    wanted_species_file: Path = typer.Argument(
        ..., help="Path to the file containing species to filter down to."
    ),
    out_path: Path = typer.Argument(
        ..., help="Path to the output directory."
    )
):
    # Retrieve clade and alignment names for output
    try:
        clade_name = wanted_species_file.stem.split('_')[1]
    except IndexError:
        typer.echo("Error: The wanted_species_file is formatted incorrectly, it needs to be taxa_*.txt")
        raise typer.Exit(code=1)
    
    alignment_name = alignment_file.stem
    alignment_format = alignment_file.suffix.lstrip('.')
    out_filename = f"{alignment_name}_{clade_name}_filtered.{alignment_format}"
    output_path_filename = sys.argv[3]

    # Read in the wanted species names
    wanted_species_names = []
    with wanted_species_file.open('r') as species_file:
        for line in species_file:
            name = line.strip()
            if name:
                wanted_species_names.append(name)
    
    # Collect matching records
    matching_records = []
    sequence_length = None
    for record in SeqIO.parse(alignment_file, 'fasta'):
        if record.id in wanted_species_names:
            matching_records.append(record)
            wanted_species_names.remove(record.id)
            sequence_length = len(record.seq)  # Get sequence length for gap sequences
    
    if sequence_length is None:
        typer.echo("Error: No matching records found in the alignment file.")
        raise typer.Exit(code=1)
    
    # Add missing species with gap sequences
    for name in wanted_species_names:
        gap_seq = '-' * sequence_length
        empty_seq = SeqRecord(Seq(gap_seq), id=name, description='')
        matching_records.append(empty_seq)
    
    # Write the filtered records to a new file
    SeqIO.write(matching_records, output_path_filename, format=alignment_format)
    typer.echo(f"Filtered alignment written to {output_path_filename}")
    
    # Delete the old alignment file
    alignment_file.unlink()

if __name__ == "__main__":
    typer.run(main)