#!/usr/bin/env python3

from pathlib import Path  # Manipulating filenames
import typer  # CLI argument handler
from Bio import SeqIO  # Reading sequences
from Bio.Data import CodonTable  # Check for oddities in the data

def check_frame(alignment_file: str):
    """
    Check if an alignment is in-frame, has no stop codons, or other erroneous characters.

    Args:
        alignment_file (str): Path to the alignment file to check.
    """

    # Use pathlib to handle file paths
    alignment_path = Path(alignment_file)
    
    # Retrieve file extension to determine the format, removing the leading dot
    alignment_format = alignment_path.suffix.lstrip('.')

    # Remove 'temp_' prefix from filename if present
    clean_alignment_file_name = alignment_path.name
    if clean_alignment_file_name.startswith("temp_"):
        clean_alignment_file_name = clean_alignment_file_name[len("temp_"):]

    # Check if the file exists
    if not alignment_path.is_file():
        typer.echo(f"File not found: {alignment_file}")
        raise typer.Exit(code=1)

    try:
        # Parse the first record to check sequence length
        record = next(SeqIO.parse(alignment_path, format=alignment_format))
        seq_len = len(record.seq)
    except Exception as e:
        typer.echo(f"Error reading {alignment_file}: {e}")
        raise typer.Exit(code=1)

    # Check if sequence length is a multiple of 3
    if seq_len % 3 != 0:
        typer.echo(f"{clean_alignment_file_name} is out of frame.")
    else:
        # Check each record for stop codons after translation
        for record in SeqIO.parse(alignment_path, format=alignment_format):
            record_N = record.seq.replace('-', 'N')
            try:
                aa_record = record_N.translate()
                if '*' in aa_record:
                    typer.echo(f"{clean_alignment_file_name} contains at least one stop codon in {record.id}.")
            except CodonTable.TranslationError as e:
                typer.echo(f"{clean_alignment_file_name} {record.id} has a Translation Error: {e}")

if __name__ == "__main__":
    typer.run(check_frame)
    