from pathlib import Path  # Manipulating filenames
from Bio import SeqIO  # Reading in sequences
import typer  # CLI argument handler

def main(
    alignment_file: Path = typer.Argument(
        ..., help="Path to the alignment file."
    )
):
    # Get the filename without extension
    alignment_name = alignment_file.stem
    # Get the file extension without the dot
    alignment_format = alignment_file.suffix[1:]
    # Output filename with '_trimmed' suffix
    alignment_outname = alignment_file.parent / f"{alignment_name}_trimmed.{alignment_format}"


    terminal_stops_trimmed_alignment = []
    for record in SeqIO.parse(alignment_file, format=alignment_format):
        # Trim the last three nucleotides (terminal stop codon)
        record.seq = record.seq[:-3]
        terminal_stops_trimmed_alignment.append(record)

    # Write the trimmed sequences to the output file
    SeqIO.write(terminal_stops_trimmed_alignment, alignment_outname, format=alignment_format)
    typer.echo(f"Trimmed alignment written to {alignment_outname}")

    # Delete the old alignment file
    alignment_file.unlink()

if __name__ == "__main__":
    typer.run(main)