'''Usage: filter-species.py <alignment> <wanted_species>'''

# Modules
import os  # Manipulating filenames
from Bio import SeqIO  # Reading in DNA sequences
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Check if running interactively in an iPython console, or in a script from the
# command line
def in_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


# Run in a script from the command line
if in_ipython() is False:
    from docopt import docopt  # Command-line argument handler
    cmdln_args = docopt(__doc__)
    alignment_file = cmdln_args.get('<alignment>')
    wanted_species_file = cmdln_args.get('<wanted_species>')
# Run interactively in an iPython console
if in_ipython() is True:
    alignment_file = '4471.fasta'
    wanted_species_file = 'taxa_Dioscoreales.txt'

# Read in the wanted_species_file as a list of lines and retain the name in the
# list wanted_species_names.
wanted_species_names = []
with open(wanted_species_file, 'r') as species_file:
    wanted_species = species_file.readlines()
    for line in wanted_species:
        line = line.rstrip('\n')
        wanted_species_names.append(line)

# Loop through each wanted name and if it matches with a name in the alignment,
# add that sequence record to the list matching_records.
matching_records = []
for record in SeqIO.parse(alignment_file, 'fasta'):
    for name in wanted_species_names:
        if name == record.id:
            matching_records.append(record)
            wanted_species_names.remove(name)

# for record in matching_records:
#     for name in wanted_species_names:
#         if name not in record.id:
#             gap_seq = '-'*len(record)
#             matching_records.append(Seq(gap_seq, id=name))
#             wanted_species_names.remove(name)

missing_taxa = []
for name in wanted_species_names:
    gap_seq = '-'*len(record)
    empty_seq = SeqRecord(Seq(gap_seq), id=name, description='')
    matching_records.append(empty_seq)


# Write matching_records to file using the original filename and appending
# .fasta to the end.
alignment_name = os.path.splitext(alignment_file)[0]
SeqIO.write(matching_records, alignment_name + '.fasta', format='fasta')
