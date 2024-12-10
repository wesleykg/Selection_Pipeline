#!/bin/bash

parallel -u --joblog paml.log 'cd {} && python3 ../../../scripts/2a_paml.py *.fasta ../*.tre branch ../test_taxa_*.txt' ::: ../paml/*/*/
