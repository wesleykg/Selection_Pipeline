#!/bin/bash

find paml/*/*/*.csv -exec cat {} \; >> combined_results.csv

python scripts/3b_paml_stats.py combined_results.csv