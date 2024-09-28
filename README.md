# Master-thesis-Wout-Vanderheyden

Spreadsheet:
https://docs.google.com/spreadsheets/d/1d01DQvNWWwdYsyP32cWNkdUFDZJ2UpF8rZ6Lsj_QUGY/edit?gid=0#gid=0

# Convert fasta to json

A small script that helped me convert the fasta files to a json file to submit to AF3.

```bash
python scripts/fasta_to_af3json.py -i data/<multi-fasta-sequence.fa> -o <outdir>
```

-   Stores json files of sequences that passed the checks:
    -   `<outdir>/<header>.json`
-   Stores json files of sequences that did not pass:
    -   `<outdir>/bad/<header>.json`
-   Creates a log file `<outdir>_log.txt` in which it stores the reason why certain sequences were excluded:

```txt
INFO - Processing file: gisaid_epiflu_sequence.fasta
WARNING - ID: A/Ishikawa/102/2002|EPI_ISL_11491 -  less then 1500 nucleotides (or equivalent for proteins): 987
WARNING - ID: A/Incheon/677/2006|EPI_ISL_20293 -  less then 1500 nucleotides (or equivalent for proteins): 1113
WARNING - ID: A/Johannesburg/15/2008|EPI_ISL_60841 -  less then 1500 nucleotides (or equivalent for proteins): 921
WARNING - ID: A/KITAKYUSHU/5/2006|EPI_ISL_17705 -  less then 1500 nucleotides (or equivalent for proteins): 987
WARNING - ID: A/Israel/26/2009|EPI_ISL_77133 -  less then 1500 nucleotides (or equivalent for proteins): 1101
WARNING - ID: A/Israel/24/2009|EPI_ISL_77132 -  less then 1500 nucleotides (or equivalent for proteins): 1101
WARNING - ID: A/IOWA/19/2010|EPI_ISL_88046 - Protein sequence contains 'X'
WARNING - ID: A/Israel/22/2009|EPI_ISL_77134 -  less then 1500 nucleotides (or equivalent for proteins): 1101
WARNING - ID: A/Kentucky/3/2006|EPI_ISL_21271 -  less then 1500 nucleotides (or equivalent for proteins): 987
WARNING - ID: A/Idaho/ID-2008/2003|EPI_ISL_17946 -  less then 1500 nucleotides (or equivalent for proteins): 987
INFO - Processing complete.
```

# Finding the codon of the "X" AA
Python script that helps finding the codon. 
Input files: protein_file.fa  nucleotide_file.fa

```bash
python scripts/find_missing_x.py 
```

```txt
Codon for 'X' at position 152: trt
Corresponding Amino Acid: X
```
