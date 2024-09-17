import json
from Bio import SeqIO
import argparse
import os
import logging


def is_nucleotide(sequence):
    nucleotides = set("ACTGN")
    nucleotide_ratio = sum(1 for char in sequence.upper() if char in nucleotides) / len(
        sequence
    )
    return nucleotide_ratio > 0.9


def check_sequence(sequence, name, is_dna):
    issues = []
    seq_length = len(sequence) if is_dna else len(sequence) * 3

    if seq_length < 1500:
        issues.append(
            f" less then 1500 nucleotides (or equivalent for proteins): {seq_length}"
        )

    if not is_dna and "X" in sequence:
        issues.append(f"Protein sequence contains 'X'")

    return issues


def create_json_data(name, sequence, is_dna):
    sequence_type = "dnaSequence" if is_dna else "proteinChain"
    return [
        {
            "name": name,
            "sequences": [{sequence_type: {"sequence": sequence, "count": 1}}],
        }
    ]


def process_sequence(record, output_dir):
    sequence = str(record.seq)
    name = record.id
    is_dna = is_nucleotide(sequence)
    file_name = name.replace("/", "_").replace("|", "_")
    output_filename = f"{output_dir}/{file_name}.json"

    issues = check_sequence(sequence, name, is_dna)
    if issues:
        if not os.path.exists(f"{output_dir}/bad"):
            os.makedirs(f"{output_dir}/bad")
        for issue in issues:
            logging.warning(f"ID: {name} - {issue}")

    output_filename = (
        f"{output_dir}/bad/{file_name}.json" if issues else output_filename
    )
    json_data = create_json_data(file_name, sequence, is_dna)

    with open(output_filename, "w") as json_file:
        json.dump(json_data, json_file, indent=2)


def fasta_to_json(input_file, output_dir):
    logging.info(f"Processing file: {input_file}")

    for record in SeqIO.parse(input_file, "fasta"):
        process_sequence(record, output_dir)

    logging.info("Processing complete.")


def main():
    parser = argparse.ArgumentParser(description="Convert FASTA file to JSON files.")
    parser.add_argument("-i", "--input", metavar="input_file", help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--outdir",
        metavar="output_dir",
        help="Prefix for output JSON files",
        default=".",
    )
    args = parser.parse_args()

    logging.basicConfig(
        filename=f"{args.outdir}_log.txt",
        level=logging.INFO,
        format="%(levelname)s - %(message)s",
        filemode="w",
    )
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    fasta_to_json(args.input, args.outdir)


if __name__ == "__main__":
    main()
