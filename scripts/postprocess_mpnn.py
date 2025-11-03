import os
import glob

import tempfile
import argparse
import subprocess

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, AlignIO

import numpy as np
import pandas as pd


def records_from_dir(path=None) -> list[SeqRecord]:
    if path is None:
        path = os.getcwd()
    else:
        path = os.path.abspath(path)

    sequence_records = []

    for name in os.listdir(path=path):
        isdir = os.path.isdir(os.path.join(path, name))
        if isdir and name != "seqs":
            sequence_records.extend(records_from_dir(path=os.path.join(path, name)))
        elif isdir and name == "seqs":
            fastas = glob.glob(os.path.join(path, "seqs/*.fa"))
            records = records_from_files(fastas)
            sequence_records.extend(records)
        else:
            continue
    return sequence_records

def records_from_files(files: list) -> list[SeqRecord]:
    sequence_records = []
    for fn in files:
        sequence_records.extend([r for r in SeqIO.parse(fn, format="fasta")])
    return sequence_records

def parse_seqrecord(record: SeqRecord) -> dict:
    """
    Store the metadata from a mpnn sequence in a dictionary
    :param record:
    :return:
    """
    data_dict = {}

    data = record.description.split(",")

    data_dict["name"] = data.pop(0)
    for info in data:
        key, value = info.split("=")
        data_dict[key.strip()] = value.strip()

    data_dict["seq"] = str(record.seq)
    return data_dict

def perform_msa(sequence_records: list[SeqRecord], output_file: str ="alignment.afa", muscle_exe: str ="muscle") -> list[SeqRecord]:
    """
    Perform Multiple Sequence Alignment (MSA) using MUSCLE on a list of amino acid sequences.

    Parameters
    ----------
    sequence_records : list of dict
        Each dict must contain:
          - 'name': unique identifier for the sequence
          - 'seq' : amino acid sequence string
    output_file : str, optional
        File path where the aligned sequences will be saved in FASTA format.
    muscle_exe : str, optional
        Path to the MUSCLE executable (default assumes it's in PATH).

    Returns
    -------
    list of SeqRecords
        Updated list of sequence_records where 'seq' values are replaced by aligned sequences.
    """

    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fasta") as tmp_in:

        input_fasta = tmp_in.name
        SeqIO.write(sequence_records, tmp_in, "fasta")

    # --- Run MUSCLE ---
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fasta") as tmp_out:

        output_fasta = tmp_out.name

        cmd = [muscle_exe, "-align", input_fasta, "-output", output_fasta]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MUSCLE alignment failed:\n{e.stderr.decode()}")

    # Parse aligned sequences
    alignment = AlignIO.read(output_fasta, "fasta")

    # Write alignment to file
    AlignIO.write(alignment, output_file, "fasta")

    os.remove(input_fasta)
    os.remove(output_fasta)

    return alignment

def amino_acid_counts(sequences: list | np.ndarray, weights: list | np.ndarray = None) -> pd.DataFrame:
    """
    Main function that optionally performs MSA and computes amino acid frequencies per position.

    Parameters
    ----------
    sequences : list of str
        Amino acid sequences.
    weights : list of float or None

    Returns
    -------
    pandas.DataFrame
        Amino acid frequency table per sequence position.
    """

    if weights is not None:
        if len(weights) != len(sequences):
            raise ValueError(f"Weights ({len(weights)}) and sequences ({len(sequences)}) must be of the same length.")
    else:
        weights = np.ones(len(sequences))

    # Get sequences length and check if uniform
    seq_len = len(sequences[0])
    if not all([seq_len == len(seq) for seq in sequences]):
        raise ValueError("Can only count amino acids if all sequences are aligned.")

    columns = ["M", "A", "V", "I", "L", "W", "F", "Y", "T", "S", "C", "D", "E", "N", "Q", "H", "K", "R", "G", "P", "-"]

    column_dict = dict(list(zip(columns, np.arange(len(columns)))))

    counts = np.zeros((seq_len, len(columns)))

    for seq, w in zip(sequences, weights):
        for i, aa in enumerate(seq):
            counts[i, column_dict[aa]] += 1*float(w)

    # We do this to ensure that the final dict contains all symbols, even ones that are not observed
    df = pd.DataFrame(counts, columns=columns)

    df.index = range(1, len(df) + 1)
    df.index.name = "Position"

    return df

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="parse the results of a ligandmpnn run")

    parser.add_argument("-d",
                        dest="indir",
                        default=None,
                        type=str,
                        metavar="DIR",
                        help="A directory containing the output of one or more MPNN runs."
                             "The script will recursively identify all subdirectories named 'seqs' and parse their contents")

    parser.add_argument("-i",
                        dest="infiles",
                        default=None,
                        nargs="+",
                        metavar="FILE",
                        help="One or more input fasta files")

    parser.add_argument("--prefix",
                        dest = "prefix",
                        default = "mpnn_out",
                        metavar = "STRING",
                        help = "The output file")

    parser.add_argument("--align",
                        dest = "align",
                        default = False,
                        action="store_true",
                        help = "Perform multisequence alignment")

    args = parser.parse_args()

    if args.indir is None and args.infiles is None:
        parser.error("Provide either input fasta files or a input directory")
    elif args.indir is not None and not os.path.isdir(args.indir):
        parser.error(f"{args.indir} doe not exist")
    elif args.infiles is not None:
        for fn in args.infiles:
            if not os.path.isfile(fn):
                parser.error(f"{fn} does not exist")

    return args

def main():
    args = parse_args()

    records = []
    if args.indir is not None:
        records.extend(records_from_dir(args.indir))
    if args.infiles is not None:
        records.extend(records_from_files(args.infiles))

    if args.align:
        records = perform_msa(records, output_file = f"{args.prefix}.afa")

    data = pd.DataFrame([parse_seqrecord(record) for record in records])

    # For batch size >1 Every run has one entry for run parameters [batchsize] entrys with sequences
    sequences = data.loc[~data["id"].isna(), "seq"].values
    counts = amino_acid_counts(sequences)
    counts_weighted = amino_acid_counts(sequences,
                                                weights=data.loc[~data["id"].isna(), "overall_confidence"].values)

    counts.to_csv(f"{args.prefix}_aa_counts.csv")
    counts_weighted.to_csv(f"{args.prefix}_aa_counts_weighted.csv")
    data.to_csv(f"{args.prefix}_output.csv")

if __name__=="__main__":
    main()
