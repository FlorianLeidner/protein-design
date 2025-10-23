import argparse
import os
import glob

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import numpy as np
import pandas as pd


def records_from_dir(path=None) -> list:
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

def records_from_files(files: list) -> list:
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
        data_dict[key] = value

    data_dict["seq"] = str(record.seq)
    return data_dict

def perform_msa(sequences, msa_tool="clustalo"):
    raise NotImplementedError("Not implemented yet")


def calculate_position_frequencies(sequences: list | np.ndarray) -> pd.DataFrame:
    """
    Compute amino acid frequencies per position from aligned sequences.

    Parameters
    ----------
    sequences : list of str
        Aligned amino acid sequences (same length).

    Returns
    -------
    pandas.DataFrame
        Rows = sequence positions (1-indexed)
        Columns = amino acids
        Values = frequencies
    """

    max_len = max(len(seq) for seq in sequences)

    columns = ["M", "A", "V", "I", "L", "W", "F", "Y", "T", "S", "C", "D", "E", "N", "Q", "H", "K", "R", "G", "P", "-"]

    df = pd.DataFrame(columns=columns)

    aa_counts = []
    for pos in range(max_len):
        aa_list = [seq[pos] if pos < len(seq) else '-' for seq in sequences]
        count_dict = dict(list(zip(*np.unique(aa_list, return_counts=True))))
        aa_counts.append(count_dict)

    # We do this to ensure that the final dict contains all symbols, even ones that are not observed
    df = pd.concat((df, pd.DataFrame(aa_counts)))

    df.index = range(1, len(df) + 1)
    df.index.name = "Position"

    df = df.fillna(0.)

    df = df / len(sequences)

    return df


def amino_acid_frequencies(sequences, do_msa=False, msa_tool="clustalo"):
    """
    Main function that optionally performs MSA and computes amino acid frequencies per position.

    Parameters
    ----------
    sequences : list of str
        Amino acid sequences.
    do_msa : bool
        Whether to perform multiple sequence alignment first.
    msa_tool : str
        Which MSA tool to use if do_msa=True ('clustalo' or 'muscle').

    Returns
    -------
    pandas.DataFrame
        Amino acid frequency table per sequence position.
    """
    if do_msa:
        sequences = perform_msa(sequences, msa_tool=msa_tool)

    df = calculate_position_frequencies(sequences)
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
                        dest="prefix",
                        default="mpnn_out",
                        metavar="STRING",
                        help="The output file")

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
    if args.infile is not None:
        records.extend(records_from_files(args.infiles))

    data = pd.DataFrame([parse_seqrecord(record) for record in records])

    # For batch size >1 Every run has one entry for run parameters [batchsize] entrys with sequences
    sequences = data.loc[~data["id"].isna(), "seq"].values
    frequencies = calculate_position_frequencies(sequences)

    frequencies.to_csv(f"{args.prefix}_aa_frequencies.csv")
    data.to_csv(f"{args.prefix}_output.csv")

if __name__=="__main__":
    main()