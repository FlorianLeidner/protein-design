import os
import re
import glob

import tempfile
import argparse
import subprocess

from copy import copy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices

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
        records = list(SeqIO.parse(fn, format="fasta"))
        sequence_records.extend(records[1:])
    return sequence_records

def parse_contigmap(contigs: str) -> tuple[int, int, int, int]:
    pattern = r"([0-9]+)-([0-9]+)"

    head = None
    designed = None
    tail = None

    contigs = contigs.lstrip("[").rstrip("]")

    elements = contigs.split(",")

    for element in elements:

        match = re.search(pattern, element)
        if match is None:
            raise ValueError(f"{element} did not match expected format.")

        start = int(match.group(1))
        stop = int(match.group(2))

        if element[0].isalpha():
            if head is None and designed is None:
                head = stop - start
            elif tail is None:
                tail = start - stop
        else:
            designed = (start, stop)

    if designed is None:
        print(contigs)
        print(head, designed, tail)
        raise RuntimeError("Did not find any variable elements. This is an indicator that something went wrong")

    return head, designed[0], designed[1], tail

def mask_seq(record: SeqRecord, contigs: str) -> SeqRecord:

    head, min_var, max_var, tail = parse_contigmap(contigs)

    seq = record.seq
    designed = seq[head:tail]

    if len(designed) > max_var or len(designed) < min_var:
        raise ValueError(f"Expected sequence length {min_var}-{max_var} but got {len(designed)}")

    record.seq = designed
    return record

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
    data_dict["seq_len"] = len(record.seq.replace("-", ""))

    return data_dict

def perform_msa(sequence_records: list[SeqRecord], muscle_exe: str ="muscle") -> list[SeqRecord]:
    """
    Perform Multiple Sequence Alignment (MSA) using MUSCLE on a list of amino acid sequences.

    Parameters
    ----------
    sequence_records : list of dict
        Each dict must contain:
          - 'name': unique identifier for the sequence
          - 'seq' : amino acid sequence string
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
        if len(sequence_records) < 1000:
            cmd = [muscle_exe, "-align", input_fasta, "-output", output_fasta]
        else:
            cmd = [muscle_exe, "-super5", input_fasta, "-output", output_fasta]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MUSCLE alignment failed:\n{e.stderr.decode()}")

    # Parse aligned sequences
    alignment = SeqIO.read(output_fasta, "fasta")

    os.remove(input_fasta)
    os.remove(output_fasta)

    return alignment

def pairwise_alignment(target: SeqRecord, query: SeqRecord) -> SeqRecord:

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -20
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(target.seq, query.seq)

    alignments_scored = []
    for alignment in alignments:
        alignments_scored.append((alignment, alignment.score))

    alignments_scored = sorted(alignments_scored, key=lambda x: x[1], reverse=True)

    alignment = alignments_scored[0][0]

    new_record = copy(query)
    new_record.seq = Seq(alignment[1])

    return new_record

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
                        nargs="+",
                        metavar="DIR",
                        help="One or more directories containing the output of one or more MPNN runs."
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

    parser.add_argument("--contigs",
                        dest = "contigs",
                        default = None,
                        metavar = "STRING",
                        help = "String defining the constant and variable parts of a design. Similar to rfdiffusion contigmaps. \
                         If this option is provided only the variable part of the design will be processed.")

    parser.add_argument("--write_fasta",
                        dest = "write_fasta",
                        default = False,
                        action = "store_true",
                        help = "If true write sequences to fasta")

    parser.add_argument("--align",
                        dest = "align",
                        default = False,
                        action = "store_true",
                        help = "Perform multisequence alignment")

    parser.add_argument("--target",
                        dest = "target",
                        default = None,
                        type = str,
                        metavar = "STRING",
                        help = "Align to target sequence. Can be sequence string or fasta file."
                        )

    parser.add_argument("--count",
                        dest = "count",
                        default= False,
                        action = "store_true",
                        help = "If true count amino occurrence. Requires --align if working with variable length designs")

    args = parser.parse_args()

    if args.indir is None and args.infiles is None:
        parser.error("Provide either input fasta files or a input directory")
    if args.indir is not None:
        for d in args.indir:
            if not os.path.isdir(d):
                parser.error(f"{d} doe not exist")
    if args.infiles is not None:
        for fn in args.infiles:
            if not os.path.isfile(fn):
                parser.error(f"{fn} does not exist")

    if args.target is not None:
        if os.path.isfile(args.target):
            try:
                args.target = SeqIO.read(args.target, format="fasta")
            except Exception as e:
                print(e)
                raise IOError("If query is file it must in in fasta format")
        else:
            seq = args.target.upper()
            seq = seq.strip()
            if not all([symbol.isalpha() for symbol in seq]):
                raise IOError(f"Sequence contains unknown characters: {seq}")
            else:
                args.target = SeqRecord(Seq(seq), id="query")

    return args

def main():
    args = parse_args()

    records = []
    if args.indir is not None:
        for d in args.indir:
            records.extend(records_from_dir(d))
    if args.infiles is not None:
        records.extend(records_from_files(args.infiles))

    target = args.target

    if args.contigs is not None:
        records = [mask_seq(record, args.contigs) for record in records]
        if target is not None:
            target = mask_seq(target, args.contigs)

    if args.align:
        if target is None:
            records = perform_msa(records)
        else:
            records = [pairwise_alignment(target, r) for r in records]

    data = pd.DataFrame([parse_seqrecord(record) for record in records])
    data.to_csv(f"{args.prefix}_output.csv")

    if args.write_fasta:
        fasta_file = f"{args.prefix}_sequences.fa"
        SeqIO.write(records, fasta_file, "fasta")

    # For batch size >1 Every run has one entry for run parameters [batchsize] entries with sequences
    sequences = data.loc[~data["id"].isna(), "seq"].values

    if args.count:
        counts = amino_acid_counts(sequences)
        counts_weighted = amino_acid_counts(sequences,
                                                    weights=data.loc[~data["id"].isna(), "overall_confidence"].values)
    
        counts.to_csv(f"{args.prefix}_aa_counts.csv")
        counts_weighted.to_csv(f"{args.prefix}_aa_counts_weighted.csv")

if __name__ == "__main__":
    main()
