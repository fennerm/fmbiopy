#!/usr/bin/env python
"""Extract random N bp subsequences from a fasta file."""
from __future__ import print_function
import sys

from Bio import SeqIO
import click
from numpy.random import randint
from pandas import read_csv


def select_rand_pos(interval_length, subinterval_length):
    return randint(0, interval_length - subinterval_length + 1)


def write_subseq(df_row, seq_record, subseq_length):
    while True:
        startpoint = select_rand_pos(df_row["length"], subseq_length)
        endpoint = startpoint + subseq_length
        subseq = str(seq_record.seq[startpoint:endpoint])
        if "N" not in subseq:
            break
    print(
        ">{}:{}-{}".format(df_row["chrom"], startpoint, endpoint),
        file=sys.stdout,
    )
    print(subseq, file=sys.stdout)


def extract_random_subseqs(fasta_filename, num_subseqs, subseq_length):
    faidx = read_csv(
        fasta_filename + ".fai",
        names=["chrom", "length"],
        sep="\t",
        usecols=[0, 1],
    )
    rand_csomes = faidx.sample(num_subseqs, weights="length", replace=True)
    rand_csomes["startpoint"] = rand_csomes["length"].apply(
        select_rand_pos, subinterval_length=subseq_length
    )
    for seq_record in SeqIO.parse(fasta_filename, format="fasta"):
        subseqs = rand_csomes[rand_csomes["chrom"] == seq_record.id]
        subseqs.apply(
            write_subseq,
            axis=1,
            seq_record=seq_record,
            subseq_length=subseq_length,
        )


@click.command()
@click.option(
    "-n", "--num-sequences", type=int, help="Number of subsequences to extract"
)
@click.option("-b", "--num-bases", type=int, help="Subsequence length")
@click.argument("fasta")
def cli(num_bases, num_sequences, fasta):
    """Extract random N bp subsequences from a fasta file."""
    extract_random_subseqs(fasta, num_sequences, num_bases)


if __name__ == "__main__":
    cli()
