#!/usr/bin/env python3

#  Copyright (c) 2022. Dennis Wiersma.
#  Licensed under GPLv3. See LICENSE file.

"""
Module description
"""

# METADATA #

# IMPORTS #
import sys
import os
import argparse
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# CODE #
class MSAHandler:

    def __init__(self, path_to_protein_family):
        self.path_to_protein_family = path_to_protein_family

    def msa_from_fasta(self):
        """
        Use Clustalo on the command line to open a fasta file and align it.

        :return: MSA data in fasta format
        """
        try:
            msa_file = os.popen("/Users/denniswiersma/clustalo -i " + self.path_to_protein_family)
            return msa_file
        except ChildProcessError:
            print("Clustalo does not function correctly. Please check your install.")

    def read_msa(self):
        """
        Read the MSA file using BioPython's AlignIO object.

        :return: AlignIO object containing the MSA data.
        """
        try:
            return AlignIO.read(self.msa_from_fasta(), "fasta")
        except FileNotFoundError:
            print("Could not read MSA.")


class SNPHandler:
    def __init__(self, dna_seq_or_file):
        self.dna_seq_or_file = dna_seq_or_file

    def read_dna(self):
        """
        Takes either a DNA sequence directly or a fasta file and extracts it.

        :return: A SecRecord object.
        """

        # Check if sequence or filename was entered
        if (i.upper() in "ATCG" for i in self.dna_seq_or_file):

            # Prompt user for said parameters
            print("Looks like you've entered a DNA sequence. "
                  "You can optionally provide some more data. You may leave these empty if you wish.")
            seq_id = input("Please enter a sequence id: ")
            seq_name = input("Please enter a sequence name: ")
            seq_desc = input("Please enter a sequence description: ")

            # Return a SeqRecord instance
            return SeqRecord(
                Seq(self.dna_seq_or_file),
                id=seq_id,
                name=seq_name,
                description=seq_desc
            )

        else:
            try:
                # Open fasta file containing DNA sequence
                with open(self.dna_seq_or_file) as dna_file:
                    return SeqIO.parse(dna_file, "fasta")
            except FileNotFoundError:
                print("Could not read DNA sequence file.")


def get_arguments():
    """
    Creates arguments and handles interfacing with the cli using an ArgumentParser object.

    :return: parsed arguments from the cli.
    """
    parser = argparse.ArgumentParser(prog="Smoking Gun",
                                     description="Determine SNP severity for a given protein family.")
    parser.add_argument("ProteinFile",
                        type=str,
                        help="Path to the fasta file containing the protein family.")
    parser.add_argument("-p", "--snp",
                        required=True,
                        nargs=2,
                        help="Location of the SNP in the DNA sequence and the nucleotide to change to.\n"
                             "Example: 6 A\n"
                             "This would change position 6 to an A.")
    parser.add_argument("-s", "--sequence",
                        required=True,
                        type=str,
                        help="DNA sequence containing the SNP."
                             "Can be either the path to a FASTA file or a DNA sequence directly.")
    return parser.parse_args()


def main(args):
    """Function description"""
    arguments = get_arguments()

    msa_handler = MSAHandler(arguments.ProteinFile)
    aligned_msa = msa_handler.read_msa()

    snp_handler = SNPHandler(arguments.sequence)
    protein_dna = snp_handler.read_dna()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
